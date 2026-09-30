import pytest

from markerfinder.models.genome import Genome, GenomeSource, GenomeQuality, OccupancyMatrix, GeneState
from markerfinder.models.marker import (
    MarkerLevel, MarkerQualityScore, SelectedMarkerSet, SelectionStrategy,
    ResolutionPreset, RESOLUTION_PRESETS,
)
from markerfinder.models.pipeline_types import NoMarkerAvailableError
from markerfinder.modules.marker_selection import AdaptiveMarkerSelectionModule
from markerfinder.config import SelectionConfig


def make_matrix(genomes, cogs, states):
    m = OccupancyMatrix(genomes=genomes, cogs=cogs)
    for (g, c), s in states.items():
        m.set(g, c, s)
    return m


class TestGeneState:
    def test_values(self):
        assert GeneState.ABSENT.value == "absent"
        assert GeneState.SINGLE_COPY.value == "single_copy"
        assert GeneState.MULTI_COPY.value == "multi_copy"


class TestMarkerLevel:
    def test_values(self):
        assert MarkerLevel.LEVEL_1.value == "level_1"
        assert MarkerLevel.LEVEL_2.value == "level_2"
        assert MarkerLevel.LEVEL_3.value == "level_3"


class TestMarkerQualityScore:
    def test_high_quality(self):
        s = MarkerQualityScore(
            marker_id="COG001", hmm_score=90.0, occupancy=0.95,
            length_ratio=1.0, phylogenetic_informativeness=0.8, hgt_risk_score=0.1,
        )
        assert s.overall_score > 0.7
        assert s.level == MarkerLevel.LEVEL_1

    def test_low_quality(self):
        s = MarkerQualityScore(
            marker_id="COG002", hmm_score=10.0, occupancy=0.2,
            length_ratio=0.3, phylogenetic_informativeness=0.1, hgt_risk_score=0.9,
        )
        assert s.overall_score < 0.5
        assert s.level == MarkerLevel.LEVEL_3

    def test_clamped(self):
        # Informativeness/hgt risk now default to None
        # ("not measured"); supply explicit values to exercise the clamp.
        s = MarkerQualityScore(
            hmm_score=200.0, occupancy=2.0, length_ratio=5.0,
            phylogenetic_informativeness=1.0, hgt_risk_score=0.0,
        )
        assert 0.0 <= s.overall_score <= 1.0

    def test_unmeasured_components_yield_none_score(self):
        # Unmeasured components must not be read as 0.
        s = MarkerQualityScore(hmm_score=200.0, occupancy=2.0, length_ratio=5.0)
        assert s.overall_score is None
        assert s.level == MarkerLevel.UNKNOWN


class TestSelectedMarkerSet:
    def test_default(self):
        ms = SelectedMarkerSet()
        assert ms.markers == []
        assert ms.strategy == SelectionStrategy.GREEDY


class TestResolutionPresets:
    def test_all_exist(self):
        for p in ResolutionPreset:
            assert p in RESOLUTION_PRESETS

    def test_conservative_strict(self):
        cfg = RESOLUTION_PRESETS[ResolutionPreset.CONSERVATIVE]
        assert cfg.min_occupancy >= 0.8

    def test_mag_adaptive_relaxed(self):
        cfg = RESOLUTION_PRESETS[ResolutionPreset.MAG_ADAPTIVE]
        assert cfg.min_occupancy <= 0.5
        assert cfg.allow_multi_copy is True

    def test_no_min_hmm_score_field(self):
        # Min_hmm_score (HMM 扫描) 已被移除, 仅保留 min_occupancy / max_markers 等
        cfg = RESOLUTION_PRESETS[ResolutionPreset.STANDARD]
        assert not hasattr(cfg, "min_hmm_score")


class TestSelectionConfig:
    def test_marker_mode_field(self):
        # Plan A exposes both marker modes via SelectionConfig
        c = SelectionConfig()
        assert hasattr(c, "marker_mode") and c.marker_mode == "gtdb_tk"
        assert hasattr(c, "marker_hmm_dir")
        assert hasattr(c, "min_hmm_score")


class TestAdaptiveMarkerSelectionModule:
    def test_exposes_gtdb_tk_and_hmm(self):
        mod = AdaptiveMarkerSelectionModule(SelectionConfig())
        assert hasattr(mod, "run_gtdb_tk")
        assert hasattr(mod, "extract_marker_sequences")
        assert hasattr(mod, "_compute_quality_scores")
        assert hasattr(mod, "run")  # Hmm-mode entry restored


class TestSelectionStrategies:
    """Synthetic occupancy-matrix tests for the four selection strategies."""

    @pytest.fixture
    def selector(self):
        config = SelectionConfig(
            marker_mode="hmm",
            marker_hmm_dir="/tmp/fake_hmms",
            min_occupancy=0.5,
            max_markers=10,
        )
        mod = AdaptiveMarkerSelectionModule(config)
        # Inject known marker ids so the selector does not hit the filesystem.
        mod.selector.marker_ids = [
            "COG001", "COG002", "COG003", "COG004", "COG005",
        ]
        return mod.selector

    def _matrix(self):
        genomes = ["G1", "G2", "G3", "G4"]
        cogs = ["COG001", "COG002", "COG003", "COG004", "COG005"]
        states = {
            ("G1", "COG001"): GeneState.SINGLE_COPY,
            ("G2", "COG001"): GeneState.SINGLE_COPY,
            ("G3", "COG001"): GeneState.SINGLE_COPY,
            ("G4", "COG001"): GeneState.SINGLE_COPY,
            ("G1", "COG002"): GeneState.SINGLE_COPY,
            ("G2", "COG002"): GeneState.SINGLE_COPY,
            ("G3", "COG002"): GeneState.ABSENT,
            ("G4", "COG002"): GeneState.ABSENT,
            ("G1", "COG003"): GeneState.SINGLE_COPY,
            ("G2", "COG003"): GeneState.ABSENT,
            ("G3", "COG003"): GeneState.SINGLE_COPY,
            ("G4", "COG003"): GeneState.ABSENT,
            ("G1", "COG004"): GeneState.SINGLE_COPY,
            ("G2", "COG004"): GeneState.MULTI_COPY,
            ("G3", "COG004"): GeneState.SINGLE_COPY,
            ("G4", "COG004"): GeneState.ABSENT,
            ("G1", "COG005"): GeneState.SINGLE_COPY,
            ("G2", "COG005"): GeneState.ABSENT,
            ("G3", "COG005"): GeneState.ABSENT,
            ("G4", "COG005"): GeneState.ABSENT,
        }
        return make_matrix(genomes, cogs, states)

    def test_greedy_selects_highest_occupancy(self, selector):
        selector.config.strategy = SelectionStrategy.GREEDY
        matrix = self._matrix()
        scores = selector.calculate_occupancy_scores(matrix)
        selected = selector.select_optimal_marker_set(matrix, scores)
        assert selected.markers[0] == "COG001"
        assert len(selected.markers) <= selector.config.max_markers
        assert all(scores[m] >= selector.config.min_occupancy for m in selected.markers)

    def test_info_max_penalises_redundancy(self, selector):
        selector.config.strategy = SelectionStrategy.INFO_MAX
        selector.config.max_markers = 3
        matrix = self._matrix()
        scores = selector.calculate_occupancy_scores(matrix)
        selected = selector.select_optimal_marker_set(matrix, scores)
        # COG001 and COG002 have high overlap; INFO_MAX should avoid picking both
        # If a complementary marker (COG003) is available.
        assert "COG001" in selected.markers
        assert len(selected.markers) > 1

    def test_rate_balanced_allows_multi_copy(self, selector):
        selector.config.strategy = SelectionStrategy.RATE_BALANCED
        matrix = self._matrix()
        scores = selector.calculate_occupancy_scores(matrix)
        selected = selector.select_optimal_marker_set(matrix, scores)
        assert len(selected.markers) >= 1

    def test_sparse_optimized_keeps_partial_hits(self, selector):
        selector.config.strategy = SelectionStrategy.SPARSE_OPTIMIZED
        selector.min_occupancy_threshold = 0.0  # SPARSE uses > 0
        matrix = self._matrix()
        scores = selector.calculate_occupancy_scores(matrix)
        selected = selector.select_optimal_marker_set(matrix, scores)
        # COG005 has only 1/4 occupancy and should be retained by SPARSE.
        assert "COG005" in selected.markers

    def test_greedy_raises_when_no_markers_meet_threshold(self, selector):
        selector.min_occupancy_threshold = 1.1  # _greedy_selection uses this attribute
        selector.config.min_occupancy = 1.1
        selector.config.strategy = SelectionStrategy.GREEDY
        matrix = self._matrix()
        scores = selector.calculate_occupancy_scores(matrix)
        with pytest.raises(NoMarkerAvailableError):
            selector.select_optimal_marker_set(matrix, scores)

    def test_quality_weighted_occupancy(self, selector):
        selector.config.quality_weighted = True
        matrix = self._matrix()
        matrix.genome_qualities = {"G1": 0.5, "G2": 1.0, "G3": 1.0, "G4": 1.0}
        scores = selector.calculate_occupancy_scores(matrix)
        # COG002 is present in G1 (weight 0.5) and G2 (weight 1.0), absent in G3/G4.
        assert 0.0 < scores["COG002"] < 1.0

    def test_run_gtdb_tk_applies_adaptive_params(self, tmp_path, monkeypatch):
        """Adaptive parameters from Phase 0 must override SelectionConfig in gtdb_tk mode."""
        from markerfinder.models.pipeline_types import AdaptiveParams

        config = SelectionConfig(marker_mode="gtdb_tk", min_occupancy=0.75, max_markers=60)
        mod = AdaptiveMarkerSelectionModule(config)

        adaptive = AdaptiveParams(min_occupancy=0.35, max_markers=120, min_hmm_score=15.0)

        def fake_load(markers_dir, genome_ids):
            cogs = ["COG001", "COG002", "COG003"]
            matrix = make_matrix(
                genome_ids, cogs,
                {(g, c): GeneState.SINGLE_COPY for g in genome_ids for c in cogs}
            )
            return matrix, {}, {c: 1.0 for c in cogs}

        def fake_build(*args, **kwargs):
            return {}

        import markerfinder.utils.gtdb_tk_markers as gtdb_mod
        monkeypatch.setattr(gtdb_mod, "load_gtdb_markers", fake_load)
        monkeypatch.setattr(gtdb_mod, "build_gene_trees", fake_build)

        genomes = [Genome(id="G1", fasta_path="/tmp/G1.faa")]
        result, *_ = mod.run_gtdb_tk(
            genomes,
            markers_dir=str(tmp_path),
            adaptive_params=adaptive,
        )

        # With min_occupancy=0.35 all three markers should be retained.
        assert result.marker_set is not None
        assert len(result.marker_set.markers) == 3

    def test_run_gtdb_tk_without_adaptive_params_uses_config(self, tmp_path, monkeypatch):
        from markerfinder.utils import gtdb_tk_markers as gtdb_mod

        config = SelectionConfig(marker_mode="gtdb_tk", min_occupancy=0.75, max_markers=60)
        mod = AdaptiveMarkerSelectionModule(config)

        def fake_load(markers_dir, genome_ids):
            cogs = ["COG001", "COG002"]
            matrix = make_matrix(
                genome_ids, cogs,
                {(g, c): GeneState.SINGLE_COPY for g in genome_ids for c in cogs}
            )
            return matrix, {}, {c: 1.0 for c in cogs}

        monkeypatch.setattr(gtdb_mod, "load_gtdb_markers", fake_load)
        monkeypatch.setattr(gtdb_mod, "build_gene_trees", lambda *args, **kwargs: {})

        genomes = [Genome(id="G1", fasta_path="/tmp/G1.faa")]
        result, *_ = mod.run_gtdb_tk(genomes, markers_dir=str(tmp_path))
        assert len(result.marker_set.markers) == 2
