"""Unit tests for HGT-aware marker filtering module."""

import pytest
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.tree import Tree
from markerfinder.models.pipeline_types import PhyloStepResult
from markerfinder.models.pipeline_types import HGTEvaluation, HGTReport
from markerfinder.config import HGTConfig
from markerfinder.modules.hgt_filter import (
    PhylogeneticHGTDetector,
    HGTDecisionEngine,
)


def _ete3_available() -> bool:
    """MAD rooting needs ete3; monophyly-screen assertions are
    environment-aware (documented NOT_MEASURABLE behaviour without it).
    Uses a real import probe, NOT find_spec — a broken install still has a
    spec. `markerfinder` is imported first because it is what installs the
    stdlib `cgi` stand-in that ete3 needs on Python 3.13+; probing ete3 in a
    bare namespace would report every supported interpreter as unusable."""
    try:
        import markerfinder  # Noqa: F401  -- installs the cgi shim when needed

        import ete3  # Noqa: F401
        return True
    except Exception:
        return False


@pytest.fixture
def default_hgt_config():
    return HGTConfig()


class TestPhylogeneticHGTDetector:
    def test_creation(self, default_hgt_config):
        detector = PhylogeneticHGTDetector(default_hgt_config)
        assert detector is not None

    def test_detect_with_trees(self, default_hgt_config):
        detector = PhylogeneticHGTDetector(default_hgt_config)
        gene_tree = Tree(newick="((A,B),(C,D));")
        species_tree = Tree(newick="((A,C),(B,D));")

        result = detector.detect("gene1", gene_tree, species_tree)
        assert isinstance(result, PhyloStepResult)
        assert 0 <= result.overall_risk <= 1


    def test_detect_monophyly_clean_marker(self, default_hgt_config):
        detector = PhylogeneticHGTDetector(default_hgt_config)
        # Gene tree matches taxonomy: B1/B2 monophyletic at genus level.
        gene_tree_newick = "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));"
        taxonomy_map = {
            "B1": {"domain": "Bacteria", "phylum": "P1", "genus": "B"},
            "B2": {"domain": "Bacteria", "phylum": "P1", "genus": "B"},
            "A1": {"domain": "Bacteria", "phylum": "P1", "genus": "A"},
            "A2": {"domain": "Bacteria", "phylum": "P1", "genus": "A"},
        }
        result = detector.detect_monophyly(
            "m1", gene_tree_newick, taxonomy_map, level="genus", threshold=0.5
        )
        if _ete3_available():
            assert result is not None
            assert result.monophyly_proportion == 1.0
            assert result.overall_risk == 0.0
            assert result.is_suspicious is False
        else:
            # MAD rooting still requires
            # Ete3; without it the monophyly screen is NOT_MEASURABLE and
            # Detect_monophyly returns None (→ UNKNOWN once lands).
            assert result is None

    def test_detect_monophyly_incongruent_marker(self, default_hgt_config):
        detector = PhylogeneticHGTDetector(default_hgt_config)
        # Gene tree places one A with B, breaking genus monophyly.
        gene_tree_newick = "((B1:0.1,A1:0.2),(B2:0.3,A2:0.4));"
        taxonomy_map = {
            "B1": {"domain": "Bacteria", "phylum": "P1", "genus": "B"},
            "B2": {"domain": "Bacteria", "phylum": "P1", "genus": "B"},
            "A1": {"domain": "Bacteria", "phylum": "P1", "genus": "A"},
            "A2": {"domain": "Bacteria", "phylum": "P1", "genus": "A"},
        }
        result = detector.detect_monophyly(
            "m1", gene_tree_newick, taxonomy_map, level="genus", threshold=0.5
        )
        if _ete3_available():
            assert result is not None
            assert result.monophyly_proportion == 0.0
            assert result.overall_risk == 1.0
            assert result.is_suspicious is True
        else:
            # MAD rooting unavailable without ete3 → None (NOT_MEASURABLE).
            assert result is None

    def test_quartet_consistency_perfect(self, default_hgt_config):
        detector = PhylogeneticHGTDetector(default_hgt_config)
        gene_tree = Tree(newick="((A,B),(C,D));")
        species_tree = Tree(newick="((A,B),(C,D));")
        q, reason = detector._calculate_quartet_consistency(gene_tree, species_tree)
        assert q == 1.0
        assert reason == ""

    def test_quartet_consistency_inconsistent(self, default_hgt_config):
        detector = PhylogeneticHGTDetector(default_hgt_config)
        gene_tree = Tree(newick="((A,C),(B,D));")
        species_tree = Tree(newick="((A,B),(C,D));")
        q, reason = detector._calculate_quartet_consistency(gene_tree, species_tree)
        assert q < 1.0
        assert reason == ""

    def test_far_distance_relaxation(self, default_hgt_config):
        from markerfinder.config import HGTConfig
        config = HGTConfig(adaptive_far_thresholds=True)
        config._far_distance_active = True
        engine = HGTDecisionEngine(config)
        result = engine.evaluate_marker(
            "COG003",
            phylogenetic_result=PhyloStepResult(overall_risk=0.7),
        )
        # With far-distance relaxation, level2_max is raised to 0.95,
        # So risk=0.7 stays in Level 2 rather than Level 3.
        assert result.level == MarkerLevel.LEVEL_2

    def test_no_phylogenetic_result_is_unknown(self, default_hgt_config):
        # 无参考物种树/无 taxonomy 表时系统发育筛查无法运行, 必须显式
        # 标记为 UNKNOWN, 不能再冒充 Level 3(高 HGT 风险并排除).
        engine = HGTDecisionEngine(default_hgt_config)
        result = engine.evaluate_marker("COG001")
        assert result.level == MarkerLevel.UNKNOWN
        assert result.confidence == "unknown"
        assert "UNKNOWN" in result.notes

    def test_phylogenetic_low_risk_level1(self, default_hgt_config):
        engine = HGTDecisionEngine(default_hgt_config)
        result = engine.evaluate_marker(
            "COG001",
            phylogenetic_result=PhyloStepResult(overall_risk=0.1),
        )
        assert result.level == MarkerLevel.LEVEL_1
        # Risk=0.1 sits 0.15 from the L1/L2 boundary
        # (relative distance 0.25) => medium, even though level is LEVEL_1.
        assert result.confidence == "medium"

    def test_phylogenetic_medium_risk_level2(self, default_hgt_config):
        engine = HGTDecisionEngine(default_hgt_config)
        result = engine.evaluate_marker(
            "COG003",
            phylogenetic_result=PhyloStepResult(overall_risk=0.4),
        )
        assert result.level == MarkerLevel.LEVEL_2
        # Risk=0.4 is mid-band: relative distance 0.33 => medium.
        assert result.confidence == "medium"

    def test_phylogenetic_high_risk_level3(self, default_hgt_config):
        engine = HGTDecisionEngine(default_hgt_config)
        result = engine.evaluate_marker(
            "COG002",
            phylogenetic_result=PhyloStepResult(overall_risk=0.9),
        )
        assert result.level == MarkerLevel.LEVEL_3
        # Risk=0.9 is far from both boundaries => high under semantics
        # (distance-based, no longer a synonym of level).
        assert result.confidence == "high"

    def test_generate_report(self, default_hgt_config):
        engine = HGTDecisionEngine(default_hgt_config)
        evaluations = [
            HGTEvaluation(marker_id="c1", overall_risk=0.1, level=MarkerLevel.LEVEL_1),
            HGTEvaluation(marker_id="c2", overall_risk=0.5, level=MarkerLevel.LEVEL_2),
            HGTEvaluation(marker_id="c3", overall_risk=0.8, level=MarkerLevel.LEVEL_3),
        ]
        report = engine.generate_hgt_report(evaluations)

        assert isinstance(report, HGTReport)
        assert report.total_markers == 3
        assert report.level1_count == 1
        assert report.level2_count == 1
        assert report.level3_count == 1
