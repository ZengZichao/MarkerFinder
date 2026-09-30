"""Additional unit tests for marker selection helpers and module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from markerfinder.config import SelectionConfig
from markerfinder.exceptions import ExternalToolError
from markerfinder.models.genome import Genome, GeneState, OccupancyMatrix
from markerfinder.models.marker import SelectionStrategy
from markerfinder.models.pipeline_types import NoMarkerAvailableError
from markerfinder.modules.marker_selection import (
    AdaptiveMarkerSelectionModule,
    AdaptiveMarkerSelector,
    _parse_fasta,
    _run_hmmsearch,
    gather_marker_ids,
    hmm_path_for,
)


class TestGatherMarkerIds:
    def test_gathers_hmm_profiles(self, tmp_path):
        (tmp_path / "TIGR00001.hmm").write_text("")
        (tmp_path / "Pfam-A.hmm").write_text("")  # Merged lib, should be ignored
        ids = gather_marker_ids(str(tmp_path))
        assert ids == ["TIGR00001"]

    def test_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            gather_marker_ids(str(tmp_path / "missing"))

    def test_no_profiles(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            gather_marker_ids(str(tmp_path))


class TestHmmPathFor:
    def test_find_direct(self, tmp_path):
        (tmp_path / "TIGR00001.HMM").write_text("")
        assert hmm_path_for(str(tmp_path), "TIGR00001") == str(tmp_path / "TIGR00001.HMM")

    def test_find_in_individual_hmms(self, tmp_path):
        sub = tmp_path / "individual_hmms"
        sub.mkdir()
        (sub / "TIGR00001.hmm").write_text("")
        # No direct profile in tmp_path, so the individual_hmms subdirectory is used.
        found = hmm_path_for(str(tmp_path), "TIGR00001")
        assert "individual_hmms/TIGR00001" in found.replace("\\", "/")

    def test_not_found(self, tmp_path):
        assert hmm_path_for(str(tmp_path), "MISSING") is None


class TestRunHmmsearch:
    def _mock_run(self, hits):
        def runner(cmd, **kwargs):
            domtblout = cmd[cmd.index("--domtblout") + 1]
            lines = ["# header\n"]
            for h in hits:
                # 22 columns domtblout line
                line = (
                    f"{h['target']}\t-\t{h['target']}\t-\t0\t0\t{h.get('evalue', 1e-10)}\t{h['score']}\t"
                    f"0\t0\t{h.get('domain_evalue', 1e-10)}\t{h.get('domain_score', h['score'])}\t"
                    f"1\t100\t{h.get('ali_from', 1)}\t{h.get('ali_to', 100)}\t"
                    f"{h.get('env_from', 1)}\t{h.get('env_to', 100)}\t0\t0\t0\t-\n"
                )
                lines.append(line)
            Path(domtblout).resolve().write_text("".join(lines), encoding="utf-8")
            return MagicMock(returncode=0)
        return runner

    def test_returns_hits(self, tmp_path):
        hmm = tmp_path / "m.hmm"
        hmm.write_text("")
        fasta = tmp_path / "g.faa"
        fasta.write_text(">p1\nACDEFG\n")
        hits = [{"target": "p1", "score": 100.0, "ali_from": 1, "ali_to": 6}]
        with patch("subprocess.run", side_effect=self._mock_run(hits)):
            result = _run_hmmsearch(str(fasta), str(hmm), cpus=1, tmp_dir=str(tmp_path))
        assert "p1" in result
        assert result["p1"][0]["score"] == 100.0

    def test_missing_hmm_returns_empty(self, tmp_path):
        fasta = tmp_path / "g.faa"
        fasta.write_text(">p1\nACGT\n")
        assert _run_hmmsearch(str(fasta), str(tmp_path / "missing.hmm")) == {}

    def test_hmmsearch_not_found(self, tmp_path):
        # Hmmsearch 缺失必须硬失败, 不能再静默返回 {} 把"工具挂了"伪装成
        # "该 marker 在所有基因组中 ABSENT".
        hmm = tmp_path / "m.hmm"
        hmm.write_text("")
        fasta = tmp_path / "g.faa"
        fasta.write_text(">p1\nACGT\n")
        with patch("subprocess.run", side_effect=FileNotFoundError):
            with pytest.raises(ExternalToolError):
                _run_hmmsearch(str(fasta), str(hmm), tmp_dir=str(tmp_path))


class TestParseFasta:
    def test_parse_simple(self, tmp_path):
        fa = tmp_path / "s.faa"
        fa.write_text(">p1 first\nACGT\n>p2\nTGCA\n")
        seqs = _parse_fasta(str(fa))
        assert seqs == {"p1": "ACGT", "p2": "TGCA"}

    def test_missing_file(self, tmp_path):
        seqs = _parse_fasta(str(tmp_path / "missing.faa"))
        assert seqs == {}


class TestAdaptiveMarkerSelector:
    @pytest.fixture
    def selector(self):
        config = SelectionConfig(
            marker_mode="hmm",
            marker_hmm_dir="/tmp/fake_hmms",
            min_occupancy=0.5,
            max_markers=10,
        )
        mod = AdaptiveMarkerSelector(config, marker_ids=["M1", "M2"])
        return mod

    def test_scan_all_candidates_mocked(self, selector, tmp_path):
        from markerfinder.models.genome import Genome
        genomes = [Genome(id="G1", fasta_path=str(tmp_path / "G1.faa"))]
        (tmp_path / "G1.faa").write_text(">p1\nACGT\n")
        (tmp_path / "M1.hmm").write_text("")
        selector.config.marker_hmm_dir = str(tmp_path)

        def runner(cmd, **kwargs):
            domtblout = cmd[cmd.index("--domtblout") + 1]
            Path(domtblout).resolve().write_text(
                "p1\t-\tp1\t-\t0\t0\t1e-10\t50\t0\t0\t1e-10\t50\t"
                "1\t100\t1\t4\t1\t4\t0\t0\t0\t-\n",
                encoding="utf-8",
            )
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            matrix = selector.scan_all_candidates(genomes)
        assert matrix.get("G1", "M1") == GeneState.SINGLE_COPY

    def test_scan_target(self, selector, tmp_path):
        (tmp_path / "M1.hmm").write_text("")
        selector.config.marker_hmm_dir = str(tmp_path)

        def runner(cmd, **kwargs):
            domtblout = cmd[cmd.index("--domtblout") + 1]
            Path(domtblout).resolve().write_text(
                "p1\t-\tp1\t-\t0\t0\t1e-10\t50\t0\t0\t1e-10\t50\t"
                "1\t100\t1\t4\t1\t4\t0\t0\t0\t-\n",
                encoding="utf-8",
            )
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            hits = selector.scan_target(str(tmp_path / "g.faa"), "M1")
        assert len(hits) == 1

    def test_select_optimal_set_raises_for_unknown_strategy(self, selector):
        selector.config.strategy = "not_a_strategy"  # Type: ignore
        matrix = OccupancyMatrix(genomes=["G1"], cogs=["M1"])
        with pytest.raises(ValueError):
            selector.select_optimal_marker_set(matrix, {"M1": 1.0})

    def test_information_maximization(self, selector):
        selector.config.strategy = SelectionStrategy.INFO_MAX
        matrix = OccupancyMatrix(genomes=["G1", "G2"], cogs=["M1", "M2"])
        matrix.set("G1", "M1", GeneState.SINGLE_COPY)
        matrix.set("G2", "M1", GeneState.SINGLE_COPY)
        matrix.set("G1", "M2", GeneState.SINGLE_COPY)
        matrix.set("G2", "M2", GeneState.ABSENT)
        selected = selector.select_optimal_marker_set(matrix, {"M1": 1.0, "M2": 0.5})
        assert "M1" in selected.markers

    def test_rate_balanced_selection(self, selector):
        selector.config.strategy = SelectionStrategy.RATE_BALANCED
        matrix = OccupancyMatrix(genomes=["G1", "G2"], cogs=["M1"])
        matrix.set("G1", "M1", GeneState.SINGLE_COPY)
        matrix.set("G2", "M1", GeneState.MULTI_COPY)
        selected = selector.select_optimal_marker_set(matrix, {"M1": 0.5})
        assert "M1" in selected.markers

    def test_sparse_optimized_selection(self, selector):
        selector.config.strategy = SelectionStrategy.SPARSE_OPTIMIZED
        matrix = OccupancyMatrix(genomes=["G1", "G2"], cogs=["M1"])
        matrix.set("G1", "M1", GeneState.SINGLE_COPY)
        matrix.set("G2", "M1", GeneState.ABSENT)
        selected = selector.select_optimal_marker_set(matrix, {"M1": 0.5})
        assert "M1" in selected.markers


class TestAdaptiveMarkerSelectionModule:
    def test_module_discovers_marker_ids(self, tmp_path):
        (tmp_path / "TIGR00001.hmm").write_text("")
        config = SelectionConfig(marker_mode="hmm", marker_hmm_dir=str(tmp_path))
        mod = AdaptiveMarkerSelectionModule(config)
        assert mod.selector.marker_ids == ["TIGR00001"]

    def test_run_with_adaptive_params(self, tmp_path):
        from markerfinder.models.pipeline_types import AdaptiveParams
        from markerfinder.models.genome import Genome
        config = SelectionConfig(
            marker_mode="hmm", marker_hmm_dir=str(tmp_path),
            min_occupancy=0.75, max_markers=10,
        )
        (tmp_path / "M1.hmm").write_text("")
        (tmp_path / "G1.faa").write_text(">p1\nACGT\n")
        mod = AdaptiveMarkerSelectionModule(config)
        mod.selector.marker_ids = ["M1"]

        def runner(cmd, **kwargs):
            domtblout = cmd[cmd.index("--domtblout") + 1]
            Path(domtblout).resolve().write_text(
                "p1\t-\tp1\t-\t0\t0\t1e-10\t50\t0\t0\t1e-10\t50\t"
                "1\t100\t1\t4\t1\t4\t0\t0\t0\t-\n",
                encoding="utf-8",
            )
            return MagicMock(returncode=0)

        ctx = type("Ctx", (), {"adaptive_params": AdaptiveParams(min_occupancy=0.1, max_markers=5, min_hmm_score=10)})()
        with patch("subprocess.run", side_effect=runner):
            result = mod.run([Genome(id="G1", fasta_path=str(tmp_path / "G1.faa"))], phase_context=ctx)
        assert result.marker_set is not None
        assert "M1" in result.marker_set.markers

    def test_extract_marker_sequences(self, tmp_path):
        from markerfinder.models.genome import Genome
        config = SelectionConfig(marker_mode="hmm", marker_hmm_dir=str(tmp_path))
        (tmp_path / "M1.hmm").write_text("")
        (tmp_path / "G1.faa").write_text(">p1\nACDEFGHIKLMNPQRSTVWY\n>p2\nACDEFGHIKLMNPQRSTVWY\n")
        mod = AdaptiveMarkerSelectionModule(config)
        mod.selector.marker_ids = ["M1"]
        mod._hmm_cache = {
            "G1": {
                "M1": {
                    "p1": [{"target": "p1", "score": 100, "evalue": 1e-10, "length": 20, "ali_start": 1, "ali_end": 20}],
                    "p2": [{"target": "p2", "score": 50, "evalue": 1e-5, "length": 20, "ali_start": 1, "ali_end": 20}],
                }
            }
        }
        ms = type("MS", (), {"markers": ["M1"], "occupancy_scores": {"M1": 1.0}})()
        matrix = OccupancyMatrix(genomes=["G1"], cogs=["M1"])
        matrix.set("G1", "M1", GeneState.SINGLE_COPY)
        result = mod.extract_marker_sequences(
            [Genome(id="G1", fasta_path=str(tmp_path / "G1.faa"))],
            ms, matrix, tmp_dir=str(tmp_path)
        )
        assert len(result["M1"]) == 1
        assert result["M1"][0]["id"] == "G1"
