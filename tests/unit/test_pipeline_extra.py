"""Additional unit tests for MarkerFinderPipeline orchestration."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from markerfinder.config import PipelineConfig, SelectionConfig
from markerfinder.models.genome import Genome
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.report import PhaseContext, RuntimeInfo
from markerfinder.models.pipeline_types import (
    AdaptiveParams,
    HGTReport,
    HGTEvaluation,
    MarkerSelectionResult,
    PhylogeneticResult,
    PreprocessingResult,
    QualityData,
    SupermatrixResult,
    CoalescentResult,
)
from markerfinder.models.tree import Tree
from markerfinder.pipeline import MarkerFinderPipeline, _gene_tree_mean_support


class TestPipelineHelpers:
    def test_gene_tree_mean_support_no_support(self):
        """G1/: a tree with no evaluable support is UNMEASURED, not 0.5.

        This assertion used to lock in ``0.5`` — the same neutral-placeholder
        defect as baseline re-introduced on the Phase 3 back-fill path.
        Changed (not deleted) per.
        """
        tree = Tree(newick="((A,B),(C,D));")
        assert _gene_tree_mean_support(tree) is None

    def test_gene_tree_mean_support_with_values(self):
        """The back-fill must actually measure, not abstain by accident.

        This used to read ``val is None or 0 <= val <= 1``, which stayed green
        while the helper returned ``None`` for *every* marker (it walked
        ``.traverse`` on a model object that has no such method). The contract
        is now stated per interpreter, so neither a dead feature nor a silent
        abstention can pass as a measurement.
        """
        from markerfinder.utils.etree import require_ete3

        tree = Tree(newick="((A:0.1,B:0.2)0.9,(C:0.1,D:0.2)0.8);")
        val = _gene_tree_mean_support(tree)
        try:
            require_ete3()
        except Exception:
            assert val is None, (
                "ete3 absent: support reading is NOT_MEASURABLE, and must never "
                "be a placeholder"
            )
        else:
            assert val == pytest.approx(0.85), val


class TestPipelineSteps:
    def test_run_scan_hmm_mode(self, tmp_path):
        cfg = PipelineConfig(
            input_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            tmp_dir=str(tmp_path / "tmp"),
            selection_config=SelectionConfig(marker_mode="hmm", marker_hmm_dir=str(tmp_path), min_occupancy=0.0, max_markers=5),
        )
        (tmp_path / "G1.faa").write_text(">p1\nACGT\n")
        (tmp_path / "M1.hmm").write_text("")
        pipeline = MarkerFinderPipeline(cfg)

        def runner(cmd, **kwargs):
            if cmd[0].lower() == "hmmsearch":
                domtblout = cmd[cmd.index("--domtblout") + 1]
                Path(domtblout).resolve().write_text(
                    "p1\t-\tp1\t-\t0\t0\t1e-10\t50\t0\t0\t1e-10\t50\t"
                    "1\t100\t1\t4\t1\t4\t0\t0\t0\t-\n",
                    encoding="utf-8",
                )
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            scan = pipeline.run_scan([Genome(id="G1", fasta_path=str(tmp_path / "G1.faa"))], PhaseContext())
        assert "marker_selection" in scan
        assert "marker_sequences" in scan

    def test_run_filter_backfills_quality(self, tmp_path):
        cfg = PipelineConfig(
            output_dir=str(tmp_path / "out"),
            tmp_dir=str(tmp_path / "tmp"),
        )
        pipeline = MarkerFinderPipeline(cfg)
        marker_selection = MagicMock()
        marker_selection.marker_set.markers = ["M1"]
        marker_selection.quality_scores = {"M1": MagicMock(hgt_risk_score=0.0)}
        marker_selection.occupancy_matrix = MagicMock()
        scan_data = {
            "marker_selection": marker_selection,
            "marker_sequences": {"M1": [{"id": "G1", "seq": "ACGT"}]},
            "rank_map": {"M1": 1},
        }
        hgt_report = HGTReport(
            total_markers=1,
            marker_evaluations=[HGTEvaluation(marker_id="M1", overall_risk=0.1, level=MarkerLevel.LEVEL_1)],
        )
        with patch.object(pipeline.hgt_filter, "run", return_value=(hgt_report, {"M1": MarkerLevel.LEVEL_1}, False)):
            result = pipeline.run_filter([Genome(id="G1", fasta_path="x")], scan_data, PhaseContext())
        assert result["hgt_report"] == hgt_report
        assert "M1" in result["marker_genes"]

    def test_run_report(self, tmp_path):
        cfg = PipelineConfig(output_dir=str(tmp_path / "out"), tmp_dir=str(tmp_path / "tmp"))
        pipeline = MarkerFinderPipeline(cfg)
        ms = MagicMock()
        ms.marker_set.markers = ["M1"]
        ms.marker_set.occupancy_scores = {"M1": 1.0}
        ms.marker_set.mean_occupancy = 1.0
        ms.quality_scores = {}
        ms.occupancy_matrix = MagicMock(n_genomes=1)
        scan_data = {"marker_selection": ms, "rank_map": {"M1": 1}}
        hgt_report = HGTReport(total_markers=1, marker_evaluations=[])
        phylo = PhylogeneticResult(
            supermatrix=SupermatrixResult(tree=Tree(newick="(G1);"), partition=MagicMock()),
            coalescent=CoalescentResult(species_tree=Tree(newick="(G1);"), species_tree_source="astral"),
        )
        filter_data = {"hgt_report": hgt_report, "marker_genes": {"M1": []}, "precomputed_levels": {}}
        infer_data = {"phylo_result": phylo}
        runtime = RuntimeInfo(duration=1.0)
        report = pipeline.run_report(scan_data, filter_data, infer_data, PhaseContext(), runtime)
        assert "html_output" in report
        assert "text_output" in report

    def test_write_run_config(self, tmp_path):
        cfg = PipelineConfig(output_dir=str(tmp_path / "out"))
        pipeline = MarkerFinderPipeline(cfg)
        pipeline._write_run_config(0.0, 1.0)
        run_config = Path(tmp_path / "out" / "Phase5_metadata" / "run_config.json")
        assert run_config.exists()
        import json
        data = json.loads(run_config.read_text())
        assert "parameters" in data

    def test_recorded_paths_are_resolved_not_relativized(self, tmp_path):
        """A snapshot must say where the run looked, in a form that can be read.

        Relative-to-output paths looked portable but were uninterpretable: a
        replay resolved them against the current working directory and landed in
        the previous run's output tree (or on a read-only ``/tmp``). MUST-FAIL
        CONTROL: the same input under the old rule produced ``../..``-shaped
        values, which is exactly what this pins out.
        """
        cfg = PipelineConfig(output_dir=str(tmp_path / "out"), input_dir="/abs/in")
        pipeline = MarkerFinderPipeline(cfg)
        resolved = pipeline._record_paths_resolved(
            {"input_dir": "/abs/in", "output_dir": str(tmp_path / "out"),
             "selection_config": {"gtdb_markers_dir": "relative/markers"}})
        assert Path(resolved["input_dir"]).is_absolute()
        assert resolved["input_dir"] == "/abs/in"
        # A relative value the user typed is recorded as what it meant then.
        assert Path(resolved["selection_config"]["gtdb_markers_dir"]).is_absolute()
        assert not resolved["selection_config"]["gtdb_markers_dir"].startswith("..")

    def test_save_intermediates(self, tmp_path):
        cfg = PipelineConfig(
            output_dir=str(tmp_path / "out"),
            tmp_dir=str(tmp_path / "tmp"),
            output_prefix="markerfinder",
            save_intermediates=True,
        )
        pipeline = MarkerFinderPipeline(cfg)
        tmp_dir = Path(tmp_path / "tmp")
        tmp_dir.mkdir()
        (tmp_dir / "M1.faa").write_text(">G1\nACGT\n")
        (tmp_dir / "M1.aln").write_text(">G1\nACGT\n")
        (tmp_dir / "markerfinder.concat.fasta").write_text(">G1\nACGT\n")
        (tmp_dir / "checkm.tsv").write_text("header\n")
        pipeline._save_intermediates()
        assert (tmp_path / "out" / "Phase4_intermediate").exists()


class TestPipelineStepSubcommands:
    def test_run_step_scan(self, tmp_path):
        cfg = PipelineConfig(
            input_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            tmp_dir=str(tmp_path / "tmp"),
            selection_config=SelectionConfig(marker_mode="hmm", marker_hmm_dir=str(tmp_path), min_occupancy=0.0, max_markers=5),
        )
        (tmp_path / "G1.faa").write_text(">p1\nACGT\n")
        (tmp_path / "M1.hmm").write_text("")
        pipeline = MarkerFinderPipeline(cfg)

        def runner(cmd, **kwargs):
            if cmd[0].lower() == "hmmsearch":
                domtblout = cmd[cmd.index("--domtblout") + 1]
                Path(domtblout).resolve().write_text(
                    "p1\t-\tp1\t-\t0\t0\t1e-10\t50\t0\t0\t1e-10\t50\t"
                    "1\t100\t1\t4\t1\t4\t0\t0\t0\t-\n",
                    encoding="utf-8",
                )
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = pipeline.run_step("scan", genomes=[Genome(id="G1", fasta_path=str(tmp_path / "G1.faa"))])
        assert "completed_steps" in result or "marker_selection" in result

    def test_run_step_missing_prerequisite(self, tmp_path):
        cfg = PipelineConfig(output_dir=str(tmp_path / "out"), tmp_dir=str(tmp_path / "tmp"))
        pipeline = MarkerFinderPipeline(cfg)
        with pytest.raises(RuntimeError):
            pipeline.run_step("filter")
