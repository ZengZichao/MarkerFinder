import pytest
from pathlib import Path

from markerfinder.config import ReportConfig
from markerfinder.models.genome import OccupancyMatrix
from markerfinder.models.marker import SelectedMarkerSet, MarkerQualityScore, MarkerLevel
from markerfinder.models.pipeline_types import (
    HGTReport, HGTEvaluation, MarkerSelectionResult, PhylogeneticResult,
    ConflictReport,
)
from markerfinder.models.report import RuntimeInfo, ReportOutput, PlainTextReportOutput
from markerfinder.modules.report_generator import (
    InteractiveReportGenerator, PlainTextReportGenerator, ReportGeneratorModule,
)


def make_marker_result():
    ms = SelectedMarkerSet(markers=["C1", "C2"], occupancy_scores={"C1": 0.9, "C2": 0.7})
    return MarkerSelectionResult(
        marker_set=ms,
        quality_scores={"C1": MarkerQualityScore(marker_id="C1", hmm_score=80, occupancy=0.9)},
        occupancy_matrix=OccupancyMatrix(genomes=["g1"], cogs=["C1", "C2"]),
    )


def make_hgt_report():
    return HGTReport(
        total_markers=2, level1_count=1, level2_count=1,
        marker_evaluations=[
            HGTEvaluation(marker_id="C1", overall_risk=0.1, level=MarkerLevel.LEVEL_1),
            HGTEvaluation(marker_id="C2", overall_risk=0.4, level=MarkerLevel.LEVEL_2),
        ],
    )


class TestInteractiveReportGenerator:
    def test_creates_html(self, tmp_path):
        cfg = ReportConfig(output_dir=str(tmp_path), output_prefix=str(tmp_path / "t"))
        gen = InteractiveReportGenerator(cfg)
        out = gen.generate(make_marker_result(), make_hgt_report(), PhylogeneticResult(), None, RuntimeInfo(duration=5.0))
        assert out.html_path.endswith(".html")
        assert Path(out.html_path).exists()

    def test_rf_distance_uses_conflict_report(self, tmp_path):
        cfg = ReportConfig(output_dir=str(tmp_path), output_prefix=str(tmp_path / "t"))
        gen = InteractiveReportGenerator(cfg)
        phylo = PhylogeneticResult(
            conflict_report=ConflictReport(
                rf_distance=2,
                normalized_rf=0.3333,
                n_quartet_conflicts=0,
                quartet_conflicts=[],
                n_conflicting_branches=0,
                conflicting_branches=[],
                gene_tree_agreement=None,
                conflict_summary={},
            )
        )
        out = gen.generate(make_marker_result(), make_hgt_report(), phylo, None, RuntimeInfo(duration=5.0))
        html = Path(out.html_path).read_text(encoding="utf-8")
        assert "0.3333" in html
        assert "Normalized RF distance" in html


class TestPlainTextReportGenerator:
    def test_creates_files(self, tmp_path):
        cfg = ReportConfig(output_dir=str(tmp_path), output_prefix=str(tmp_path / "t"))
        gen = PlainTextReportGenerator(cfg)
        out = gen.generate(make_marker_result(), make_hgt_report(), PhylogeneticResult(), None, RuntimeInfo(duration=5.0))
        assert len(out.file_paths) > 0
        for p in out.file_paths:
            assert Path(p).exists()


class TestReportGeneratorModule:
    def test_run(self, tmp_path):
        cfg = ReportConfig(output_dir=str(tmp_path), output_prefix=str(tmp_path / "t"))
        mod = ReportGeneratorModule(cfg)
        h, t = mod.run(make_marker_result(), make_hgt_report(), PhylogeneticResult(), None, RuntimeInfo(duration=5.0))
        assert isinstance(h, ReportOutput)
        assert isinstance(t, PlainTextReportOutput)
