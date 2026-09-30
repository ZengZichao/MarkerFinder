""" Second bullet: the independence verdict must reach a product.

 measured how independent the concatenation and coalescent legs are and
stored it on ``ConflictReport.independence`` — and nothing ever read the field.
An independence caveat that lives only in memory is the same defect the
the review of the pipeline listed as baseline ("computed, stored, then dropped"),
so this test pins the rendering.
"""

from __future__ import annotations

from pathlib import Path

from markerfinder.config import ReportConfig
from markerfinder.models.marker import SelectedMarkerSet
from markerfinder.models.pipeline_types import (
    ConflictReport,
    HGTEvaluation,
    HGTReport,
    IndependenceReport,
    MarkerSelectionResult,
    PhylogeneticResult,
)
from markerfinder.modules.report_generator import PlainTextReportGenerator


class _Runtime:
    duration = 1.0


def _generate(tmp_path, independence):
    conflict = None
    if independence is not None:
        conflict = ConflictReport(independence=independence)
    phylo = PhylogeneticResult(conflict_report=conflict)

    marker = "M1"
    marker_set = SelectedMarkerSet(
        markers=[marker], occupancy_scores={marker: 0.9}, quality_scores={},
    )
    marker_result = MarkerSelectionResult(
        marker_set=marker_set, quality_scores={},
    )
    hgt_report = HGTReport(
        marker_evaluations=[HGTEvaluation(
            marker_id=marker, overall_risk=0.1,
        )],
        total_markers=1,
        evidence_coverage=1.0,
    )

    config = ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
    generator = PlainTextReportGenerator(config)
    out_dir = tmp_path / "Phase5_reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    generator._write_pipeline_summary(
        marker_result, hgt_report, phylo, _Runtime(), out_dir, "mf",
    )
    return (out_dir / "mf.pipeline_summary.txt").read_text(encoding="utf-8")


class TestIndependenceIsSurfaced:
    # The R2 independence record (Jaccard, shared gene-tree ratio,
    # Reference-built-from-tested-markers) must reach a product.
    def test_shared_legs_are_declared_in_the_summary(self, tmp_path):
        body = _generate(tmp_path, IndependenceReport(
            marker_set_jaccard=1.0,
            shared_gene_tree_ratio=0.83,
            same_trimming_regime=False,
            ref_built_from_tested_markers=True,
        ))
        assert "Leg independence:" in body
        assert "shared gene trees 83%" in body
        # The conclusion a reader must not miss: the legs are not independent.
        assert "NOT independent upstream" in body
        assert "capped at medium" in body

    def test_absent_verdict_invents_nothing(self, tmp_path):
        body = _generate(tmp_path, None)
        assert "Leg independence:" not in body

    def test_independent_legs_do_not_carry_the_caveat(self, tmp_path):
        body = _generate(tmp_path, IndependenceReport(
            marker_set_jaccard=0.10,
            shared_gene_tree_ratio=0.0,
            same_trimming_regime=False,
            ref_built_from_tested_markers=False,
        ))
        assert "Leg independence:" in body
        assert "NOT independent upstream" not in body

    def test_field_names_match_the_producer_exactly(self):
        # The consumer reads attributes the producer sets; a rename on either
        # Side would otherwise silently drop the whole section again.
        import inspect

        from markerfinder.modules import phylogenetic_inference

        assert list(inspect.signature(IndependenceReport).parameters) == [
            "marker_set_jaccard", "shared_gene_tree_ratio",
            "same_trimming_regime", "ref_built_from_tested_markers",
        ]
        assert "independence = IndependenceReport(" in inspect.getsource(
            phylogenetic_inference
        )
