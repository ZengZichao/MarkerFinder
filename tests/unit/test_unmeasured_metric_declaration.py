""" (P0): a run must declare WHICH metrics it did not measure, and why.

The requirement's own evidence site is ``hgt_filter.py:52-54`` — the branch that
used to write ``normalized_rf = 0.5`` when the RF distance could not be computed.
The placeholder is gone, but the other half of was never tested at
all: "日志与报告必须显式声明哪些指标因此未被测量". Before ``utils/unmeasured.py`` the
run said "evidence_coverage: 0.30" (a fraction that names no metric and no cause)
and rendered per-marker ``NA`` cells. Nobody could tell that the reason was
"ete3 unusable here" rather than "these markers have no gene tree" — two different
remedies, one undecidable summary.

These tests drive the real measurement functions, so the ledger is populated by
the code paths that actually fail, not by a test calling ``record`` by hand.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from markerfinder.config import HGTConfig, ReportConfig
from markerfinder.models.marker import SelectedMarkerSet
from markerfinder.models.pipeline_types import (
    HGTEvaluation,
    HGTReport,
    MarkerSelectionResult,
    PhylogeneticResult,
)
from markerfinder.models.tree import Tree
from markerfinder.modules import hgt_filter as hgt_filter_module
from markerfinder.modules.report_generator import PlainTextReportGenerator
from markerfinder.utils import etree as etree_module
from markerfinder.utils import unmeasured


class _Runtime:
    duration = 1.0


@pytest.fixture(autouse=True)
def _clean_ledger():
    """The ledger is process-global; a leak between tests would fake evidence."""
    unmeasured.clear()
    yield
    unmeasured.clear()


def _summary_text(tmp_path) -> str:
    marker = "M1"
    marker_result = MarkerSelectionResult(
        marker_set=SelectedMarkerSet(
            markers=[marker], occupancy_scores={marker: 0.9}, quality_scores={},
        ),
        quality_scores={},
    )
    hgt_report = HGTReport(
        marker_evaluations=[HGTEvaluation(marker_id=marker, overall_risk=0.1)],
        total_markers=1,
        evidence_coverage=0.0,
    )
    generator = PlainTextReportGenerator(
        ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
    )
    out_dir = tmp_path / "Phase5_reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    generator._write_pipeline_summary(
        marker_result, hgt_report, PhylogeneticResult(), _Runtime(), out_dir, "mf",
    )
    return (out_dir / "mf.pipeline_summary.txt").read_text(encoding="utf-8")


def _unavailable_ete3():
    from markerfinder.exceptions import PhyloToolUnavailable

    def _raise():
        raise PhyloToolUnavailable("simulated: ete3 not importable here")

    return _raise


# ── ledger mechanics (with planted controls) ──────────────────────────────

def test_control_ledger_counts_groups_collapses_and_clears():
    for marker in ["m2", "m1", "m3", "m4"]:
        unmeasured.record("metric_a", "cause one", marker)
    unmeasured.record("metric_b", "cause two")  # Run-level, no marker
    rows = unmeasured.rollup()
    assert [(r[0], r[2]) for r in rows] == [("metric_a", 4), ("metric_b", 0)]
    # Marker ids are named in sorted order, not insertion order.
    assert rows[0][3] == ["m1", "m2", "m3", "m4"]
    text = "\n".join(unmeasured.render_lines())
    assert "metric_a: cause one [4 marker(s): m1, m2, m3, m4]" in text
    assert "metric_b: cause two [whole run]" in text
    assert "\n".join(unmeasured.render_lines()) == text, "render must be stable"
    unmeasured.clear()
    assert unmeasured.rollup() == []
    assert unmeasured.render_lines() == []


def test_control_long_marker_lists_are_truncated_not_dropped():
    for index in range(unmeasured.MAX_NAMED_MARKERS + 3):
        unmeasured.record("metric_a", "cause", f"m{index:02d}")
    text = "\n".join(unmeasured.render_lines())
    assert "(+3 more)" in text, text
    assert "m00" in text and "m07" in text


# ── the real call sites declare themselves ────────────────────────────────

def test_rf_distance_failure_is_declared_by_name_and_cause():
    from markerfinder.utils.tree_utils import calculate_rf_distance

    rf, norm = calculate_rf_distance("not a tree", "also not a tree")
    assert (rf, norm) == (None, None), "an unmeasured quantity stays None"
    metrics = {row[0] for row in unmeasured.rollup()}
    assert "rf_distance (Robinson-Foulds)" in metrics, metrics


def test_quartet_topology_failure_is_declared():
    from markerfinder.utils.tree_utils import get_quartet_topology

    assert get_quartet_topology("broken((", ("A", "B", "C", "D")) is None
    metrics = {row[0] for row in unmeasured.rollup()}
    assert "quartet_topology" in metrics, metrics


def test_monophyly_screen_names_the_metric_and_the_marker(monkeypatch):
    monkeypatch.setattr(etree_module, "require_ete3", _unavailable_ete3())
    detector = hgt_filter_module.PhylogeneticHGTDetector(HGTConfig())
    taxonomy = {
        "A": {"genus": "G1"}, "B": {"genus": "G1"},
        "C": {"genus": "G2"}, "D": {"genus": "G2"},
    }
    assert detector.detect_monophyly("marker_7", "((A,B),(C,D));", taxonomy,
                                    level="genus", threshold=0.5) is None
    rows = unmeasured.rollup()
    named = {row[0]: row for row in rows}
    key = "hgt_monophyly_proportion (MAD-rooted)"
    assert key in named, rows
    assert "ete3 unusable" in named[key][1]
    assert named[key][3] == ["marker_7"]


def test_gene_tree_support_declares_why_the_pis_input_is_missing():
    from markerfinder.pipeline import _gene_tree_mean_support

    assert _gene_tree_mean_support(Tree(newick="((A,B),(C,D));")) is None
    rows = unmeasured.rollup()
    metric = "gene_tree_mean_support (PIS back-fill input)"
    assert any(row[0] == metric for row in rows), rows
    assert any("no readable internal support label" in row[1] for row in rows), rows


# ── the declaration reaches both products ─────────────────────────────────

def test_control_summary_stays_clean_when_nothing_failed(tmp_path):
    text = _summary_text(tmp_path)
    assert "METRICS NOT MEASURED THIS RUN" not in text, (
        "a declaration invented out of thin air is worse than none"
    )


def test_summary_product_declares_missing_metrics_with_causes(tmp_path):
    from markerfinder.utils.tree_utils import calculate_rf_distance

    calculate_rf_distance("not a tree", "also not a tree")
    unmeasured.record("checkm_completeness", "no CheckM results supplied", "M1")

    text = _summary_text(tmp_path)
    assert "METRICS NOT MEASURED THIS RUN" in text
    assert "rf_distance (Robinson-Foulds)" in text
    assert "checkm_completeness: no CheckM results supplied [1 marker(s): M1]" in text
    # The declaration is a statement about absence: nothing here may hand back a
    # Neutral number for the missing metrics.
    assert "NA" in text or "0.30" in text


def test_pipeline_log_carries_the_same_declaration(caplog):
    """ Names the log as a carrier too: a summary nobody opens is not a
    declaration. The helper is what run calls, so this tests the real
    emission path instead of grepping for it."""
    unmeasured.record("quartet_topology", "tip cap exceeded", "M1")
    with caplog.at_level(logging.WARNING, logger="markerfinder.utils.unmeasured"):
        count = unmeasured.declare_in_log(logging.getLogger("markerfinder.pipeline"))
    assert count == 4, count  # Header, sub-header, one metric row, trailing blank
    messages = [rec.message for rec in caplog.records]
    assert any("METRICS NOT MEASURED THIS RUN" in m for m in messages), messages
    assert any("quartet_topology: tip cap exceeded [1 marker(s): M1]" in m for m in messages)
    assert all(rec.levelno == logging.WARNING for rec in caplog.records), messages


def test_run_still_declares_through_the_helper():
    """Wiring lock: run must not go back to rendering only into the report."""
    source = (
        Path(__file__).resolve().parents[2] / "markerfinder" / "pipeline.py"
    ).read_text(encoding="utf-8")
    assert "unmeasured.declare_in_log(logger)" in source
    assert "unmeasured.clear()" in source
