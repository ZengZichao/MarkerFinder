"""Evidence coverage must be computed, printed, and warn below the floor.

The field and the warning existed (``HGTReport.n_actually_measured``,
``evidence_coverage``, the ``EVIDENCE COVERAGE LOW`` block at the top of
``pipeline_summary.txt``, the ``--require-evidence-coverage`` flag) — and no test
named or asserted any of it. For a metric whose entire purpose is to stop a
reader from trusting a mostly-unmeasured run, "nobody checks it fires" is the same
position the baseline was in when 785 tests passed over a wrong conclusion.

The interesting failure mode is not a crash, it is a number that quietly means
something else, so the controls here check both directions: the warning must
appear below the floor AND must not appear above it, and the coverage must count
markers that were *actually measured*, not markers that were looked at.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from markerfinder.config import ReportConfig
from markerfinder.models.marker import MarkerLevel, SelectedMarkerSet
from markerfinder.models.pipeline_types import (
    HGTEvaluation,
    HGTReport,
    MarkerSelectionResult,
    PhylogeneticResult,
)
from markerfinder.modules.hgt_filter import HGTDecisionEngine


class _Runtime:
    duration = 1.0


def _engine() -> HGTDecisionEngine:
    from markerfinder.config import HGTConfig

    engine = HGTDecisionEngine.__new__(HGTDecisionEngine)
    engine.config = HGTConfig()
    return engine


def _evaluations(n_measured: int, n_unknown: int):
    evs = [
        HGTEvaluation(marker_id=f"measured_{i}", overall_risk=0.1,
                      level=MarkerLevel.LEVEL_1)
        for i in range(n_measured)
    ]
    evs += [
        HGTEvaluation(
            marker_id=f"unknown_{i}", overall_risk=0.0, level=MarkerLevel.UNKNOWN,
            step_details={"unscreened": {"reason": "no_gene_tree"}},
        )
        for i in range(n_unknown)
    ]
    return evs


def _summary(tmp_path, report: HGTReport, config: ReportConfig) -> str:
    from markerfinder.modules.report_generator import PlainTextReportGenerator

    marker_result = MarkerSelectionResult(
        marker_set=SelectedMarkerSet(
            markers=[e.marker_id for e in report.marker_evaluations] or ["M1"],
            occupancy_scores={}, quality_scores={},
        ),
        quality_scores={},
    )
    generator = PlainTextReportGenerator(config)
    out_dir = tmp_path / "Phase5_reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    generator._write_pipeline_summary(
        marker_result, report, PhylogeneticResult(), _Runtime(), out_dir, "mf",
    )
    return (out_dir / "mf.pipeline_summary.txt").read_text(encoding="utf-8")


def _report(n_measured: int, n_unknown: int) -> HGTReport:
    return _engine().generate_hgt_report(_evaluations(n_measured, n_unknown))


# ── the number itself ─────────────────────────────────────────────────────

def test_control_coverage_counts_only_actually_measured_markers():
    report = _report(n_measured=2, n_unknown=6)
    assert report.total_markers == 8
    assert report.n_actually_measured == 2, (
        "UNKNOWN markers were looked at but not measured; counting them would "
        "inflate the one number that tells a reader how much to trust the run"
    )
    assert report.evidence_coverage == pytest.approx(0.25)
    assert report.unknown_reason_counts.get("no_gene_tree") == 6


def test_control_perfect_coverage_is_one_and_empty_run_is_zero():
    assert _report(3, 0).evidence_coverage == pytest.approx(1.0)
    empty = _report(0, 0)
    assert empty.total_markers == 0
    assert empty.evidence_coverage == 0.0


# ── the product ───────────────────────────────────────────────────────────

def test_coverage_leads_the_hgt_section(tmp_path):
    config = ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
    text = _summary(tmp_path, _report(2, 6), config)
    assert "Evidence coverage: 0.25 (2/8 markers truly measured)" in text


def test_low_coverage_warns_at_the_top_of_the_summary(tmp_path):
    config = ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
    text = _summary(tmp_path, _report(2, 6), config)
    assert "EVIDENCE COVERAGE LOW" in text
    assert "only 25%" in text
    # "顶部" is observable, not adjectival: the warning must precede the
    # Input/Output block that a reader would otherwise take as the report body.
    assert text.index("EVIDENCE COVERAGE LOW") < text.index("Input/Output:")


def test_control_sufficient_coverage_does_not_warn(tmp_path):
    config = ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
    text = _summary(tmp_path, _report(9, 1), config)
    assert "EVIDENCE COVERAGE LOW" not in text
    assert "Evidence coverage: 0.90 (9/10 markers truly measured)" in text


def test_floor_is_configurable_and_default_is_half(tmp_path):
    assert ReportConfig().require_evidence_coverage == 0.5
    strict = ReportConfig(
        output_dir=str(tmp_path), output_prefix="mf", require_evidence_coverage=0.95
    )
    text = _summary(tmp_path, _report(9, 1), strict)
    assert "EVIDENCE COVERAGE LOW" in text, "0.90 < 0.95 must warn once the floor moves"


def test_control_empty_run_does_not_warn_about_coverage(tmp_path):
    """No markers means no HGT conclusions to distrust; the warning is for
    graded-but-unmeasured runs, and firing here would teach readers to ignore it."""
    config = ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
    text = _summary(tmp_path, _report(0, 0), config)
    assert "EVIDENCE COVERAGE LOW" not in text


# ── the CLI wiring ────────────────────────────────────────────────────────

def test_flag_reaches_the_report_config_end_to_end(tmp_path):
    """CLI -> config -> report, not just a flag that parses."""
    from markerfinder.cli.config_build import _build_pipeline_config
    from markerfinder.cli.parser import _build_parser

    parser = _build_parser()
    inbox = tmp_path / "genomes"
    inbox.mkdir()
    args = parser.parse_args([
        "-i", str(inbox), "-o", str(tmp_path / "out"),
        "--require-evidence-coverage", "0.8",
    ])

    assert args.require_evidence_coverage == 0.8
    config = _build_pipeline_config(args, parser)
    assert config.report_config.require_evidence_coverage == 0.8, (
        "the flag is only real if the object that renders the warning reads it"
    )
