"""Threshold plumbing + scan tests."""

import json

import pytest

from markerfinder.cli.self_test import interpreter_in_range  # Noqa: F401 (module import sanity)
from markerfinder.config import HGTConfig
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.pipeline_types import HGTEvaluation
from markerfinder.modules.hgt_scan import DEFAULT_BANDS, run_threshold_scan


def _evals(risks):
    return [
        HGTEvaluation(marker_id=f"m{i}", overall_risk=r, level=MarkerLevel.LEVEL_1)
        for i, r in enumerate(risks)
    ]


class TestCliThresholdSplit:
    def _build(self, **kwargs):
        """Build a config through the real parser + config_build path."""
        from markerfinder.cli.config_build import _build_pipeline_config
        from markerfinder.cli.parser import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--output", "out"])
        for key, value in kwargs.items():
            setattr(args, key, value)
        return _build_pipeline_config(args)

    def test_l1l2_cli_moves_only_l1l2(self):
        cfg = self._build(hgt_threshold_l1l2=0.40)
        thresholds = cfg.hgt_config.level_thresholds
        assert thresholds["level1_max"] == 0.40
        assert thresholds["level2_max"] == 0.60  # Untouched

    def test_l2l3_cli_moves_only_l2l3(self):
        cfg = self._build(hgt_threshold_l2l3=0.80)
        thresholds = cfg.hgt_config.level_thresholds
        assert thresholds["level1_max"] == 0.25  # Untouched
        assert thresholds["level2_max"] == 0.80

    def test_legacy_alias_maps_to_l1l2(self):
        cfg = self._build(hgt_threshold=0.30)
        assert cfg.hgt_config.level_thresholds["level1_max"] == 0.30
        assert cfg.hgt_config.level_thresholds["level2_max"] == 0.60

    def test_explicit_new_flag_wins_over_alias(self):
        cfg = self._build(hgt_threshold=0.30, hgt_threshold_l1l2=0.45)
        assert cfg.hgt_config.level_thresholds["level1_max"] == 0.45

    def test_max_markers_honored(self):
        cfg = self._build(max_markers=17)
        assert cfg.selection_config.max_markers == 17

    def test_hgt_mode_default_risk(self):
        cfg = self._build()
        assert cfg.hgt_config.hgt_mode == "risk"  #


class TestThresholdScan:
    def test_scan_output_schema(self, tmp_path):
        evals = _evals([0.05, 0.2, 0.3, 0.5, 0.9])
        results = run_threshold_scan(evals, str(tmp_path), "mf")
        tsv = tmp_path / "Phase5_reports" / "mf.threshold_scan.tsv"
        assert tsv.exists()
        header = tsv.read_text(encoding="utf-8").splitlines()[0]
        for col in ("band_level1_max", "band_level2_max", "n_L1", "n_L2",
                    "n_L3", "n_UNKNOWN", "n_selected", "n_at_boundary", "mean_risk"):
            assert col in header
        assert len(results) == len(DEFAULT_BANDS)

    def test_scan_regrades_without_external_tools(self, tmp_path):
        # Re-grading is pure arithmetic over stored risks —
        # The same risks must grade differently under different bands.
        evals = _evals([0.3])
        results = run_threshold_scan(evals, str(tmp_path), "mf")
        by_band = {(b.level1_max, b.level2_max): b for b in results}
        assert by_band[(0.15, 0.60)].n_L2 == 1
        assert by_band[(0.40, 0.60)].n_L1 == 1
        assert by_band[(0.25, 0.95)].n_L2 == 1

    def test_scan_does_not_rebuild_trees(self, tmp_path):
        """Guard: the scan must not invoke any external tool (counter stub)."""
        from unittest.mock import patch

        evals = _evals([0.1, 0.7])
        with patch("subprocess.run", side_effect=AssertionError("external tool called")):
            run_threshold_scan(evals, str(tmp_path), "mf")

    def test_unknown_untouched_by_bands(self, tmp_path):
        evals = _evals([0.1])
        unknown = HGTEvaluation(marker_id="u", overall_risk=0.0,
                                level=MarkerLevel.UNKNOWN)
        results = run_threshold_scan(evals + [unknown], str(tmp_path), "mf")
        for band in results:
            assert band.n_unknown == 1
            assert band.n_selected == 1  # Only the graded marker

    def test_boundary_counted(self, tmp_path):
        evals = _evals([0.25, 0.6])
        results = run_threshold_scan(evals, str(tmp_path), "mf")
        band = {(b.level1_max, b.level2_max): b for b in results}[(0.25, 0.60)]
        # Strict-less-than: risk==boundary lands in the upper band and is counted
        assert band.n_at_boundary == 2
        assert band.n_L1 == 0 and band.n_L2 == 1 and band.n_L3 == 1


class TestFarModeVisibility:
    def test_far_declaration_present_in_report(self):
        from markerfinder.models.pipeline_types import HGTReport
        from markerfinder.modules.report_generator import PlainTextReportGenerator

        report = HGTReport(total_markers=1, far_active=True)
        assert report.far_active is True

    def test_far_note_on_affected_marker_card(self):
        # Far mode active => affected marker notes record it.
        from markerfinder.config import HGTConfig
        from markerfinder.modules.hgt_filter import HGTDecisionEngine
        from markerfinder.models.pipeline_types import PhyloStepResult

        cfg = HGTConfig(adaptive_far_thresholds=True)
        cfg._far_distance_active = True
        engine = HGTDecisionEngine(cfg)
        ev = engine.evaluate_marker(
            "m", PhyloStepResult(overall_risk=0.7), far_active=True,
        )
        assert "far-distance" in ev.notes
        assert ev.decision_card["far_active"] is True
