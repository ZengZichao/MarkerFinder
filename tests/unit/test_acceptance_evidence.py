"""Acceptance evidence that needs no external tools —.

Three items on the checklist are provable in-process but had no live
demonstration yet:

* "扫描不重跑外部工具（桩计数器证明）" — the scan must be pure re-grading.
* "旧 SCHEMA_VERSION=1 的 state 在 --resume 时明确失败并给指引".
* "stringency 1..5 给出单调趋势" — the ladder must tighten.

Each test refuses to pass vacuously: where a check depends on a measurement
that this environment cannot make, it reports NOT EXECUTED rather than green.
"""

from __future__ import annotations

import inspect
import subprocess

import pytest

from markerfinder.models.marker import MarkerLevel
from markerfinder.modules.consistency_screen import BOUNDS, screen
from markerfinder.modules.hgt_scan import run_threshold_scan
from markerfinder.utils.state_codec import _assert_compatible


def _ev(marker_id, risk, level=MarkerLevel.LEVEL_2):
    from types import SimpleNamespace
    return SimpleNamespace(marker_id=marker_id, overall_risk=risk, level=level)


class TestScanMakesNoExternalCalls:
    """The threshold scan re-grades stored risks, nothing else."""

    def test_scan_runs_with_every_external_tool_call_boobytrapped(self, tmp_path):
        evaluations = [
            _ev("M1", 0.10), _ev("M2", 0.30), _ev("M3", 0.55),
            _ev("M4", 0.95, level=MarkerLevel.UNKNOWN),
        ]

        def _explode(*args, **kwargs):
            raise AssertionError(
                f"threshold scan invoked an external tool: {args[:2]!r}"
            )

        original = subprocess.run
        subprocess.run = _explode
        try:
            bands = run_threshold_scan(evaluations, str(tmp_path), "mf")
        finally:
            subprocess.run = original

        # The scan must still have produced its own product.
        assert (tmp_path / "Phase5_reports" / "mf.threshold_scan.tsv").exists()
        assert len(bands) >= 2
        # UNKNOWN markers are counted, never re-graded into a band.
        assert all(band.n_unknown == 1 for band in bands)


class TestOldStateIsNeverSilentlyAccepted:
    """Schema v1 state must fail loudly with a way forward."""

    def test_v1_state_raises_and_tells_the_user_how_to_proceed(self):
        with pytest.raises(Exception) as exc:
            _assert_compatible({"schema_version": 1, "completed_steps": ["scan"]})
        message = str(exc.value)
        assert "schema" in message.lower(), message
        assert any(hint in message for hint in ("--redo", "--resume", "rerun", "重跑")), (
            f"v1 state fails but gives no actionable way forward: {message!r}"
        )

    def test_current_version_is_accepted(self):
        from markerfinder.utils.state_codec import SCHEMA_VERSION
        assert SCHEMA_VERSION == 2
        _assert_compatible({"schema_version": SCHEMA_VERSION})


class TestStringencyLadder:
    """ Must move in one direction, monotonically."""

    def test_bounds_ladder_tightens_from_th1_to_th5(self):
        assert set(BOUNDS) == {1, 2, 3, 4, 5}, BOUNDS
        rf_bounds = [BOUNDS[i][0] for i in range(1, 6)]
        quartet_bounds = [BOUNDS[i][1] for i in range(1, 6)]
        # Is the strictest rung: smallest RF tolerance, highest agreement
        # Requirement. Both edges must move away from strictness monotonically.
        assert rf_bounds == sorted(rf_bounds) and len(set(rf_bounds)) == 5
        assert quartet_bounds == sorted(quartet_bounds, reverse=True)
        assert len(set(quartet_bounds)) == 5

    def test_grading_a_real_set_is_monotone_in_stringency(self):
        tips = ["A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2"]
        concat = "(((A1,A2),(B1,B2)),((C1,C2),(D1,D2)));"
        astral = "(((A1,A2),(C1,C2)),((B1,B2),(D1,D2)));"
        gene_trees = {
            "matches_concat": concat,
            "matches_astral": astral,
            "in_between": "(((A1,A2),(B1,B2)),((C1,D1),(D2,C2)));",
        }
        counts = []
        for stringency in range(1, 6):
            results = screen(
                gene_trees, concat, astral, stringency=stringency
            )
            measured = [
                r for r in results
                if r.concat_stance is not None
                and getattr(r.concat_stance, "state", None) is not None
                and str(getattr(r.concat_stance.state, "value", "")).lower()
                == "measured"
            ]
            if not measured:
                pytest.skip(
                    "NOT EXECUTED: RF/quartet unmeasurable in this interpreter "
                    "(needs ete3, which is importable on every supported "
                    "interpreter once `markerfinder` is loaded — so this means a "
                    "broken install) — a vacuous monotonicity check would be "
                    "meaningless."
                )
            counts.append(
                sum(1 for r in results
                    if str(getattr(r.grade, "value", r.grade)) == "consistent")
            )
        # A looser rung can never certify fewer markers as consistent.
        assert counts == sorted(counts), (
            f"consistent counts must be non-decreasing as stringency loosens: {counts}"
        )
