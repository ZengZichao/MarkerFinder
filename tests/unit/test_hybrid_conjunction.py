"""``--hgt-mode hybrid`` is the conjunction of both legs.

The specification fixes the semantics that the closure report had left as the
one open methodological decision: ``risk``（默认，行为不变）/ ``consistency``
（以 ``grade == CONSISTENT`` 为入选条件）/ ``hybrid``（**两者都要过**）. A
marker is includable only when it clears the risk level screen AND grades
CONSISTENT. These tests pin the three user-visible consequences: the per-marker
verdict names the leg that decided it, every verdict class is reachable, and
the gate still applies to hybrid exactly as it does to consistency.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from markerfinder.models.marker import MarkerLevel
from markerfinder.pipeline import _annotate_hybrid_verdicts
from markerfinder.modules.consistency_screen import ConsistencyGrade, guard
from markerfinder.exceptions import UnsupportedCriterion


@dataclass
class _Eval:
    marker_id: str
    level: MarkerLevel = MarkerLevel.LEVEL_1
    notes: str = ""


@dataclass
class _Report:
    marker_evaluations: list = field(default_factory=list)


def _graded_notes(evals, grades, passed):
    report = _Report(marker_evaluations=evals)
    _annotate_hybrid_verdicts(report, grades, passed)
    return {e.marker_id: e.notes for e in report.marker_evaluations}


def test_every_verdict_class_names_the_leg_that_decided_it():
    evals = [
        _Eval("m_pass", MarkerLevel.LEVEL_1),
        _Eval("m_level2", MarkerLevel.LEVEL_2),
        _Eval("m_inconsistent", MarkerLevel.LEVEL_1),
        _Eval("m_inconclusive", MarkerLevel.LEVEL_1),
        _Eval("m_ungraded", MarkerLevel.LEVEL_1),
        _Eval("m_risk_excluded", MarkerLevel.LEVEL_3),
    ]
    notes = _graded_notes(
        evals,
        grades={
            "m_pass": "consistent",
            "m_level2": "consistent",
            "m_inconsistent": "inconsistent",
            "m_inconclusive": "inconclusive",
        },
        passed={"m_pass", "m_level2"},
    )
    assert "passes both legs" in notes["m_pass"], notes["m_pass"]
    assert "passes both legs" in notes["m_level2"], notes["m_level2"]
    assert "EXCLUDED by the conjunction" in notes["m_inconsistent"], notes
    assert "grade=inconsistent" in notes["m_inconsistent"], notes
    assert "EXCLUDED by the conjunction" in notes["m_inconclusive"], notes
    assert "no consistency grade" in notes["m_ungraded"], notes["m_ungraded"]
    assert "risk leg (LEVEL_3)" in notes["m_risk_excluded"], notes
    assert "cannot rescue" in notes["m_risk_excluded"], notes


def test_a_risk_level3_marker_is_not_rescued_by_a_consistent_grade():
    """The conjunction cuts one way only: risk exclusion stands regardless of
    the consistency leg — both must pass, not either."""
    notes = _graded_notes(
        [_Eval("m", MarkerLevel.LEVEL_3)],
        grades={"m": "consistent"},
        passed={"m"},
    )
    assert "cannot rescue" in notes["m"], notes["m"]
    assert "passes both legs" not in notes["m"], notes["m"]


def test_existing_notes_are_preserved_not_replaced():
    report = _Report(marker_evaluations=[
        _Eval("m", MarkerLevel.LEVEL_1, notes="far: level2_max relaxed"),
    ])
    _annotate_hybrid_verdicts(report, {"m": "consistent"}, {"m"})
    note = report.marker_evaluations[0].notes
    assert note.startswith("far: level2_max relaxed; hybrid:"), note


def test_gate_still_governs_hybrid(monkeypatch):
    """Hybrid rides the same structural gate as consistency: with one
    prerequisite blocked, hybrid must refuse and name it by its own label."""
    import markerfinder.modules.consistency_screen as cs

    blocked = cs.BENCHMARK_EXPECTED_TRUTH
    real_probe = cs._probe
    monkeypatch.setattr(
        cs, "_probe",
        # One prerequisite is blocked; the rest are probed for real.
        lambda requirement: (
            requirement != blocked and real_probe(requirement)
        ),
    )
    with pytest.raises(UnsupportedCriterion) as err:
        guard("hybrid", None)
    assert blocked in str(err.value), err.value


def test_spec_decision_is_recorded_where_the_refusal_used_to_be():
    """The closure report's open item was 'hybrid 语义待定案'. The spec fixes it
    (两者都要过), so pipeline.py must cite the decision instead of
    carrying a refusal — a reader of the code should meet the rationale where
    the dead end used to be."""
    import inspect

    import markerfinder.pipeline as pipeline_module

    src = inspect.getsource(pipeline_module)
    assert "两者都要过" in src, "the spec citation for the hybrid semantics is gone"
    assert "语义尚未定案" not in src, "the old undecided-refusal text is still there"
