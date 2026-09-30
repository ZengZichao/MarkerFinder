"""Output-numeric assertion layer.

The assertion layer is deliberately separate from ``validation.py`` (AD):
input validation rejects bad *inputs*; this layer observes computed *outputs*
and aborts (exit code 4) when a number is impossible or meaningless. It is
observe-only: it never changes a marker's grade.

Every registered assertion carries a must-fail fixture: a negative
state that MUST trip it. ``--check`` runs every negative fixture and requires
each assertion to fire exactly once — an assertion that cannot go red is
treated as absent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from markerfinder._version import __version__
from markerfinder.models.evidence import MeasureState
from markerfinder.models.pipeline_types import PhyloStepResult


@dataclass(frozen=True)
class AssertionSpec:
    """One registered output assertion."""

    id: str                                   # E.g. "A-01"
    name: str
    applies_to: str                           # "PhyloStepResult" | "run" |...
    check: Callable[[Dict[str, Any]], Optional[str]]    # None=pass, str=violation text
    severity: str                             # "warn" | "fail"
    rationale: str                            # Why the check exists
    must_fail_fixture: str                    # Tests/fixtures/must_fail/<id>.json
    provisional: bool = False                 # Threshold pending calibration


@dataclass
class AssertionResult:
    assertion_id: str
    severity: str
    passed: bool
    detail: str = ""
    provisional: bool = False


@dataclass
class AssertionReport:
    results: List[AssertionResult] = field(default_factory=list)
    mode: str = "run"
    code_version: str = ""
    # Set when a FAIL was suppressed by --allow-assertion-failure;
    # The bypass must be visible in the written report.
    bypassed: bool = False

    @property
    def fail_count(self) -> int:
        return sum(1 for r in self.results if not r.passed and r.severity == "fail")

    @property
    def warn_count(self) -> int:
        return sum(1 for r in self.results if not r.passed and r.severity == "warn")

    @property
    def pass_count(self) -> int:
        return sum(1 for r in self.results if r.passed)

    def to_tsv_rows(self) -> List[List[str]]:
        header = ["assertion_id", "name", "severity", "result", "detail", "provisional"]
        rows = [header]
        name_by_id = {spec.id: spec.name for spec in REGISTRY}
        for r in self.results:
            rows.append([
                r.assertion_id,
                name_by_id.get(r.assertion_id, ""),
                r.severity,
                "PASS" if r.passed else ("FAIL" if r.severity == "fail" else "WARN"),
                r.detail.replace("\t", " ").replace("\n", " "),
                "true" if r.provisional else "false",
            ])
        return rows


# ---------------------------------------------------------------------------
# Check helpers
# ---------------------------------------------------------------------------

def _range_check(value: Optional[float], low: float, high: float, label: str) -> Optional[str]:
    """Interval check that rejects NaN/inf explicitly ( trap)."""
    if value is None:
        return None
    if not math.isfinite(value):
        return f"{label} is not finite: {value!r}"
    if not (low <= value <= high):
        return f"{label}={value!r} outside [{low}, {high}]"
    return None


def _steps(state: Dict[str, Any]) -> List[PhyloStepResult]:
    return list(state.get("phylo_steps") or [])


# ---------------------------------------------------------------------------
# Individual checks (A-01.. A-14)
# ---------------------------------------------------------------------------

def _check_a01(state: Dict[str, Any]) -> Optional[str]:
    for item in state.get("rendered_fields") or []:
        if item.get("value") is None:
            rendered = str(item.get("rendered"))
            if rendered not in ("NA", "N/A"):
                return (
                    f"optional field {item.get('field')!r} rendered as {rendered!r}; "
                    "missing values must render NA (TSV) / N/A (HTML)"
                )
    return None


def _check_a02(state: Dict[str, Any]) -> Optional[str]:
    for step in _steps(state):
        v = step.normalized_rf
        if v is None:
            continue
        problem = _range_check(v, 0.0, 1.0, f"normalized_rf[{step.gene_id}]")
        if problem:
            return problem
    return None


def _check_a03(state: Dict[str, Any]) -> Optional[str]:
    for step in _steps(state):
        n_shared = state.get("shared_tips", {}).get(step.gene_id)
        if n_shared is not None and n_shared < 4 and step.quartet_consistency is not None:
            return (
                f"quartet_consistency has a value for {step.gene_id} with "
                f"{n_shared} shared tips (<4)"
            )
        # Also treat a not-measurable state carrying a value as a violation
        if (
            step.quartet_state in (MeasureState.NOT_MEASURABLE, MeasureState.REJECTED)
            and step.quartet_consistency is not None
        ):
            return f"quartet_consistency has a value although state={step.quartet_state.value} for {step.gene_id}"
    return None


def _check_a04(state: Dict[str, Any]) -> Optional[str]:
    for rec in state.get("monophyly_records") or []:
        if rec.get("proportion") is None:
            continue
        if rec.get("n_total", 0) < 2:
            return (
                f"monophyly proportion {rec.get('proportion')!r} reported with "
                f"n_total={rec.get('n_total')} (<2)"
            )
    return None


def _check_a05(state: Dict[str, Any]) -> Optional[str]:
    for ev in state.get("hgt_evaluations") or []:
        problem = _range_check(ev.get("overall_risk"), 0.0, 1.0, f"overall_risk[{ev.get('marker_id')}]")
        if problem:
            return problem
    return None


def _check_a06(state: Dict[str, Any]) -> Optional[str]:
    values = [
        s.normalized_rf for s in _steps(state)
        if s.normalized_rf is not None
    ]
    if len(values) < 5:
        return None  # Distribution probe needs a minimal sample
    saturated = sum(1 for v in values if v in (0.0, 0.5, 1.0))
    if saturated / len(values) > 0.20:
        return (
            f"{saturated}/{len(values)} normalized_rf values are exactly "
            f"0.0/0.5/1.0 (>20%) — placeholder or saturated distribution"
        )
    return None


def _check_a07(state: Dict[str, Any]) -> Optional[str]:
    for aln in state.get("alignments") or []:
        if aln.get("effective_columns", 1) <= 0:
            return f"alignment for {aln.get('marker_id')} is empty after trimming"
    return None


def _check_a08(state: Dict[str, Any]) -> Optional[str]:
    supports = [s for s in (state.get("astral_supports") or []) if s is not None]
    bad = [s for s in supports if not (0.0 < s <= 1.0)]
    if bad:
        return f"ASTRAL LPP values outside (0,1]: {bad[:3]}"
    if supports:
        constant = sum(1 for s in supports if s in (0.0, 1.0))
        if constant / len(supports) > 0.9:
            return f"{constant}/{len(supports)} LPP values are exactly 0 or 1 — degenerate supports"
    return None


def _check_a09(state: Dict[str, Any]) -> Optional[str]:
    for v in state.get("ufboot_values") or []:
        if v is None:
            continue
        if not (0.0 <= v <= 100.0):
            return f"UFBoot value {v!r} outside [0, 100]"
    return None


def _check_a10(state: Dict[str, Any]) -> Optional[str]:
    for diff in state.get("lnl_diffs") or []:
        if diff is not None and abs(diff) > 30:
            return (
                f"|dlnL|={abs(diff):.1f} exceeds 30 (retraction Error-B scale: "
                f"typical <10, thousands = wrong tree compared)"
            )
    return None


def _check_a11(state: Dict[str, Any]) -> Optional[str]:
    coverage = state.get("evidence_coverage")
    n_markers = state.get("n_markers", 0)
    if coverage is None:
        if n_markers == 0:
            return None  # No markers: nothing to assert, and never NaN
        return "evidence_coverage missing although markers were screened"
    if not math.isfinite(coverage):
        return f"evidence_coverage is not finite: {coverage!r}"
    if coverage < 0.5:
        return f"evidence coverage {coverage:.2f} below 0.5 (provisional floor)"
    return None


def _check_a12(state: Dict[str, Any]) -> Optional[str]:
    probe = state.get("discrimination_probe") or {}
    conflict = probe.get("conflict_risk")
    identical = probe.get("identical_risk")
    if conflict is None or identical is None:
        return None  # Probe not run in this mode
    if conflict == identical:
        return (
            f"zero discrimination: conflicting and identical references both "
            f"produce overall_risk={conflict!r} (zero-discrimination)"
        )
    return None


def _check_a13(state: Dict[str, Any]) -> Optional[str]:
    thresholds = state.get("level_thresholds") or {}
    l1 = thresholds.get("level1_max")
    l2 = thresholds.get("level2_max")
    if l1 is None or l2 is None:
        return "level thresholds missing"
    if not (0.0 <= l1 < l2 <= 1.0):
        return f"invalid thresholds: 0 <= level1_max({l1}) < level2_max({l2}) <= 1 violated"
    return None


def _check_a14(state: Dict[str, Any]) -> Optional[str]:
    if not state.get("far_active"):
        return None
    entries = state.get("far_card_entries") or []
    if not entries:
        return "far-distance mode active but no affected marker card records it"
    return None


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

REGISTRY: List[AssertionSpec] = [
    AssertionSpec(
        "A-01", "None rendered as number", "report",
        _check_a01, "fail",
        "None must never render as a number in a product", "A-01_none_rendered_as_number.json",
    ),
    AssertionSpec(
        "A-02", "normalized_rf in [0,1] and finite", "PhyloStepResult",
        _check_a02, "fail",
        "NaN passes a naive range check, so finiteness is checked first", "A-02_rf_out_of_range.json",
    ),
    AssertionSpec(
        "A-03", "no quartet value below 4 shared tips", "PhyloStepResult",
        _check_a03, "fail",
        "a quartet is only defined on four shared tips", "A-03_quartet_on_three_tips.json",
    ),
    AssertionSpec(
        "A-04", "monophyly n_total >= 2", "PhyloStepResult",
        _check_a04, "fail",
        "a proportion over one taxon says nothing", "A-04_monophyly_singleton_taxon.json",
    ),
    AssertionSpec(
        "A-05", "overall_risk in [0,1]", "HGTEvaluation",
        _check_a05, "fail",
        "risk is an interval quantity, negative is impossible", "A-05_negative_overall_risk.json",
    ),
    AssertionSpec(
        "A-06", "placeholder/saturation distribution probe", "run",
        _check_a06, "fail",
        ">20% exact 0.0/0.5/1.0 means the number was not measured", "A-06_all_half_rf.json",
    ),
    AssertionSpec(
        "A-07", "trimmed alignment non-empty", "alignment",
        _check_a07, "fail",
        "the gtdb_tk_markers fallback can yield an empty alignment", "A-07_empty_alignment.json",
    ),
    AssertionSpec(
        "A-08", "ASTRAL LPP in (0,1], not degenerate", "run",
        _check_a08, "warn",
        "a constant LPP carries no information", "A-08_constant_lpp.json",
    ),
    AssertionSpec(
        "A-09", "UFBoot in [0,100]", "run",
        _check_a09, "fail",
        "ufboot out of its bootstrap range is a misread column", "A-09_ufboot_120.json",
    ),
    AssertionSpec(
        "A-10", "|dlnL| magnitude guard (reserved)", "run",
        _check_a10, "fail",
        "provisional: no likelihood path exists in this version",
        "A-10_delta_lnL_4000.json", provisional=True,
    ),
    AssertionSpec(
        "A-11", "evidence coverage >= 0.5 (provisional)", "run",
        _check_a11, "warn",
        "isfinite-guarded; the coverage floor is not yet calibrated",
        "A-11_low_coverage.json", provisional=True,
    ),
    AssertionSpec(
        "A-12", "zero-discrimination probe", "run",
        _check_a12, "fail",
        "an output that never varies is not a measurement", "A-12_zero_discrimination.json",
    ),
    AssertionSpec(
        "A-13", "thresholds 0 <= L1 < L2 <= 1", "config",
        _check_a13, "fail",
        "thresholds must be ordered and inside [0,1]", "A-13_bad_thresholds.json",
    ),
    AssertionSpec(
        "A-14", "far mode recorded on affected cards", "run",
        _check_a14, "fail",
        "far-mode relaxation must be visible on the affected cards", "A-14_far_not_on_cards.json",
    ),
]


def _good_state() -> Dict[str, Any]:
    """A canonical, healthy run state that every assertion must accept."""
    steps = [
        PhyloStepResult(gene_id="m1", rf_distance=0, normalized_rf=0.0,
                        quartet_consistency=1.0),
        PhyloStepResult(gene_id="m2", rf_distance=1, normalized_rf=0.25,
                        quartet_consistency=0.7),
    ]
    return {
        "phylo_steps": steps,
        "shared_tips": {"m1": 4, "m2": 8},
        "monophyly_records": [{"proportion": 0.8, "n_total": 6, "n_mono": 5}],
        "hgt_evaluations": [{"marker_id": "m1", "overall_risk": 0.1}],
        "alignments": [{"marker_id": "m1", "effective_columns": 120}],
        "astral_supports": [0.95, 0.87, 0.7],
        "ufboot_values": [95.0, 88.0],
        "lnl_diffs": [],
        "evidence_coverage": 0.9,
        "n_markers": 10,
        "discrimination_probe": {"conflict_risk": 1.0, "identical_risk": 0.0},
        "level_thresholds": {"level1_max": 0.25, "level2_max": 0.60},
        "far_active": False,
        "far_card_entries": [],
        "rendered_fields": [{"field": "avg_ufboot", "value": None, "rendered": "NA"}],
    }


def run_assertions(state: Dict[str, Any], *, mode: str = "run") -> AssertionReport:
    """Run every registered assertion against ``state``.

    An assertion that raises is a violation, never a silent pass.
    """
    results: List[AssertionResult] = []
    for spec in REGISTRY:
        try:
            violation = spec.check(state)
        except Exception as e:  # Noqa: BLE001 — assertion crash must be visible
            violation = f"assertion raised {type(e).__name__}: {e}"
        results.append(
            AssertionResult(
                assertion_id=spec.id,
                severity=spec.severity,
                passed=violation is None,
                detail=violation or "",
                provisional=spec.provisional,
            )
        )
    return AssertionReport(results=results, mode=mode, code_version=__version__)


# ---------------------------------------------------------------------------
# Must-fail control: build a violating state per assertion id
# ---------------------------------------------------------------------------

def build_negative_state(assertion_id: str) -> Dict[str, Any]:
    """Construct a state that MUST trip ``assertion_id``.

    Used by ``--check`` and the CI test suite: an assertion whose negative
    fixture does not fire is treated as absent. JSON fixtures under
    ``tests/fixtures/must_fail/`` document each case; this builder materialises
    the equivalent in-memory state (trees/NaN cannot live in JSON).
    """
    state = _good_state()  # Start healthy; overlay only the poisoning fields

    if assertion_id == "A-01":
        state["rendered_fields"] = [
            {"field": "avg_ufboot", "value": None, "rendered": "0.0000"},
        ]
    elif assertion_id == "A-02":
        s = PhyloStepResult(gene_id="m1", rf_distance=5,
                            normalized_rf=1.3, quartet_consistency=0.5)
        state["phylo_steps"] = [s]
        nan_step = PhyloStepResult(gene_id="m2", rf_distance=5,
                                   normalized_rf=float("nan"), quartet_consistency=0.5)
        state["phylo_steps"].append(nan_step)
        state["phylo_steps"].append(
            PhyloStepResult(gene_id="m3", rf_distance=5,
                            normalized_rf=float("inf"), quartet_consistency=0.5))
    elif assertion_id == "A-03":
        s = PhyloStepResult(gene_id="m1", rf_distance=0, normalized_rf=0.0,
                            quartet_consistency=1.0)
        state["phylo_steps"] = [s]
        state["shared_tips"] = {"m1": 3}
    elif assertion_id == "A-04":
        state["monophyly_records"] = [{"proportion": 1.0, "n_total": 1, "n_mono": 1}]
    elif assertion_id == "A-05":
        state["hgt_evaluations"] = [{"marker_id": "m1", "overall_risk": -0.1}]
    elif assertion_id == "A-06":
        state["phylo_steps"] = [
            PhyloStepResult(gene_id=f"m{i}", normalized_rf=0.5)
            for i in range(10)
        ]
    elif assertion_id == "A-07":
        state["alignments"] = [{"marker_id": "m1", "effective_columns": 0}]
    elif assertion_id == "A-08":
        state["astral_supports"] = [1.0, 1.0, 1.0, 1.0]
    elif assertion_id == "A-09":
        state["ufboot_values"] = [120.0]
    elif assertion_id == "A-10":
        state["lnl_diffs"] = [4000.0]
    elif assertion_id == "A-11":
        state["evidence_coverage"] = 0.2
        state["n_markers"] = 10
    elif assertion_id == "A-12":
        state["discrimination_probe"] = {"conflict_risk": 0.25, "identical_risk": 0.25}
    elif assertion_id == "A-13":
        state["level_thresholds"] = {"level1_max": 0.7, "level2_max": 0.60}
    elif assertion_id == "A-14":
        state["far_active"] = True
        state["far_card_entries"] = []
    else:
        raise KeyError(f"no negative fixture registered for {assertion_id!r}")
    return state
