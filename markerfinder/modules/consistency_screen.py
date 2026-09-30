"""Marker-level cross-consistency screen.

Thesis: grade each marker by whether the concatenation and coalescent legs
point the same way — "both sides agree and both clear the threshold" is
CONSISTENT, opposite sides is INCONSISTENT, either side under threshold is
INCONCLUSIVE. Never an average (AD).

GATE (AD): ``--hgt-mode consistency|hybrid`` refuses to
start unless the prerequisite work packages are observable as REAL ARTIFACTS
(fields, registry entries, files). No completion flags — a flag can be set
true early, an artifact cannot. 's expected-truth YAML is one of
the probes, so the gate unlocks automatically once that data lands.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import logging

logger = logging.getLogger(__name__)

from markerfinder.exceptions import UnsupportedCriterion
from markerfinder.models.evidence import MeasureState, Measurement

# Ladder migrated from the retracted paper (§B.2): (max normalized
# |ΔRF|, min quartet agreement). PROVISIONAL: must be calibrated on
# Truth data before consistency mode is used for real conclusions.
BOUNDS: Dict[int, Tuple[float, float]] = {
    1: (0.05, 0.90),
    2: (0.10, 0.80),
    3: (0.15, 0.70),
    4: (0.20, 0.60),
    5: (0.25, 0.50),
}


class ConsistencyGrade(Enum):
    CONSISTENT = "consistent"
    INCONSISTENT = "inconsistent"
    INCONCLUSIVE = "inconclusive"


@dataclass
class MarkerConsistency:
    marker_id: str
    concat_stance: Optional[Measurement]
    coalescent_stance: Optional[Measurement]
    grade: ConsistencyGrade
    stringency: int
    reasons: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Gate
# ---------------------------------------------------------------------------

def _probe(requirement: str) -> bool:
    """Evaluate one prerequisite by observing a real artifact."""
    if requirement == UFBOOT_OPTIONAL:
        import inspect
        from markerfinder.models.pipeline_types import SupermatrixResult

        hints = {f.name: f for f in __import__("dataclasses").fields(SupermatrixResult)}
        return "avg_ufboot" in hints and hints["avg_ufboot"].default is None
    if requirement == MEASURE_STATE:
        from dataclasses import fields
        from markerfinder.models.pipeline_types import PhyloStepResult

        names = {f.name for f in fields(PhyloStepResult)}
        return {"rf_state", "quartet_state", "monophyly_state"} <= names
    if requirement == REFERENCE_VALIDATOR:
        try:
            from markerfinder.utils.reference import validate_reference_tree  # Noqa: F401
            return True
        except ImportError:
            return False
    if requirement == ASSERTIONS_REGISTRY:
        try:
            from markerfinder.assertions import REGISTRY

            ids = {s.id for s in REGISTRY}
            return {f"A-{i:02d}" for i in range(1, 15)} <= ids
        except ImportError:
            return False
    if requirement == DECISION_CARD_SCHEMA2:
        from dataclasses import fields
        from markerfinder.models.pipeline_types import HGTEvaluation

        return "decision_card" in {f.name for f in fields(HGTEvaluation)}
    if requirement == DUAL_THRESHOLDS:
        from dataclasses import fields as dc_fields
        from markerfinder.config import HGTConfig

        field_types = {f.name for f in dc_fields(HGTConfig)}
        return "hgt_mode" in field_types and "level_thresholds" in field_types
    if requirement == PIS_COLUMN:
        from dataclasses import fields
        from markerfinder.models.marker import MarkerQualityScore

        return "pis" in {f.name for f in fields(MarkerQualityScore)}
    if requirement == BENCHMARK_EXPECTED_TRUTH:
        expected_dir = (
            Path(__file__).resolve().parent.parent.parent
            / "tests" / "benchmark" / "expected"
        )
        if not expected_dir.is_dir():
            return False
        return any(expected_dir.glob("*.yaml"))
    return False


# One name per prerequisite: the label the refusal shows and the probe key
# are the same string, so a message can never name a thing the code cannot
# look up.
UFBOOT_OPTIONAL = "ufboot_is_optional"
MEASURE_STATE = "measure_state_recorded"
REFERENCE_VALIDATOR = "reference_validator_exists"
ASSERTIONS_REGISTRY = "assertions_registry_exists"
DECISION_CARD_SCHEMA2 = "decision_card_schema2"
DUAL_THRESHOLDS = "dual_threshold_config"
PIS_COLUMN = "pis_column_written"
BENCHMARK_EXPECTED_TRUTH = "benchmark_expected_truth"

REQUIRED_PREREQS = (
    UFBOOT_OPTIONAL,
    MEASURE_STATE,
    REFERENCE_VALIDATOR,
    ASSERTIONS_REGISTRY,
    DECISION_CARD_SCHEMA2,
    DUAL_THRESHOLDS,
    PIS_COLUMN,
    BENCHMARK_EXPECTED_TRUTH,
)


def guard(mode: str, cfg: Optional[Any] = None) -> None:
    """Reject consistency/hybrid until every prerequisite artifact exists."""
    if mode == "risk":
        return
    missing = [name for name in REQUIRED_PREREQS if not _probe(name)]
    if missing:
        raise UnsupportedCriterion(
            f"--hgt-mode {mode} 需要以下条件先满足: {', '.join(missing)}. "
            "未经校准的交叉一致性判据会产生'统计显著、方向错误'的结果，"
            "因此缺少任一前置条件时拒绝启动。"
        )


# ---------------------------------------------------------------------------
# Screen
# ---------------------------------------------------------------------------

def _stance(
    gene_newick: str,
    ref_newick: str,
    rf_bound: float,
    quartet_bound: float,
) -> Tuple[Optional[bool], "Measurement"]:
    """Position of the gene tree relative to a reference tree.

    Strength = quartet agreement - normalized |ΔRF| (AD: no
    likelihood computation, no new external tool). ``agrees`` is None when
    the reference is illegal or metrics are unmeasurable.

    The second element is a full ``Measurement``, not a bare float: it carries
    the evidence state, the abstention reason, and ``n_obs`` -- how many quartets
    the agreement was actually computed from. ``screen`` used to rebuild a
    ``Measurement`` from the float alone, which forced ``n_obs=0`` and an empty
    reason on every leg, so the product could not tell a 1-quartet agreement from
    a 200-quartet one ("computed, stored, then dropped" again).
    """
    def _abstain(reason: str) -> Tuple[Optional[bool], Measurement]:
        return None, Measurement(
            None, MeasureState.NOT_MEASURABLE,
            method="quartet_minus_rf", reason=reason,
        )

    from markerfinder.utils.reference import validate_reference_tree

    verdict = validate_reference_tree(ref_newick, gene_newick)
    if not verdict.ok:
        return _abstain(f"reference_illegal [{verdict.code}]")

    from markerfinder.utils.tree_utils import calculate_rf_distance, get_quartet_topology

    rf, norm_rf = calculate_rf_distance(gene_newick, ref_newick)
    if rf is None or norm_rf is None:
        # Both must be present: the strength below is ``agreement - norm_rf``,
        # And ``float - None`` is a TypeError, not an unmeasurable result. The
        # Pair is produced by one call but only ``rf`` was checked, so a
        # Provider returning (0, None) would crash the run instead of reporting
        # NOT_MEASURABLE (; found by the strict type check).
        return _abstain("rf_unmeasurable")

    # Quartet agreement between the pair, sampled deterministically
    from itertools import combinations

    tips = sorted(_tips(gene_newick) & _tips(ref_newick))
    if len(tips) < 4:
        return _abstain("fewer than 4 shared tips")

    quartets = list(combinations(tips, 4))[:200]  # Deterministic prefix sample
    consistent = total = 0
    for quartet in quartets:
        q1 = get_quartet_topology(gene_newick, quartet)
        q2 = get_quartet_topology(ref_newick, quartet)
        if q1 and q2:
            total += 1
            if q1 == q2:
                consistent += 1
    if total == 0:
        return _abstain("no quartet evaluated")
    agreement = consistent / total

    strength = agreement - norm_rf
    agrees = (norm_rf <= rf_bound) and (agreement >= quartet_bound)
    return agrees, Measurement(
        float(strength),
        MeasureState.MEASURED,
        method="quartet_minus_rf",
        n_obs=total,
    )


def _tips(newick: str) -> Set[str]:
    from markerfinder.utils import etree as _etree

    return _etree.tip_set(newick)


def screen(
    marker_genes: Dict[str, str],
    concat_tree: str,
    astral_tree: str,
    *,
    stringency: int = 1,
    min_sites: int = 0,
    pis_map: Optional[Dict[str, int]] = None,
) -> List[MarkerConsistency]:
    """Grade every marker's cross-framework consistency.

     states three rules as an ordered list: "两侧同侧且都过阈 →
    consistent；两侧异侧 → inconsistent；任一侧不过阈 → inconclusive". The
    list overlaps on one input -- one leg above its threshold and on the
    opposite side, the other leg below its threshold satisfies both "两侧异侧"
    and "任一侧不过阈". The overlap resolves by the spec's own ordering:
    evaluated first-match-wins, the marker fails rule 1 (both-above does not
    hold), matches rule 2, and is graded INCONSISTENT. INCONCLUSIVE is
    reserved for "neither leg clears" or "a leg could not be measured at
    all" -- a marker the two frameworks actively disagree about is evidence
    against it, not absence of evidence. This is a reading of the spec's rule
    order, not an invention; re-ordering the rules would be a grading change,
    not a bug fix.
    """
    if stringency not in BOUNDS:
        raise ValueError(f"stringency must be 1..5, got {stringency!r}")
    rf_bound, quartet_bound = BOUNDS[stringency]
    pis_map = pis_map or {}

    results: List[MarkerConsistency] = []
    for marker_id, gene_newick in marker_genes.items():
        reasons: List[str] = []
        if min_sites and pis_map.get(marker_id, min_sites) < min_sites:
            results.append(
                MarkerConsistency(
                    marker_id, None, None, ConsistencyGrade.INCONCLUSIVE,
                    stringency, ["insufficient_sites"],
                )
            )
            continue

        concat_agrees, concat_m = _stance(
            gene_newick, concat_tree, rf_bound, quartet_bound
        )
        coal_agrees, coal_m = _stance(
            gene_newick, astral_tree, rf_bound, quartet_bound
        )
        reasons.extend(
            reason for reason in (concat_m.reason, coal_m.reason) if reason
        )

        if concat_agrees is None or coal_agrees is None:
            grade = ConsistencyGrade.INCONCLUSIVE
            reasons.append("either side unmeasurable")
        elif concat_agrees and coal_agrees:
            grade = ConsistencyGrade.CONSISTENT
        elif not concat_agrees and not coal_agrees:
            # Both sides individually weak — cannot distinguish direction
            grade = ConsistencyGrade.INCONCLUSIVE
            reasons.append("neither side above threshold")
        else:
            grade = ConsistencyGrade.INCONSISTENT
            reasons.append("sides disagree")

        results.append(
            MarkerConsistency(
                marker_id, concat_m, coal_m, grade, stringency, reasons,
            )
        )
    return results


def load_cog_category_map(path: Optional[str]) -> Dict[str, str]:
    """Read a ``marker_id<TAB>functional_category`` TSV.

    Returns ``{}`` when no path is given or the file is unreadable, and the
    caller renders that column as ``NA`` with a report note. Deliberately does
    NOT invent categories: filling this column from nothing is how a provenance
    column becomes a fabrication. Blank lines, ``#`` comments and a
    ``marker_id``/``category`` header row are skipped; single-column rows are
    ignored rather than mapped to an empty label.
    """
    if not path:
        return {}
    file_path = Path(path)
    if not file_path.exists():
        logger.warning(
            f"  --cog-category-map not found: {path} — functional_category "
            f"will be reported as NA (D-04)."
        )
        return {}
    mapping: Dict[str, str] = {}
    for raw in file_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 2 or not parts[0] or not parts[1]:
            continue
        if parts[0].lower() in ("marker_id", "marker"):
            continue
        mapping[parts[0]] = parts[1]
    if not mapping:
        logger.warning(
            f"  --cog-category-map {path} yielded no usable rows — column "
            f"stays NA (D-04)."
        )
    return mapping


def write_excluded_profile(
    results: List[MarkerConsistency],
    output_dir: str,
    prefix: str,
    *,
    pis_map: Optional[Dict[str, int]] = None,
    category_map: Optional[Dict[str, str]] = None,
) -> Path:
    """R4: profile every marker the consistency screen removed.

    Explicitly answers "did I filter out real signal?" — the retracted paper
    read its 95% removal as CONFIRMATION; the profile makes that visible.
    Missing COG category renders NA with a note.
    """
    out = Path(output_dir) / "Phase5_reports"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{prefix}.excluded_profile.tsv"

    pis_map = pis_map or {}
    category_map = category_map or {}

    def _num(stance: Optional[Measurement]) -> str:
        if stance is None or stance.value is None:
            return "NA"
        return f"{stance.value:.4f}"

    def _count(stance: Optional[Measurement]) -> str:
        if stance is None or not stance.n_obs:
            return "NA"
        return str(stance.n_obs)
    lines = [
        "marker_id\tgrade\tstringency\treasons\tpis\tfunctional_category\t"
        "category_note\tconcat_strength\tconcat_quartets\t"
        "coalescent_strength\tcoalescent_quartets\n"
    ]
    for r in results:
        if r.grade is ConsistencyGrade.CONSISTENT:
            continue
        category = category_map.get(r.marker_id)
        note = "NA (no --cog-category-map supplied, D-04)" if category is None else ""
        lines.append(
            "\t".join([
                r.marker_id,
                r.grade.value,
                str(r.stringency),
                ";".join(r.reasons) or "NA",
                str(pis_map.get(r.marker_id, "NA")),
                category or "NA",
                note,
                # The two legs' strengths and their observation base.
                # Absent is NA, never 0 -- "no quartet could be evaluated" and
                # "agreement over 0 quartets" are different statements.
                _num(r.concat_stance),
                _count(r.concat_stance),
                _num(r.coalescent_stance),
                _count(r.coalescent_stance),
            ]) + "\n"
        )
    target = path.resolve()
    target.write_text("".join(lines), encoding="utf-8", newline="\n")
    return target


@dataclass(frozen=True)
class ProfileArtifact:
    """The excluded-marker profile: where it landed and how many rows it holds."""

    path: Path
    rows_written: int


@dataclass(frozen=True)
class ScreenOutcome:
    """What the consistency stage produced (stages report their own work).

    ``profile.rows_written`` is counted here rather than by the caller because
    "a marker belongs in the excluded profile" is this module's rule -- reaching
    into ``r.grade`` from the orchestrator would silently fork that definition.
    """

    results: List[MarkerConsistency]
    grades: Dict[str, int]
    profile: ProfileArtifact


def summarize_grades(results: List[MarkerConsistency]) -> Dict[str, int]:
    """Tally per grade, always carrying all three keys (absent grades read 0)."""
    counts: Dict[str, int] = {grade.value: 0 for grade in ConsistencyGrade}
    for result in results:
        counts[result.grade.value] = counts.get(result.grade.value, 0) + 1
    return counts


def run_stage(
    results: List[MarkerConsistency],
    output_dir: str,
    prefix: str,
    *,
    pis_map: Optional[Dict[str, int]] = None,
    category_map: Optional[Dict[str, str]] = None,
) -> ScreenOutcome:
    """Write the excluded-marker profile and return it with the stage's summary."""
    path = write_excluded_profile(
        results, output_dir, prefix, pis_map=pis_map, category_map=category_map
    )
    return ScreenOutcome(
        results=results,
        grades=summarize_grades(results),
        profile=ProfileArtifact(
            path=path,
            rows_written=sum(
                1 for r in results if r.grade is not ConsistencyGrade.CONSISTENT
            ),
        ),
    )
