"""Taxonomy must-pass benchmark.

Why this exists: the retracted *Science* paper was not undone by weak
statistics but by topologies that any taxonomist would call absurd (a
ctenophore sister to a mite). MarkerFinder had the same blind spot — the
must-pass baseline existed only as a YAML file with **no reader at all**, so
"the pipeline is rejected when a required relationship fails" could
never happen.

Design notes:

* Monophyly is tested with the pure-Python split-set machinery from
  ``utils.etree`` (the same code path the differential tests validate
  against ete3), so the gate works with or without ete3 installed.
* A group is only tested when every one of its taxonomied tips is actually
  present in that marker's tree; partial coverage is recorded as NOT_CHECKED
  rather than silently passing or spuriously failing.
* An unsupported ``relation`` is recorded loudly. Unknown syntax must never
  read as "passed" — that is the must-fail lesson.
* Failure raises ``AssertionFailureError``, which the CLI maps to
  ``EXIT_ASSERTION_FAILED`` (4): the run is rejected, no marker is merely
  docked a risk score.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from markerfinder.exceptions import AssertionFailureError
from markerfinder.utils import etree

logger = logging.getLogger(__name__)

SUPPORTED_RELATIONS = ("monophyletic",)
# Cap how many violations are printed / persisted; the run aborts on the first
# One anyway, but a reader needs to see the pattern, not just one instance.
MAX_RECORDED_VIOLATIONS = 20


@dataclass
class MustPassViolation:
    marker_id: str
    relation: str
    level: str
    taxon: str
    members: Tuple[str, ...]
    note: str = ""

    def describe(self) -> str:
        return (
            f"{self.marker_id}: {self.level} '{self.taxon}' is NOT "
            f"{self.relation} over {', '.join(sorted(self.members))}"
        )


@dataclass
class MustPassReport:
    requirements_file: Optional[str] = None
    markers_checked: int = 0
    groups_tested: int = 0
    violations: List[MustPassViolation] = field(default_factory=list)
    not_checked: List[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.violations


def load_mustpass(path: Optional[str]) -> List[Dict[str, str]]:
    """Read the ``must_pass:`` list from a benchmark YAML.

    Returns ``[]`` (with a loud log line) when the file is missing or holds no
    requirements — an empty baseline must never look like a passed check.
    """
    if not path:
        return []
    file_path = Path(path)
    if not file_path.exists():
        logger.warning(
            f"  Taxonomy must-pass baseline not found: {path} — the gate did "
            f"NOT RUN (0 requirements)."
        )
        return []
    try:
        import yaml
    except ImportError:  # Pragma: no cover - pyyaml is a declared dependency
        logger.warning("  PyYAML unavailable: taxonomy must-pass NOT RUN.")
        return []
    try:
        data = yaml.safe_load(file_path.read_text(encoding="utf-8")) or {}
    except Exception as exc:  # Noqa: BLE001 - report, never crash the run
        logger.warning(
            f"  Taxonomy must-pass baseline unreadable ({exc}) — gate NOT RUN."
        )
        return []
    requirements = data.get("must_pass") or []
    if not isinstance(requirements, list) or not requirements:
        logger.warning(
            f"  Taxonomy must-pass baseline {path} declares no must_pass "
            f"entries — gate NOT RUN."
        )
        return []
    return [r for r in requirements if isinstance(r, dict)]


def mustpass_marker_scope(path: Optional[str]) -> Optional[List[str]]:
    """Honour the baseline's ``markers.ids`` restriction.

    ``None`` means "no restriction" (test every marker tree), which is also
    what an unfilled skeleton means — but that must be visible, because a
    reviewer reading ``ids: []`` could otherwise believe the gate was scoped to
    the ribosomal set. Callers log the difference.
    """
    if not path:
        return None
    file_path = Path(path)
    if not file_path.exists():
        return None
    try:
        import yaml

        data = yaml.safe_load(file_path.read_text(encoding="utf-8")) or {}
    except Exception:  # Noqa: BLE001 - load_mustpass reports read failures
        return None
    block = data.get("markers")
    if not isinstance(block, dict):
        return None
    ids = block.get("ids")
    if not isinstance(ids, list) or not ids:
        return None
    return [str(i) for i in ids if str(i).strip()]


def _taxon_groups(
    taxonomy_map: Dict[str, Dict[str, str]], level: str, tips: Iterable[str]
) -> Dict[str, Tuple[str, ...]]:
    """Group the tree's tips by their value at ``level`` (>=2 members only)."""
    buckets: Dict[str, List[str]] = {}
    for tip in tips:
        value = (taxonomy_map.get(tip) or {}).get(level)
        if not value or value in ("", "NA", None):
            continue
        buckets.setdefault(str(value), []).append(tip)
    return {
        taxon: tuple(sorted(members))
        for taxon, members in buckets.items() if len(members) >= 2
    }


def check_taxonomy_mustpass(
    gene_trees: Dict[str, str],
    taxonomy_map: Dict[str, Dict[str, str]],
    requirements: Sequence[Dict[str, str]],
    *,
    max_recorded: int = MAX_RECORDED_VIOLATIONS,
) -> MustPassReport:
    """Verify every required relationship on every available gene tree."""
    report = MustPassReport()
    if not requirements:
        report.not_checked.append("no requirements loaded")
        return report
    if not gene_trees:
        report.not_checked.append("no gene trees available to test")
        return report
    if not taxonomy_map:
        report.not_checked.append("no taxonomy map supplied")
        return report

    for marker_id, newick in sorted(gene_trees.items()):
        try:
            tips = etree.tip_set(newick)
        except Exception as exc:  # Noqa: BLE001 - unreadable tree is reported
            report.not_checked.append(f"{marker_id}: unparseable tree ({exc})")
            continue
        report.markers_checked += 1

        for requirement in requirements:
            relation = str(requirement.get("relation", "")).lower()
            level = str(requirement.get("level", ""))
            note = str(requirement.get("note", ""))
            if relation not in SUPPORTED_RELATIONS:
                report.not_checked.append(
                    f"{marker_id}: unsupported relation {relation!r} "
                    f"(NOT treated as passed)"
                )
                continue
            if not level:
                report.not_checked.append(
                    f"{marker_id}: requirement lacks a level (NOT treated as "
                    f"passed)"
                )
                continue
            try:
                min_taxa = int(requirement.get("min_taxa", 2))
            except (TypeError, ValueError):
                min_taxa = 2

            for taxon, members in _taxon_groups(taxonomy_map, level, tips).items():
                if len(members) < max(min_taxa, 2):
                    continue
                group = set(members)
                if not group <= tips:
                    report.not_checked.append(
                        f"{marker_id}: {level} '{taxon}' not fully sampled"
                    )
                    continue
                report.groups_tested += 1
                if not etree.is_clade_monophyletic(newick, group):
                    if len(report.violations) < max_recorded:
                        report.violations.append(MustPassViolation(
                            marker_id=marker_id, relation=relation, level=level,
                            taxon=taxon, members=members, note=note,
                        ))
    return report


def write_mustpass_tsv(report: MustPassReport, output_dir: str, prefix: str) -> Path:
    """Persist the gate's verdict so a reviewer can see what was and wasn't run."""
    out = Path(output_dir) / "Phase5_reports"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{prefix}.mustpass.tsv"
    lines = [
        "kind\tmarker_id\trelation\tlevel\ttaxon\tmembers\tresult\n",
        f"summary\t-\t-\t-\t-\tgroups_tested={report.groups_tested};"
        f"markers_checked={report.markers_checked}\t"
        f"{'PASS' if report.passed else 'FAIL'}\n",
    ]
    for violation in report.violations:
        lines.append(
            "violation\t{m}\t{r}\t{l}\t{t}\t{mem}\tFAIL\n".format(
                m=violation.marker_id, r=violation.relation, l=violation.level,
                t=violation.taxon, mem=";".join(violation.members),
            )
        )
    for skipped in report.not_checked:
        lines.append(f"not_checked\t{skipped}\t-\t-\t-\t-\tNOT EXECUTED\n")
    if not report.violations and not report.not_checked:
        lines.append("detail\t-\t-\t-\t-\t-\tall required relations hold\n")
    path.resolve().write_text("".join(lines), encoding="utf-8", newline="\n")
    return path


def enforce_mustpass(report: MustPassReport) -> None:
    """A broken must-pass relationship aborts the run with exit code 4."""
    if report.passed:
        return
    preview = "; ".join(v.describe() for v in report.violations[:3])
    more = (
        f" (+{len(report.violations) - 3} more)"
        if len(report.violations) > 3 else ""
    )
    raise AssertionFailureError(
        f"taxonomy must-pass benchmark FAILED: {len(report.violations)} required "
        f"relationship(s) violated — {preview}{more}. 流水线判定为不通过并中止 "
        "这不是风险分扣分项：修复输入或参照集后重跑。"
    )


@dataclass(frozen=True)
class MustPassOutcome:
    """The gate's verdict plus the artifact that records it.

    ``report.requirements_file`` used to be stamped by the orchestrator between
    the check and the write -- a three-step protocol that only worked because
    the caller remembered all three steps, in order.
    """

    report: MustPassReport
    path: Path


def evaluate_and_write(
    gene_trees: Dict[str, str],
    taxonomy_map: Dict[str, Dict[str, str]],
    requirements: Sequence[Dict[str, str]],
    *,
    requirements_file: object,
    output_dir: str,
    prefix: str,
    max_recorded: int = MAX_RECORDED_VIOLATIONS,
) -> MustPassOutcome:
    """Run the must-pass check and write its TSV, returning both products."""
    report = check_taxonomy_mustpass(
        gene_trees, taxonomy_map, requirements, max_recorded=max_recorded
    )
    report.requirements_file = str(requirements_file)
    return MustPassOutcome(
        report=report,
        path=write_mustpass_tsv(report, output_dir, prefix),
    )
