#!/usr/bin/env python3
"""Reverse ablation — compare rankings that differ in ONE factor.

 asks for "同一份数据只改一个因子（有无 HGT 过滤 / 串联 vs 合并），比较排序变化".
The discipline in that sentence is the whole point: a ranking difference is only
evidence about the factor you changed. Two runs that differ in a filter setting
*and* in the tree-building method say nothing attributable to either, which is how
"785 tests passed but the conclusion was wrong" survives a green suite.

So the entry point here is not "compute a correlation" — it is *refuse to compute
one* unless exactly one declared factor differs, and refuse in both directions:

* zero factors differ -> the runs are the same experiment; there is nothing to
  attribute (a real risk when a switch silently fails to take effect, which this
  three times: `--scan-stability-min`, `--min-informative-sites`,
  `composition`).
* two or more differ -> the comparison is unattributable, so the report is
  REFUSED rather than printed with a caveat.

Ranking disagreement is measured with Kendall's tau-b (tie-corrected, so a block
of markers sharing a support value does not inflate the coefficient) plus the
per-marker rank displacement, and written as a product next to the two runs.
Pure stdlib; deterministic; nothing here touches an external tool.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


class AblationRefusal(Exception):
    """Raised when the two runs may not be compared ('s one-factor rule)."""


@dataclass
class AblationReport:
    factor: str
    tau: float
    n_shared: int
    displacements: List[Tuple[str, int, int, int]] = field(default_factory=list)
    only_in_a: List[str] = field(default_factory=list)
    only_in_b: List[str] = field(default_factory=list)

    @property
    def max_displacement(self) -> int:
        return max((abs(d[3]) for d in self.displacements), default=0)

    @property
    def moved(self) -> List[Tuple[str, int, int, int]]:
        """Markers whose rank changed, biggest move first (ties by marker id)."""
        return sorted(
            [row for row in self.displacements if row[3]],
            key=lambda row: (-abs(row[3]), row[0]),
        )


def differing_factors(a: Dict[str, str], b: Dict[str, str]) -> List[str]:
    """Declared config keys whose values differ (a key missing on one side counts)."""
    keys = sorted(set(a) | set(b))
    return [k for k in keys if a.get(k, "<absent>") != b.get(k, "<absent>")]


def kendall_tau_b(x: Sequence[float], y: Sequence[float]) -> float:
    """Tie-corrected rank correlation: (C - D) / sqrt((C+D+Tx)(C+D+Ty)).

    Implemented directly from the definition (n is a marker count, so O(n^2) is
    the right trade for having no scipy dependency — ).
    """
    n = len(x)
    if n < 2:
        raise ValueError("tau needs at least two paired observations")
    concordant = discordant = ties_x = ties_y = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            if dx == 0 and dy == 0:
                continue
            if dx == 0:
                ties_x += 1
            elif dy == 0:
                ties_y += 1
            elif (dx > 0) == (dy > 0):
                concordant += 1
            else:
                discordant += 1
    denom = math.sqrt(
        (concordant + discordant + ties_x) * (concordant + discordant + ties_y)
    )
    if denom == 0.0:
        # Every pair is tied on at least one side: the ordering carries no
        # Information, which is not the same statement as "perfectly correlated".
        raise ValueError("all observations tied — tau undefined, not 1.0")
    return (concordant - discordant) / denom


def rank_by_score(scores: Dict[str, float]) -> Dict[str, int]:
    """Dense ranks, best (highest) first; ties broken by marker id.

    The id tie-break matters for reproducibility: ranking ties by dict
    insertion order would make two runs of the same data disagree.
    """
    order = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    ranks: Dict[str, int] = {}
    previous: Optional[float] = None
    current = 0
    for pos, (marker, score) in enumerate(order, start=1):
        if previous is None or score != previous:
            current = pos
            previous = score
        ranks[marker] = current
    return ranks


def compare_runs(
    ranks_a: Dict[str, int],
    ranks_b: Dict[str, int],
    factors_a: Dict[str, str],
    factors_b: Dict[str, str],
    *,
    min_shared: int = 3,
) -> AblationReport:
    """The one-factor gate plus the measurement. Raises AblationRefusal."""
    diff = differing_factors(factors_a, factors_b)
    if not diff:
        raise AblationRefusal(
            "the two runs declare identical settings, so there is nothing to "
            "attribute. Either the ablation switch did not take effect (this "
            "three dead switches of that kind) or the same "
            "experiment was run twice."
        )
    if len(diff) > 1:
        raise AblationRefusal(
            f"{len(diff)} factors differ ({', '.join(diff)}); a ranking "
            "difference cannot be attributed to any one of them. Re-run with "
            "exactly one change."
        )
    shared = sorted(set(ranks_a) & set(ranks_b))
    if len(shared) < min_shared:
        raise AblationRefusal(
            f"only {len(shared)} marker(s) are present in both runs "
            f"(need >= {min_shared}); a 'ranking' of that is not a ranking."
        )
    xa = [float(ranks_a[m]) for m in shared]
    xb = [float(ranks_b[m]) for m in shared]
    return AblationReport(
        factor=diff[0],
        tau=kendall_tau_b(xa, xb),
        n_shared=len(shared),
        displacements=[
            (m, ranks_a[m], ranks_b[m], ranks_b[m] - ranks_a[m]) for m in shared
        ],
        only_in_a=sorted(set(ranks_a) - set(ranks_b)),
        only_in_b=sorted(set(ranks_b) - set(ranks_a)),
    )


def read_ranking(path: Path, id_column: str, score_column: str) -> Dict[str, float]:
    """Load ``marker -> score`` from a TSV produced by a run."""
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise AblationRefusal(f"{path}: no data rows")
    missing = [c for c in (id_column, score_column) if c not in rows[0]]
    if missing:
        raise AblationRefusal(
            f"{path}: column(s) {missing} absent; got {sorted(rows[0])}"
        )
    out: Dict[str, float] = {}
    unparsable: List[str] = []
    for row in rows:
        marker = (row[id_column] or "").strip()
        if not marker:
            continue
        try:
            out[marker] = float(str(row[score_column]).strip())
        except ValueError:
            unparsable.append(marker)  # NA / UNKNOWN is not a score
    if unparsable:
        print(
            f"  NOT_MEASURABLE: {len(unparsable)} marker(s) in {path.name} carry no "
            f"numeric {score_column} and are excluded from the ranking "
            f"({', '.join(sorted(unparsable)[:10])}{' ...' if len(unparsable) > 10 else ''})",
            file=sys.stderr,
        )
    if not out:
        raise AblationRefusal(f"{path}: no parsable numeric {score_column} values")
    return out


def write_ablation_tsv(report: AblationReport, out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Resolve pins the target under a real directory before any write, the
    # Same convention the production tree adopted when the scanner flagged
    # Dynamic output paths.
    target = (out_dir / f"{prefix}.ablation.tsv").resolve()
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["marker_id", "rank_baseline", "rank_ablated", "displacement"])
        for marker, ra, rb, delta in sorted(report.displacements):
            writer.writerow([marker, ra, rb, delta])
    return target


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-a", type=Path, required=True, help="baseline ranking TSV")
    parser.add_argument("--run-b", type=Path, required=True, help="ablated ranking TSV")
    parser.add_argument("--id-column", default="marker_id")
    parser.add_argument("--score-column", default="overall_score")
    parser.add_argument(
        "--factor", action="append", default=[], metavar="NAME=OLD:NEW",
        help="the single changed factor, e.g. hgt_filter=on:off (repeatable, "
             "but more than one is a refusal)",
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--prefix", default="ablation")
    args = parser.parse_args(argv)

    if not args.factor:
        print(
            "NOT EXECUTED: no --factor NAME=OLD:NEW given, so the run "
            "cannot state which single factor changed. A ranking difference "
            "without a named cause is not an ablation.",
            file=sys.stderr,
        )
        return 2
    factors_a: Dict[str, str] = {}
    factors_b: Dict[str, str] = {}
    for spec in args.factor:
        name, _, pair = spec.partition("=")
        old, sep, new = pair.partition(":")
        if not sep or not name.strip():
            print(
                f"NOT EXECUTED: --factor {spec!r} is not NAME=OLD:NEW.",
                file=sys.stderr,
            )
            return 2
        factors_a[name.strip()] = old.strip()
        factors_b[name.strip()] = new.strip()

    try:
        report = compare_runs(
            rank_by_score(read_ranking(args.run_a, args.id_column, args.score_column)),
            rank_by_score(read_ranking(args.run_b, args.id_column, args.score_column)),
            factors_a,
            factors_b,
        )
    except (AblationRefusal, OSError) as exc:
        print(f"NOT EXECUTED: {exc}", file=sys.stderr)
        return 2

    print(
        f"reverse ablation — factor {report.factor}: "
        f"tau_b={report.tau:+.4f} over {report.n_shared} shared marker(s); "
        f"max |displacement|={report.max_displacement}; "
        f"{len(report.moved)} marker(s) moved"
    )
    for marker, ra, rb, delta in report.moved[:20]:
        print(f"  {marker}: {ra} -> {rb} ({delta:+d})")
    if report.only_in_a or report.only_in_b:
        print(
            f"  present in only one run: a-only={len(report.only_in_a)}, "
            f"b-only={len(report.only_in_b)} (excluded from tau)"
        )
    if args.out_dir is not None:
        print(f"  written: {write_ablation_tsv(report, args.out_dir, args.prefix)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
