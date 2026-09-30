#!/usr/bin/env python3
""" Performance budget probe.

 caps the *added pure-Python computation* at 5% of total run time for
`n_markers=60`, `n_genomes<=60`. Measuring it end-to-end needs mafft/trimal/
FastTree/IQ-TREE/ASTRAL and real genomes, which this repository deliberately
does not ship. This probe therefore measures the part it CAN measure honestly:

  * the wall-clock cost of the new pure-Python paths (split-set RF, quartet
    topology, PIS) at the shape: 60 markers x 60 taxa;
  * the same cost expressed against a reference external-tool time the user
    can read off their own run log.

It prints the seconds-per-run for the added code and the marker budget, then
tells you exactly which denominator is still missing. It exits non-zero if the
added computation alone exceeds the printed ceiling, so CI can fail loudly on
an accidental O(n^2)-or-worse regression in these paths.

Usage: python tests/benchmark/perf_probe.py [--markers 60] [--taxa 60]
                                             [--repeats 5] [--ceiling 60]
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from markerfinder.utils import etree  # Noqa: E402
from markerfinder.modules.composition import (  # Noqa: E402
    marker_composition_metrics,
)
from markerfinder.utils.informative_sites import (  # Noqa: E402
    parsimony_informative_sites,
)


def random_resolved_newick(tips, rng) -> str:
    """A uniformly random resolved binary rooted Newick string.

    Leaves must be plain strings — the join below unwraps each internal node
    as an ``(left, right)`` pair, which is exactly the bug that made the ete3
    differential tests unrunnable once before.
    """
    clades = [t for t in tips]
    while len(clades) > 2:
        i, j = sorted(rng.sample(range(len(clades)), 2), reverse=True)
        a, b = clades.pop(i), clades.pop(j)
        clades.append((a, b))

    def fmt(node) -> str:
        if isinstance(node, str):
            return node
        left, right = node
        return f"({fmt(left)},{fmt(right)})"

    a, b = clades
    return f"({fmt(a)},{fmt(b)});"


def random_alignment(tips, sites, rng) -> dict:
    residues = "ACDEFGHIKLMNPQRSTVWY"
    return {
        tip: "".join(rng.choice(residues) for _ in range(sites))
        for tip in tips
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--markers", type=int, default=60)
    ap.add_argument("--taxa", type=int, default=60)
    ap.add_argument("--sites", type=int, default=400,
                    help="alignment length per marker (columns)")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--ceiling", type=float, default=60.0,
                    help="seconds of added pure-Python work allowed per run")
    args = ap.parse_args(argv)

    if args.taxa > etree.MAX_PURE_PYTHON_TIPS:
        print(f"taxa {args.taxa} exceeds the pure-Python tip cap "
              f"({etree.MAX_PURE_PYTHON_TIPS}); the code would return NOT_MEASURABLE "
              f"instead of measuring, so the probe would not be measuring the "
              f"budget path.")
        return 2

    rng = random.Random(20260920)  # Fixed seed
    tips = [f"G{k:03d}" for k in range(args.taxa)]

    trees = [random_resolved_newick(tips, rng) for _ in range(args.markers)]
    trees_b = [random_resolved_newick(tips, rng) for _ in range(args.markers)]
    alns = [random_alignment(tips, args.sites, rng) for _ in range(args.markers)]

    best = float("inf")
    for _ in range(args.repeats):
        start = time.perf_counter()
        for i in range(args.markers):
            etree.rf_distance(trees[i], trees_b[i])
            etree.quartet_topology(trees[i], tuple(rng.sample(tips, 4)))
            parsimony_informative_sites(alns[i])
            # 's parallel-evidence path, counted even though the screen is
            # Opt-in: it runs on the same alignment, so it belongs in the budget.
            marker_composition_metrics(alns[i])
        best = min(best, time.perf_counter() - start)

    per_marker = best / args.markers
    print("pure-Python performance budget probe")
    print(f"  shape            : {args.markers} markers x {args.taxa} taxa "
          f"x {args.sites} columns")
    print(f"  repeats          : {args.repeats} (fastest kept)")
    print(f"  added computation: {best:.3f} s per run ({per_marker * 1e3:.2f} ms/marker)")
    print(f"  ceiling          : {args.ceiling:.1f} s per run")
    print(
        "  denominator      : NOT MEASURED — the 5% ratio needs total run time,\n"
        "                     which requires mafft/trimal/FastTree/IQ-TREE/ASTRAL\n"
        "                     and a real genome set. Compare the seconds above\n"
        "                     against your own run's wall time from the log."
    )
    if best > args.ceiling:
        print(f"FAIL: added computation alone exceeds the {args.ceiling}s ceiling.")
        return 1
    print("PASS: added computation is within the absolute ceiling.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
