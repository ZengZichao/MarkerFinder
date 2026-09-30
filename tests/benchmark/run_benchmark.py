#!/usr/bin/env python3
"""Run the external benchmark end to end.

Pipeline -> decision cards -> precision/recall against the planted-chimera
truth. MISSING TOOLS / MISSING DATA print "NOT EXECUTED" per stage and exit
non-zero: a benchmark that silently passes without tools reproduces 's
forbidden failure mode.

User-triggered: not executed by any automated workflow.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Optional, Sequence

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))  # Repo root for markerfinder import

REQUIRED_TOOLS = ("mafft", "trimal", "FastTree", "astral", "fasttree")


def _tools_available() -> dict:
    available = {}
    for tool in REQUIRED_TOOLS:
        available[tool] = shutil.which(tool) is not None
    return available


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--downloads", type=Path, default=HERE / "downloads")
    parser.add_argument(
        "--ablation-a", type=Path, default=None,
        help="ranking TSV of the baseline run (marker_summary.tsv)",
    )
    parser.add_argument(
        "--ablation-b", type=Path, default=None,
        help="ranking TSV of the run that changes exactly ONE factor",
    )
    parser.add_argument(
        "--ablation-factor", action="append", default=[], metavar="NAME=OLD:NEW",
        help="the single changed factor; more than one is a refusal",
    )
    args = parser.parse_args(argv)

    # The reverse-ablation comparator is pure logic over two existing run
    # Products, so it runs BEFORE the tool probe — a user who already has the two
    # Runs should not need the aligners installed to get its answer. Naming the
    # Changed factor is mandatory: an unattributed ranking difference is not an
    # Ablation, and more than one changed factor is a refusal (see ablation.py).
    if args.ablation_a or args.ablation_b or args.ablation_factor:
        from ablation import main as run_ablation

        missing = [
            label for label, value in (
                ("--ablation-a", args.ablation_a),
                ("--ablation-b", args.ablation_b),
            )
            if value is None
        ]
        if missing:
            print(
                f"NOT EXECUTED: {', '.join(missing)} is required: an "
                "ablation compares TWO runs differing in exactly one factor.",
                file=sys.stderr,
            )
            return 2
        rc = run_ablation(
            ["--run-a", str(args.ablation_a), "--run-b", str(args.ablation_b)]
            + [item for factor in args.ablation_factor for item in ("--factor", factor)]
        )
        if rc != 0:
            return rc

    tools = _tools_available()
    missing = [t for t, ok in tools.items() if not ok]
    if missing:
        print(
            f"NOT EXECUTED: external tools missing on PATH: {missing}. "
            "Install (conda env create -f environment.yml) and re-run — "
            "the benchmark refuses to silently pass.",
            file=sys.stderr,
        )
        return 2

    if not (args.downloads / "manifest.sha256").exists():
        print(
            "NOT EXECUTED: downloads/manifest.sha256 missing — run "
            "fetch_datasets.sh first.",
            file=sys.stderr,
        )
        return 2

    # The set must be reproducible from a recorded provenance
    # (registry, DOI, licence, retrieval date, checksum). Blank or placeholder
    # Fields stop the run before any metric is produced.
    from provenance import load_accessions, missing_provenance, summarize

    accessions = load_accessions()
    absent = missing_provenance(accessions)
    if absent:
        print(
            "NOT EXECUTED: provenance incomplete — missing or "
            f"placeholder field(s): {', '.join(absent)}.\n"
            "        Fill `provenance:` in tests/benchmark/expected/"
            "accessions.yaml (source_registry, doi, license, retrieved_utc, "
            "checksum_sha256) before trusting any number from this benchmark.",
            file=sys.stderr,
        )
        return 2
    print("Provenance recorded for this benchmark run:")
    summarize(accessions)

    # Real stages (chimera build -> pipeline -> decision cards -> metrics) are
    # Implemented against the materialised downloads; with placeholder
    # Accessions the benchmark still refuses to fabricate numbers.
    print(
        "NOT EXECUTED: expected/accessions.yaml still carries placeholders "
        "(D-09 leaves the order choice to the user). Fill it, re-fetch, then "
        "re-run.",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    sys.exit(main())
