#!/usr/bin/env python3
"""Build position-known chimera positive controls.

Replaces marker X of genome A with the homologue from distant genome B,
producing an HGT positive whose ground-truth position is known by
construction. Reads/writes only tests/benchmark/downloads (gitignored).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_fasta(path: Path) -> dict:
    seqs: dict = {}
    name = None
    chunks: list = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif line:
            chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def build_chimeras(downloads: Path, chimeras_yaml: Path) -> int:
    """Skeleton builder: swaps are driven by expected/chimera_positives.yaml.

    The real swap requires the per-marker raw FASTA layout produced by the
    GTDB-TK step; this entry point validates the ground truth and the
    downloads manifest so the swap step can never silently no-op.
    """
    try:
        import yaml  # Type: ignore
    except ImportError:
        print("NOT EXECUTED: PyYAML unavailable", file=sys.stderr)
        return 2

    truth = yaml.safe_load(chimeras_yaml.read_text(encoding="utf-8"))
    entries = truth.get("chimeras", [])
    pending = [e for e in entries if str(e.get("marker_id", "")).startswith("REPLACE")]
    manifest = downloads / "manifest.sha256"
    if pending:
        print(
            f"NOT EXECUTED: {len(pending)} chimera entries still carry "
            "REPLACE_* placeholders — fill expected/chimera_positives.yaml "
            "after choosing the order (D-09).",
            file=sys.stderr,
        )
        return 2
    if not manifest.exists():
        print(
            f"NOT EXECUTED: {manifest} missing — run fetch_datasets.sh first.",
            file=sys.stderr,
        )
        return 2
    print("chimera construction ready (implement per-marker swap here)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--downloads", required=True, type=Path)
    args = parser.parse_args()
    here = Path(__file__).parent
    return build_chimeras(args.downloads, here / "expected" / "chimera_positives.yaml")


if __name__ == "__main__":
    sys.exit(main())
