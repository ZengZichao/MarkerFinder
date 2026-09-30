#!/usr/bin/env python3
"""MarkerFinder end-to-end smoketest.

Exercises the full pipeline --marker-mode hmm with a real TIGRFAM/Pfam HMM
directory, from genome FASTA input through HGT screening to the HTML report.

Usage (synthetic mode — no external data required):
  python3 smoketest.py
  python3 smoketest.py --hmm-dir /path/to/GTDB-Tk-214-Markers

Usage (real data — supply your genome FASTA directory):
  python3 smoketest.py --hmm-dir /path/to/GTDB-Tk-214-Markers \
                       --genomes-dir /path/to/your/genomes/

Markers must be individual.HMM/.hmm profiles below <hmm-dir>/, e.g.:
  GTDB-Tk-214-Markers/tigrfam/individual_hmms/TIGR00006.HMM
  GTDB-Tk-214-Markers/pfam/individual_hmms/PF00380.20.hmm
"""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path


def write_fake_genomes(out_dir: str, n: int = 4) -> str:
    """Write n tiny fake genome FASTA files (proteins only)."""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    import random
    random.seed(42)
    os.makedirs(out_dir, exist_ok=True)
    for i in range(1, n + 1):
        n_proteins = 8 + i
        lines = [">genome_%d" % i]
        for j in range(n_proteins):
            length = 80 + (i * 17 + j * 31) % 120
            seq = "".join(aa[(i * 7 + j * 13 + k) % len(aa)] for k in range(length))
            lines.append(">protein_%d_%d" % (i, j))
            lines.append(seq)
        Path(out_dir, "genome_%d.faa" % i).write_text("\n".join(lines))
    return out_dir


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm-dir", required=True, help="TIGRFAM/Pfam individual_hmms dir")
    ap.add_argument("--genomes-dir", help="explicit genome FASTA dir (else synthetic)")
    ap.add_argument("--out-dir", default="./smoketest_output")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-hmm-score", type=float, default=15.0)
    ap.add_argument("--min-occupancy", type=float, default=0.3)
    args = ap.parse_args()

    if not Path(args.hmm_dir).is_dir():
        print("HMM dir missing:", args.hmm_dir, file=sys.stderr)
        return 2

    work = Path(tempfile.mkdtemp(prefix="mf_smoke_"))
    genomes_dir = args.genomes_dir or write_fake_genomes(str(work / "genomes"))
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== MarkerFinder smoketest ===")
    print("mode             : hmm")
    print("hmm dir          :", args.hmm_dir)
    print("genome dir       :", genomes_dir)
    print("output           :", out_dir)
    print("min-hmm-score    :", args.min_hmm_score)
    print("min-occupancy    :", args.min_occupancy)

    from markerfinder.__main__ import main as mf_main
    argv = [
        "-i", str(genomes_dir),
        "-o", str(out_dir),
        "-t", str(args.threads),
        "--marker-mode", "hmm",
        "--marker-hmm-dir", args.hmm_dir,
        "--min-hmm-score", str(args.min_hmm_score),
        "--min-occupancy", str(args.min_occupancy),
        "--hgt-steps", "phylogenetic",
        "--force",
        "-v",
    ]
    print("argv:", " ".join(argv))
    rc = mf_main(argv)

    print("\n=== result ===")
    print("exit code:", rc)
    outputs = sorted(p.name for p in out_dir.iterdir()) if out_dir.is_dir() else []
    print("output files:", outputs)

    # Cleanup the synthetic genome working copy only when we created it
    if not args.genomes_dir:
        shutil.rmtree(work, ignore_errors=True)

    return rc


if __name__ == "__main__":
    sys.exit(main())
