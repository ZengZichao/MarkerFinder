#!/usr/bin/env python3
"""One-command entry point for the MarkerFinder validation suite.

    python validation/run_validation.py # everything but the slow benchmark
    python validation/run_validation.py --all # including the 12-genome detection case
    python validation/run_validation.py --prepare # rebuild the test data first
    python validation/run_validation.py -n 8 # parallel workers (needs pytest-xdist)
    python validation/run_validation.py -k v06 # one case file only

Why a wrapper rather than ``pytest`` alone: the suite needs three things in a
fixed order — data present, environment sane, results archived — and each of
them is a failure mode that produces a confusing pytest error if it is skipped.
The wrapper checks them and says what to run.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
REPO = HERE.parent
CASES = HERE / "cases"
RESULTS = HERE / "results"
WORK = HERE / ".work"

REQUIRED_TOOLS = ("hmmsearch", "mafft", "trimal", "FastTree", "astral")


def _anonymise_junit(report: Path) -> None:
    """Keep the archived run record free of host identity.

    pytest writes the machine hostname into the report and, when a case fails,
    absolute paths of this checkout into the failure text. Neither belongs in a
    published artifact: a run is attributed by ``environment.txt``, not by the
    name of whoever's workstation produced it.
    """
    if not report.exists():
        return
    text = report.read_text(encoding="utf-8")
    text = re.sub(r'hostname="[^"]*"', 'hostname="validation-host"', text)
    text = text.replace(str(REPO), "<repo>").replace(str(Path.home()), "~")
    report.write_text(text, encoding="utf-8")


def _preflight(check_data: bool) -> int:
    missing_tools = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing_tools:
        print("ERROR: external tools missing from PATH: " + ", ".join(missing_tools),
              file=sys.stderr)
        print("The suite runs the real pipeline; it does not mock the tools. "
              "Activate the MarkerFinder environment (conda activate "
              "markerfinder) or install them: conda env create -f "
              "environment.yml", file=sys.stderr)
        return 2

    if check_data and not (HERE / "data" / "MANIFEST.sha256").exists():
        print("ERROR: validation/data is not built.", file=sys.stderr)
        print("Run: python validation/scripts/01_fetch_genomes.py && "
              "python validation/scripts/02_build_marker_sets.py && "
              "python validation/scripts/03_build_fixtures.py", file=sys.stderr)
        return 2
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="include the slow 12-genome planted-chimera detection case")
    ap.add_argument("--prepare", action="store_true",
                    help="run the three data-build scripts before the cases")
    ap.add_argument("-n", "--jobs", type=int, default=1,
                    help="parallel pytest workers (requires pytest-xdist)")
    ap.add_argument("-k", dest="select", default=None,
                    help="pytest expression selecting cases by name")
    ap.add_argument("--exitfirst", "-x", action="store_true",
                    help="stop at the first failure")
    args = ap.parse_args(argv)

    if args.prepare:
        for script in ("01_fetch_genomes.py", "02_build_marker_sets.py",
                       "03_build_fixtures.py"):
            cmd = [sys.executable, str(HERE / "scripts" / script)]
            print("$", " ".join(cmd))
            rc = subprocess.run(cmd, cwd=str(HERE / "scripts")).returncode
            if rc != 0:
                print(f"ERROR: {script} failed with exit code {rc}", file=sys.stderr)
                return rc

    rc = _preflight(check_data=not args.prepare)
    if rc != 0:
        return rc

    WORK.mkdir(parents=True, exist_ok=True)
    # Two scratch directories are placed under.work on purpose.
    #
    # TMPDIR: MAFFT builds its staging directories with ``mktemp -dt`` and, when
    # That fails, keeps running with an empty path — a non-writable /tmp therefore
    # Surfaces as "no marker could be aligned", i.e. a wrong scientific conclusion
    # Rather than an error.
    #
    # --basetemp: pytest's tmp_path root falls back to ``pytest-of-<user>/`` in the
    # Current directory when the system temp dir is unusable, which litters the
    # Repository root.
    env = dict(os.environ, TMPDIR=str(WORK), MAFFT_TMPDIR=str(WORK))

    RESULTS.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "pytest", str(CASES), "-v", "--tb=short",
           "-p", "no:cacheprovider", f"--basetemp={WORK / 'pytest-basetemp'}"]
    if not args.all:
        cmd += ["-m", "not slow"]
    if args.jobs > 1:
        cmd += ["-n", str(args.jobs)]
    if args.select:
        cmd += ["-k", args.select]
    if args.exitfirst:
        cmd.append("-x")
    # A machine-readable record of the run, kept beside the archived metrics.
    report = RESULTS / "junit.xml"
    cmd += [f"--junit-xml={report}"]

    print("$", " ".join(cmd))
    print(f"(cwd={REPO})")
    rc = subprocess.run(cmd, cwd=str(REPO), env=env).returncode
    _anonymise_junit(report)
    print()
    # Archive what the machine actually looked like for this run, so the shipped
    # Results can be attributed to a verifiable environment.
    subprocess.run([sys.executable, str(HERE / "scripts" / "04_record_environment.py"),
                    "--out", str(RESULTS / "environment.txt")],
                   cwd=str(REPO), env=env)
    print(f"junit report: {report}")
    print(f"metrics:      {RESULTS / 'metrics.json'}")
    print(f"coverage:     {RESULTS / 'capability_matrix.tsv'}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
