#!/usr/bin/env python3
"""Step 4: record the machine-readable environment the published run used.

A test result is only interpretable alongside the environment that produced it,
so this writes ``validation/results/environment.txt``: the interpreter, the
declared runtime dependencies with their installed versions, and the version
each external tool reports for itself.

Nothing here is transcribed from the documentation. Every line is executed
output, because the documentation is what is being audited.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent
RESULTS = ROOT / "results"
# Invoked as ``python validation/scripts/04_record_environment.py``, the
# Interpreter puts *this* directory on sys.path, not the repository root, so the
# Package would be unimportable without an editable install.
if str(ROOT.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent))

# (binary, argv that makes it print its own version)
#
# Hmmsearch has no ``--version``; its ``-h`` banner carries the release line.
# FastTree has no ``--version`` either: invoked with no arguments it prints
# Usage and exits 1, and with an unknown option it complains — neither names
# The release, so ``FastTree -help`` is captured instead.
TOOLS = [
    ("hmmsearch", ["hmmsearch", "-h"]),
    ("mafft", ["mafft", "--version"]),
    ("trimal", ["trimal", "--version"]),
    ("FastTree", ["FastTree", "-help"]),
    ("iqtree3", ["iqtree3", "--version"]),
    ("astral", ["astral", "--version"]),
    ("checkm", ["checkm", "-v"]),
    ("datasets", ["datasets", "--version"]),
]

# A version-looking line: at least one "digits.digits" run.
_VERSION_RE = re.compile(r"\d+\.\d+")


def capture(cmd: list) -> str:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except FileNotFoundError:
        return "(not installed)"
    except subprocess.TimeoutExpired:
        return "(timed out)"
    text = (proc.stdout or proc.stderr or "").strip()
    if not text:
        return f"(no output, exit {proc.returncode})"
    # Prefer a line that actually names a release, so the record reports the
    # Tool's own version string instead of its usage banner.
    for line in text.splitlines():
        if _VERSION_RE.search(line):
            return line.strip().lstrip("#").strip()
    return text.splitlines()[0].strip()


def python_packages() -> list:
    from importlib import metadata

    from markerfinder import banner

    names = banner.runtime_dependencies()
    rows = []
    for name in names:
        try:
            version = metadata.version(name)
        except metadata.PackageNotFoundError:
            version = "(not installed)"
        rows.append((name, version))
    return rows


def _redact(text: str) -> str:
    """Replace machine-specific prefixes with portable placeholders.

    The record ships with the release; the account that produced it is not part
    of the evidence. ``<repo>``, ``<env>`` and ``~`` keep each line meaningful on
    any other machine.
    """
    text = text.replace(str(ROOT.parent), "<repo>")
    text = text.replace(str(Path.home()), "~")
    return re.sub(r"~[/\\][^\s,;)]*?[/\\]envs[/\\][A-Za-z0-9_.-]+", "<env>", text)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(RESULTS / "environment.txt"),
                    help="where to write the record (default: "
                         "validation/results/environment.txt)")
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "MarkerFinder validation suite - environment record",
        f"recorded_utc: {_dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds')}",
        "",
        "[host]",
        f"platform        : {platform.platform()}",
        f"python          : {sys.version.split()[0]} ({sys.executable})",
        f"implementation  : {platform.python_implementation()}",
        f"machine         : {platform.machine()}",
        f"cpu_count       : {os.cpu_count()}",
        "",
        "[markerfinder]",
    ]

    import markerfinder

    lines.append(f"version         : {markerfinder.__version__}")
    lines.append(f"imported_from   : {Path(markerfinder.__file__).parent}")

    lines += ["", "[declared runtime dependencies (from the installed distribution)]"]
    for name, version in python_packages():
        lines.append(f"{name:<20}: {version}")

    lines += ["", "[external tools on PATH]"]
    for name, cmd in TOOLS:
        resolved = shutil.which(name)
        version = capture(cmd)
        lines.append(f"{name:<20}: {version}")
        if resolved:
            lines.append(f"{'':<20} path: {resolved}")
    lines.append("")
    lines.append(f"TMPDIR          : {os.environ.get('TMPDIR', '(unset)')}")
    lines.append(f"MAFFT_TMPDIR    : {os.environ.get('MAFFT_TMPDIR', '(unset)')}")

    path = out
    path.write_text("\n".join(_redact(line) for line in lines) + "\n", encoding="utf-8")
    print(f"[write] {path}")
    print("\n".join(_redact(line) for line in lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
