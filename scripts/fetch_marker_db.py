#!/usr/bin/env python3
"""Obtain the per-marker TIGRFAM/Pfam HMM profiles used by ``--marker-mode hmm``.

Why this is a fetch step and not a directory in the repository
-------------------------------------------------------------
``db/gtdb_markers/{bac120,ar53}`` holds one HMM profile per marker. Those
profiles are third-party work — TIGRFAM models from NCBI, Pfam models from
InterPro, repackaged by GTDB-Tk — so they are not redistributed here. This
script assembles the layout MarkerFinder discovers from a source you point it
at, and records what it did in ``db/gtdb_markers/README.md`` so a later reader
can tell which upstream files produced the profiles on disk.

What it accepts
---------------
``--from-hmm FILE`` a *concatenated* HMM library, i.e. one file holding
                           many models (GTDB-Tk ships
                           ``gtdbtk_bac120.a.hmm`` / ``gtdbtk_ar53.a.hmm``).
                           Each model is written out as its own file, named by
                           its ``ACC`` line.
``--from-dir DIR`` a directory that already holds one file per marker
                           (``*.hmm`` / ``*.HMM``), e.g. a GTDB-Tk
                           ``marker_genes`` dump or an unpacked archive.
``--auto`` look for a GTDB-Tk installation on this machine
                           (``share/gtdbtk*/marker_databases``, ``<prefix>/marker_databases``)
                           and split what is found.

``--set`` selects which marker set the input fills; the two sets are what
``marker_db_source: auto`` searches for by domain (``bac120`` Bacteria, ``ar53``
Archaea).

Integrity
---------
After writing, the script prints the combined SHA-256 of the directory, computed
exactly the way ``markerfinder --check`` computes it. ``--record`` writes that
value into ``db/expected_hashes.json``, so the profiles on disk become the
reference the next run is checked against; leave it out to keep the shipped
reference untouched and only compare later with ``markerfinder --check``.

Usage
-----
    python3 scripts/fetch_marker_db.py --auto
    python3 scripts/fetch_marker_db.py --set bac120 --from-hmm ~/gtdbtk/marker_databases/gtdbtk_bac120.a.hmm
    python3 scripts/fetch_marker_db.py --set ar53 --from-dir /path/to/ar53_marker_genes
    python3 scripts/fetch_marker_db.py --auto --record

Exit codes: 0 wrote (or verified) something, 1 nothing usable found, 2 usage error.

No third-party imports: stdlib only.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parents[1]
DB_DIR = REPO / "db"
DEFAULT_OUT = DB_DIR / "gtdb_markers"
HASHES_FILE = DB_DIR / "expected_hashes.json"

SETS = ("bac120", "ar53")

# A concatenated HMM library starts each model with its version header line.
MODEL_START = re.compile(r"^HMMER3/[fbbtrs] \[")
NAME_LINE = re.compile(r"^NAME\s+(.+)$")
ACC_LINE = re.compile(r"^ACC\s+(\S+)")
# TIGRFAM accessions are TIGR + 5 digits; Pfam are PF + 5 digits +.version.
TIGR_ID = re.compile(r"^TIGR\d{5}$")
PFAM_ID = re.compile(r"^PF\d{5}\.\d+$")


def _profile_suffix(stem: str) -> str:
    """``.hmm`` for Pfam, ``.HMM`` otherwise — the layout the scanner expects."""
    return ".hmm" if stem.startswith("PF") else ".HMM"


def split_concatenated(hmm_file: Path) -> List[Tuple[str, str]]:
    """Return ``(marker_id, model_text)`` for every model in a HMM library.

    The model text keeps its header line through the closing ``//``, which is
    what HMMER needs to read a standalone profile back.
    """
    models: List[Tuple[str, str]] = []
    block: List[str] = []
    in_model = False

    def flush(lines: List[str]) -> None:
        acc = name = ""
        for line in lines:
            match = ACC_LINE.match(line)
            if match and not acc:
                acc = match.group(1)
                continue
            match = NAME_LINE.match(line)
            if match and not name:
                name = match.group(1).strip()
        marker_id = acc or name
        if not marker_id:
            return
        models.append((marker_id, "".join(lines)))

    with hmm_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if MODEL_START.match(line):
                if in_model and block:
                    flush(block)
                block = [line]
                in_model = True
                continue
            if in_model:
                block.append(line)
                if line.rstrip() == "//":
                    flush(block)
                    block, in_model = [], False
    if in_model and block:
        flush(block)
    return models


def marker_id_from_file(path: Path) -> Optional[str]:
    """The marker a standalone profile belongs to: its ``ACC`` line, else filename."""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("//"):
                    break
                if line.startswith("HMMER"):
                    continue
                match = ACC_LINE.match(line)
                if match:
                    acc = match.group(1).split()[0]
                    if PFAM_ID.match(acc) or TIGR_ID.match(acc):
                        return acc
                match = NAME_LINE.match(line)
                if match:
                    name = match.group(1).strip()
                    if PFAM_ID.match(name) or TIGR_ID.match(name):
                        return name
    except OSError:
        return None
    stem = path.stem
    return stem if (PFAM_ID.match(stem) or TIGR_ID.match(stem)) else None


def gtdbtk_candidates() -> Dict[str, Path]:
    """Best-effort location of a GTDB-Tk marker database on this machine."""
    roots: List[Path] = []
    conda_prefix = Path(sys.prefix)
    for candidate in (
        conda_prefix / "share",
        conda_prefix / "marker_databases",
        Path.home() / ".conda" / "envs",
    ):
        if candidate.is_dir():
            roots.append(candidate)
    found: Dict[str, Path] = {}
    for root in roots:
        for hmm in sorted(root.rglob("gtdbtk_*.hmm")):
            for marker_set in SETS:
                if marker_set in hmm.name and marker_set not in found:
                    found[marker_set] = hmm
    return found


def combined_hash(target: Path) -> str:
    """Directory hash, computed the way ``DatabaseVersionManager.verify_database`` does."""
    digests = [
        hashlib.sha256(f.read_bytes()).hexdigest()
        for f in sorted(target.rglob("*")) if f.is_file()
    ]
    return hashlib.sha256("".join(digests).encode()).hexdigest()


def write_profiles(models: Dict[str, str], target: Path, overwrite: bool) -> int:
    target.mkdir(parents=True, exist_ok=True)
    written = 0
    for marker_id, text in sorted(models.items()):
        path = target / f"{marker_id}{_profile_suffix(marker_id)}"
        if path.exists() and not overwrite:
            continue
        path.write_text(text, encoding="utf-8")
        written += 1
    return written


def record_provenance(out: Path, notes: List[str]) -> None:
    """Write the mapping/provenance file ``DatabaseVersionManager`` reports."""
    lines = [
        "# GTDB-Tk marker HMM profiles (fetched, not redistributed)",
        "",
        f"- assembled by: `scripts/fetch_marker_db.py` on "
        f"{_dt.date.today().isoformat()}",
        "- contents: one HMM profile per marker; filename stem = marker "
        "accession",
        "  (`TIGR#####.HMM` for TIGRFAM, `PF#####.##.hmm` for Pfam)",
        "- upstream: TIGRFAM models (NCBI) and Pfam models (InterPro), "
        "redistributed by GTDB-Tk under its release terms",
        "- this directory is git-ignored on purpose; see `db/README.md` for how "
        "to rebuild it",
        "",
        "## This build",
        "",
    ]
    lines += [f"- {note}" for note in notes]
    lines.append("")
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")


def record_hash(marker_set: str, digest: str) -> None:
    """Store the directory hash in ``db/expected_hashes.json`` for ``--check``."""
    data: Dict[str, object] = {}
    if HASHES_FILE.exists():
        try:
            data = json.loads(HASHES_FILE.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            print(f"  {HASHES_FILE} is not valid JSON; not recording", file=sys.stderr)
            return
    data[marker_set] = digest
    data["_recorded"] = (
        f"{_dt.date.today().isoformat()}: {marker_set} recomputed by "
        f"scripts/fetch_marker_db.py from the profiles on disk"
    )
    HASHES_FILE.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"  recorded {marker_set} hash in {HASHES_FILE.relative_to(REPO)}")


def collect(models: Dict[str, str], source: str, notes: List[str]) -> int:
    """Drop entries with an unusable id, log the source, return the accepted count."""
    accepted = {
        marker_id: text for marker_id, text in models.items()
        if PFAM_ID.match(marker_id) or TIGR_ID.match(marker_id)
    }
    skipped = len(models) - len(accepted)
    models.clear()
    models.update(accepted)
    notes.append(f"source `{source}`: {len(accepted)} profiles"
                 + (f", {skipped} unidentifiable models skipped" if skipped else ""))
    return len(accepted)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assemble db/gtdb_markers from a GTDB-Tk marker database.")
    parser.add_argument("--set", choices=SETS, help="marker set to fill")
    parser.add_argument("--from-hmm", type=Path,
                        help="concatenated HMM library to split into per-marker files")
    parser.add_argument("--from-dir", type=Path,
                        help="directory that already holds one HMM file per marker")
    parser.add_argument("--auto", action="store_true",
                        help="find and split a GTDB-Tk installation on this machine")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT,
                        help=f"target directory (default: {DEFAULT_OUT.relative_to(REPO)})")
    parser.add_argument("--force", action="store_true",
                        help="overwrite profiles that are already on disk")
    parser.add_argument("--record", action="store_true",
                        help="write the resulting directory hash into "
                             "db/expected_hashes.json")
    args = parser.parse_args(argv)

    if args.auto:
        found = gtdbtk_candidates()
        if not found:
            print("no GTDB-Tk marker database found; use --from-hmm or --from-dir",
                  file=sys.stderr)
            return 1
        plans = [(marker_set, {"hmm": path}) for marker_set, path in found.items()]
    elif args.set and (args.from_hmm or args.from_dir):
        plans = [(args.set, {"hmm": args.from_hmm, "dir": args.from_dir})]
    else:
        print("give --auto, or --set together with --from-hmm/--from-dir",
              file=sys.stderr)
        return 2

    notes: List[str] = []
    filled: List[str] = []
    for marker_set, source in plans:
        models: Dict[str, str] = {}
        hmm_path: Optional[Path] = source.get("hmm")
        dir_path: Optional[Path] = source.get("dir")
        if hmm_path:
            if not hmm_path.is_file():
                print(f"  {hmm_path} is not a file", file=sys.stderr)
                continue
            for marker_id, text in split_concatenated(hmm_path):
                models[marker_id] = text
            if not collect(models, str(hmm_path), notes):
                continue
        elif dir_path:
            if not dir_path.is_dir():
                print(f"  {dir_path} is not a directory", file=sys.stderr)
                continue
            for path in sorted(dir_path.iterdir()):
                if path.suffix.lower() not in (".hmm", ".gz") or not path.is_file():
                    continue
                marker_id = marker_id_from_file(path)
                if marker_id:
                    models[marker_id] = path.read_text(encoding="utf-8",
                                                       errors="replace")
            if not collect(models, str(dir_path), notes):
                continue
        target = args.out / marker_set
        written = write_profiles(models, target, args.force)
        total = len(list(target.iterdir()))
        print(f"{marker_set}: {written} profile(s) written to {target}, "
              f"{total} on disk, combined sha256 {combined_hash(target)}")
        notes.append(f"`{marker_set}`: {total} profiles on disk at "
                     f"{target.relative_to(REPO) if target.is_relative_to(REPO) else target}")
        filled.append(marker_set)
        if args.record:
            record_hash(marker_set, combined_hash(target))

    if not filled:
        return 1
    record_provenance(args.out, notes)
    print(f"\nMarkerFinder will now find these profiles under {args.out} "
          "(marker_db_source: auto), or pass --marker-hmm-dir.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
