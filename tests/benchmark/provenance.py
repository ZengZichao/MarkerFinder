#!/usr/bin/env python3
"""External data must be reproducible from a record.

 says external data and products stay out of the repository, replaced by
"a DOI + a reproducible script", and that "GTDB version must be cited
explicitly". lists the same thing as an acceptance item: 来源 DOI / 许可证 /
获取日期 / 校验和记录在案.

Until this module existed there was nowhere in the shipped skeleton to record
any of it, so a user who filled in accessions and ran the benchmark would end
up with a result nobody else can reproduce — the exact opposite of what
is for. Here the fields are declared, and the runners refuse to start without
them.

Deliberately strict: an empty string, ``None``, or a leftover sentinel such as
``REPLACE_WITH_ORDER_NAME`` / ``TBD`` counts as MISSING. A benchmark that runs
on vague provenance reports a number nobody can re-derive.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

HERE = Path(__file__).resolve().parent
DEFAULT_ACCESSIONS = HERE / "expected" / "accessions.yaml"

# Keys required under `provenance:` in the accessions file, plus the GTDB
# Release the set was cut against ('s "GTDB 版本须显式引用").
REQUIRED_FIELDS = (
    "source_registry",
    "doi",
    "license",
    "retrieved_utc",
    "checksum_sha256",
    "gtdb_release",
)

# Values that mean "nobody filled this in yet" but would otherwise look present.
# NOTE: the empty string is handled by the blank check in _is_filled and must
# NOT live here — "" is a substring of every string, so including it made every
# Value read as missing and no record could ever validate (caught by test).
SENTINELS = ("tbd", "todo", "n/a", "na", "none", "placeholder", "replace")


def sha256_file(path) -> str:
    """Checksum a downloaded artefact, streaming (genomes are large)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_filled(value) -> bool:
    if value is None:
        return False
    text = (value if isinstance(value, str) else str(value)).strip()
    if not text:
        return False
    # Exact-match sentinels only: substring matching would reject legitimate
    # Words that merely contain them (e.g. "Nanopore" contains "na").
    if text.lower() in SENTINELS:
        return False
    if "REPLACE_WITH" in text or (text.startswith("<") and text.endswith(">")):
        return False
    return True


def load_accessions(path: Optional[Path] = None) -> Dict:
    """Read the accessions/provenance YAML (returns {} when unreadable)."""
    import yaml

    target = Path(path) if path else DEFAULT_ACCESSIONS
    if not target.exists():
        return {}
    return yaml.safe_load(target.read_text(encoding="utf-8")) or {}


def missing_provenance(data: Dict) -> List[str]:
    """Which fields are absent, blank, or still a placeholder."""
    provenance = data.get("provenance") or {}
    if not isinstance(provenance, dict):
        return list(REQUIRED_FIELDS)
    missing: List[str] = []
    for field in REQUIRED_FIELDS:
        # `gtdb_release` may also live at the top level (older skeleton shape).
        value = provenance.get(field)
        if value is None and field == "gtdb_release":
            value = data.get("gtdb_release")
        if not _is_filled(value):
            missing.append(field)
    return missing


def require_provenance(data: Dict, *, context: str = "", stream=None) -> List[str]:
    """Return the missing fields, printing a NOT-EXECUTED verdict when any.

    Raises ``SystemExit(2)`` on incomplete provenance: a benchmark whose data
    cannot be re-downloaded identically must not produce a precision/recall
    number at all.
    """
    missing = missing_provenance(data)
    if not missing:
        return missing
    out = stream or sys.stderr
    where = f" in {context}" if context else ""
    print(
        f"NOT EXECUTED{where}: provenance is incomplete — "
        f"missing or placeholder: {', '.join(missing)}.\n"
        f"        Record source_registry, doi, license, retrieved_utc and "
        f"checksum_sha256 in {DEFAULT_ACCESSIONS.name} under `provenance:` "
        f"before trusting any metric from this benchmark.",
        file=out,
    )
    raise SystemExit(2)


def summarize(data: Dict, stream=None) -> None:
    """Print the provenance block so it lands in the run record."""
    out = stream or sys.stdout
    provenance = data.get("provenance") or {}
    for field in REQUIRED_FIELDS:
        value = provenance.get(field)
        if value is None and field == "gtdb_release":
            value = data.get("gtdb_release")
        print(f"  provenance.{field} = {value if value else 'MISSING'}", file=out)


def iter_placeholder_paths(data: Dict, downloads: Path) -> Iterable[Path]:
    """Yield downloaded files a checksum should be recorded for."""
    downloads = Path(downloads)
    if not downloads.exists():
        return []
    return sorted(
        path for path in downloads.rglob("*")
        if path.is_file() and path.suffix.lower() in (".faa", ".fasta", ".fna", ".gz")
    )


if __name__ == "__main__":
    # Standalone use: check the shipped file and report what is still missing.
    loaded = load_accessions()
    absent = missing_provenance(loaded)
    if absent:
        print(f"INCOMPLETE: {', '.join(absent)}")
        raise SystemExit(2)
    summarize(loaded)
    print("COMPLETE")
