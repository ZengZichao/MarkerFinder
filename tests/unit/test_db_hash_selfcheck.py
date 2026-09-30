"""A wrong database hash must make --check red.

The hash COMPARISON exists in ``verify_database``, and before this test the run
path warned, but ``--check`` never looked at it, and
``db/expected_hashes.json`` shipped with empty values — so 's second
bullet ("数据库哈希不符时 --check 红") could not happen. Each case below is
two-sided: the same helper must report PASS for a matching pin, FAIL for a
stale one, and must NOT manufacture a pass when nothing was pinned.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from markerfinder.cli.self_test import _test_database_hashes


def _combined_hash(directory: Path) -> str:
    """Reproduce DatabaseVersionManager's directory hashing rule."""
    per_file = [
        hashlib.sha256(f.read_bytes()).hexdigest()
        for f in sorted(directory.rglob("*")) if f.is_file()
    ]
    return hashlib.sha256("".join(per_file).encode()).hexdigest()


def _make_db(tmp_path, payload="alpha"):
    db = tmp_path / "db"
    marker_dir = db / "gtdb_markers"
    marker_dir.mkdir(parents=True)
    (marker_dir / "PF00410.20.hmm").write_text(payload, encoding="utf-8")
    return db


def _write_expected(db, values):
    (db / "expected_hashes.json").write_text(
        json.dumps(values), encoding="utf-8", newline="\n",
    )


def _status(results, prefix):
    return [
        status for name, status, _detail in results
        if name.startswith(prefix)
    ]


def test_matching_hash_passes(tmp_path):
    db = _make_db(tmp_path)
    correct = _combined_hash(db / "gtdb_markers")
    _write_expected(db, {"gtdb_markers": correct})

    statuses = _status(_test_database_hashes(str(db)), "Database hash [gtdb_markers]")
    assert statuses == ["PASS"], statuses


def test_stale_hash_fails_check(tmp_path):
    """The acceptance case: pinned hash no longer matches the bytes."""
    db = _make_db(tmp_path, payload="alpha")
    _write_expected(db, {"gtdb_markers": _combined_hash(db / "gtdb_markers")})
    # Database content changes after the pin (a deliberate or accidental update).
    (db / "gtdb_markers" / "PF00410.20.hmm").write_text("mutated", encoding="utf-8")

    results = _test_database_hashes(str(db))
    statuses = _status(results, "Database hash [gtdb_markers]")
    assert statuses == ["FAIL"], statuses
    detail = next(
        d for n, s, d in results if n.startswith("Database hash [gtdb_markers]")
    )
    assert "!=" in detail and "database changed" in detail


def test_unpopulated_pin_is_reported_as_not_run(tmp_path):
    """Shipped skeleton state: never a fabricated green."""
    db = _make_db(tmp_path)
    _write_expected(db, {"gtdb_markers": "", "marker_hmm": ""})

    results = _test_database_hashes(str(db))
    assert _status(results, "Database SHA-256 comparison") == ["INFO"]
    detail = results[0][2]
    assert "NOT RUN" in detail and "not a pass" in detail
    # No PASS may appear while nothing is pinned.
    assert "PASS" not in [status for _n, status, _d in results]


def test_declared_pin_for_absent_database_is_not_silently_passed(tmp_path):
    db = _make_db(tmp_path)
    _write_expected(db, {"gtdb_markers": "0" * 64})
    (db / "gtdb_markers" / "PF00410.20.hmm").unlink()
    shutil_rmdir = db / "gtdb_markers"
    shutil_rmdir.rmdir()

    results = _test_database_hashes(str(db))
    statuses = _status(results, "Database hash [gtdb_markers]")
    assert statuses == ["INFO"], statuses
    assert "NOT CHECKED" in next(
        d for n, s, d in results if n.startswith("Database hash [gtdb_markers]")
    )


def test_the_check_is_registered_in_the_self_test_suite():
    """A check nobody calls is not a check (baseline lesson).

    Tightened when ``--check --db-dir X`` was found to ignore X and verify the
    repository's default ``db/`` anyway: the lock now also requires the call site
    to hand the user's directory to the check, so a silent revert to the
    wrong-artifact verdict fails here rather than in someone's review.
    """
    base = Path(__file__).resolve().parents[2]
    self_test = (base / "markerfinder/cli/self_test.py").read_text(encoding="utf-8")
    assert re.search(
        r"results \+= _test_database_hashes\(\s*db_dir\s*\)", self_test
    ), "the database-hash check is defined but not run by --check"
    main = (base / "markerfinder/cli/main.py").read_text(encoding="utf-8")
    assert re.search(r"_run_self_test\(\s*\w+\s*\)", main), (
        "--check must forward an explicit --db-dir; calling _run_self_test() with "
        "no argument silently checks the default database instead"
    )
