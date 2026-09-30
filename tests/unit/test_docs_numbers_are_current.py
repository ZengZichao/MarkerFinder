"""The user-facing documents must not keep quoting retired measurements.

The docs contract requires the README/MANUAL pair to stay in sync with the code.
That was honoured by hand, and the hand-copy had rotted in seven places at once: three
``tests-401 passed`` badges, "41 test modules / 401 unit tests", "16 modules,
174 unit tests", and "integration/benchmark tests are placeholders" — the last
one contradicting the same file's own directory listing two paragraphs later.

This gate encodes the rule that keeps that from recurring, which is a
single-source rule rather than a copy-editing one (same shape as for
thresholds): **suite counts are stated in the READMEs only**, because they move
every time a test file lands; the manuals describe what the suite is and point
at the README instead of repeating a number that is stale by the next commit.
The manuals may still state the ``--check`` item count, which is stable.

Every check is two-sided: the retired-token scan is proved able to fire on a
planted string, and the "manuals carry no counts" rule is proved able to fire
on a planted sentence, before either is allowed to report clean.

Naming note: the ``test_control_*`` prefix marks planted positives; the
un-prefixed tests are the real checks.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
READMES = ["README.CN.md", "README.EN.md", "README.md"]
MANUALS = ["MANUAL.CN.md", "MANUAL.EN.md"]

# Wording that described the repository before the self-check layer landed.
# Listed in both languages on purpose: a one-language token list passes on the
# File it was written against and hides the other half of the documentation.
RETIRED = [
    "401 unit tests", "401 个单元测试", "41 test modules", "41 个测试模块",
    "174 unit tests", "174 tests", "174 个单元测试",
    "16 modules, 174", "16 个模块", "tests-401",
    "are placeholders", "占位状态", "(placeholder)", "（占位）",
]

# Any "<n> test cases / 用例 / collected cases / test modules / 测试模块" claim.
SUITE_CLAIM = re.compile(
    r"\d+\s*(?:个测试模块|test modules|个用例|collected cases|个单元测试|unit tests)"
)


def scan_retired(text: str) -> list[str]:
    return [token for token in RETIRED if token in text]


def _run(args: list[str], timeout: int = 600) -> subprocess.CompletedProcess:
    return subprocess.run(args, cwd=REPO, capture_output=True, text=True,
                          timeout=timeout)


def measured_counts() -> dict[str, int]:
    """(modules, unit, integration, benchmark, total) straight from the tree.

    ``-o addopts=""`` is required, not cosmetic: the project's pytest config
    sets ``addopts = "-v --tb=short"``, which overrides ``-q`` into a tree
    listing with no ``::`` node ids — the first version of this helper counted
    zero tests and only survived because an inner assertion caught it.
    """
    modules = len(list((REPO / "tests").rglob("test_*.py")))
    proc = _run([sys.executable, "-m", "pytest", "tests", "--collect-only",
                 "-q", "-p", "no:randomly", "-o", "addopts="])
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    counts = {"unit": 0, "integration": 0, "benchmark": 0}
    for line in proc.stdout.splitlines():
        for bucket in counts:
            if line.startswith(f"tests/{bucket}/") and "::" in line:
                counts[bucket] += 1
    total = sum(counts.values())
    assert total > 0, proc.stdout[:2000]
    # Pytest -q prints "N tests collected in X.XXs"; the non-quiet format is
    # "collected N items". Accept either so the gate survives an addopts change.
    declared = (re.search(r"(\d+) tests collected", proc.stdout)
                or re.search(r"collected (\d+) items?", proc.stdout))
    assert declared, proc.stdout[:2000]
    assert int(declared.group(1)) == total, (declared.group(1), counts)
    return {"modules": modules, "total": total, **counts}


def check_item_count() -> int:
    """How many scored items ``--check`` emits (INFO lines are commentary)."""
    proc = _run([sys.executable, "-m", "markerfinder", "--check"])
    # The exit code is deliberately ignored: on an out-of-range interpreter the
    # Version and ete3 items fail *by design*, and the item count still holds.
    return sum(proc.stdout.count(f"[{state}]") for state in ("PASS", "FAIL"))


REFERENCE_DB = REPO / "db" / "gtdb_markers"


def legitimate_item_counts() -> set[int]:
    """The two scored-item counts ``--check`` may print, depending on the install.

    The database-hash comparison is scored only where the reference databases are
    on disk. On a source-only checkout the profiles were never fetched — they are
    third-party models, assembled locally by ``scripts/fetch_marker_db.py`` — and
    the item reports INFO, because a comparison that could not run is not a pass.
    So the documented count is a pair, and a document that states only one of the
    two is wrong for half the installs.
    """
    live = check_item_count()
    with_db, without_db = (live, live - 1) if REFERENCE_DB.is_dir() else (live + 1, live)
    return {with_db, without_db}


def count_guard_reason() -> str | None:
    """Why the live-count comparison cannot run here, or None when it can.

    The badge quotes the FULL-suite collection count, which only exists where
    every module collects: inside the declared range (>=3.10,<3.13) with ete3
    importable. Anywhere else the ete3-dependent modules skip and the live
    total legitimately differs from the badge -- failing there would let an
    interpreter that cannot even measure the number veto a doc sync. Per the
    suite-wide convention that is NOT EXECUTED, and never a silent pass.
    """
    if sys.version_info >= (3, 13):
        return (f"NOT EXECUTED: interpreter {sys.version_info[:3]} is outside "
                "the declared range (>=3.10,<3.13), so ete3 cannot import and "
                "the suite cannot collect fully here")
    try:
        import ete3  # Noqa: F401
    except Exception:
        return ("NOT EXECUTED: ete3 is not importable here, so the "
                "ete3-dependent modules skip and the live total differs "
                "from the badge")
    return None


def test_control_retired_scan_fires_on_planted_text():
    """A negative scan without a planted positive is indistinguishable from dead."""
    planted = "41 个测试模块共 401 个单元测试；集成测试与基准测试为占位状态"
    hits = scan_retired(planted)
    assert hits == ["401 个单元测试", "41 个测试模块", "占位状态"], hits
    assert "401 unit tests" not in hits, "English tokens must not match Chinese prose"
    assert scan_retired("81 个测试模块，863 个用例") == []


def test_control_the_suite_claim_pattern_fires():
    planted = ["863 collected cases across 81 test modules", "81 个测试模块，863 个用例"]
    for text in planted:
        assert SUITE_CLAIM.search(text), text
    assert SUITE_CLAIM.search("运行 48 项检查") is None, "--check count is not a suite count"


def test_no_document_quotes_a_retired_measurement():
    for name in READMES + MANUALS:
        text = (REPO / name).read_text(encoding="utf-8")
        assert scan_retired(text) == [], f"{name} still advertises {scan_retired(text)}"


def test_manuals_do_not_restate_suite_counts():
    """The single-source half of: numbers live in the README, prose here."""
    for name in MANUALS:
        text = (REPO / name).read_text(encoding="utf-8")
        offenders = [m.group(0) for m in SUITE_CLAIM.finditer(text)]
        assert offenders == [], f"{name} restates suite counts: {offenders}"


def test_control_the_count_guard_skips_outside_its_measurable_envelope():
    """Positive control for the NOT EXECUTED branch: on an out-of-range
    interpreter the comparison must skip with a named reason, not fail and not
    pretend the numbers agree."""
    import unittest.mock as mock

    with mock.patch.object(sys, "version_info", (3, 14, 6, "final", 0)):
        out_of_range = count_guard_reason()
    assert out_of_range is not None and "NOT EXECUTED" in out_of_range


def test_readme_badges_quote_the_measured_total():
    skip_reason = count_guard_reason()
    if skip_reason:
        pytest.skip(skip_reason)
    counts = measured_counts()
    seen = 0
    for name in READMES:
        text = (REPO / name).read_text(encoding="utf-8")
        # The badge states the COLLECTED count: what a reader can reproduce with
        # ``pytest --collect-only``. Passed/skipped splits depend on the machine
        # (optional external tools and artefacts), so a badge would go stale.
        for badge in re.findall(r"tests-(\d+)%20collected", text):
            seen += 1
            assert int(badge) == counts["total"], f"{name} badge {badge} != {counts}"
    assert seen == len(READMES), f"expected one badge per README, saw {seen}"
    for name in READMES:
        text = (REPO / name).read_text(encoding="utf-8")
        assert "%20passed" not in text, (
            f"{name} still carries a 'passed' badge; state the collected count "
            "and the pass/skip split in prose, where it can be qualified"
        )


def test_readme_suite_numbers_match_the_live_collection():
    skip_reason = count_guard_reason()
    if skip_reason:
        pytest.skip(skip_reason)
    counts = measured_counts()
    for name in READMES[:2]:  # The two localised READMEs carry the prose
        text = (REPO / name).read_text(encoding="utf-8")
        assert f"{counts['total']}" in text, f"{name} never mentions {counts['total']}"
        assert f"{counts['modules']}" in text, f"{name} never mentions its module count"


def test_documented_check_item_count_matches_what_check_emits():
    live = check_item_count()
    legit = legitimate_item_counts()
    assert live in legit, f"--check emitted {live} scored items, outside {sorted(legit)}"
    for name in READMES[:2] + MANUALS:
        text = (REPO / name).read_text(encoding="utf-8")
        claimed = {int(n) for n in
                   re.findall(r"(\d+) (?:项检查|项；|items;|items\b)", text)}
        stale = claimed - legit
        assert not stale, (
            f"{name} claims {sorted(stale)} items, which is neither of the "
            f"{sorted(legit)} counts --check can emit")
        assert live in claimed, (
            f"{name} never states the {live} items this install reports "
            f"(it claims {sorted(claimed)})")
        assert claimed, f"{name} should state the {live}-item self-check"


# --- the badge and its own breakdown must agree ----------------------------
#
# ``test_readme_badges_quote_the_measured_total`` compares the badge with a LIVE
# Collection count, and ``..._suite_numbers_match_the_live_collection`` only
# Checks that the total appears somewhere in the prose. Neither of them reads the
# Per-directory breakdown, so a README could say "897 cases -- 831 unit, 14
# Integration, 51 benchmark" (which is what it said: 831+14+51 = 896) and every
# Existing gate would report clean in an environment where the numbers cannot be
# Measured, and still report clean in one where they can.
#
# These two checks close that hole: the sentence must be internally consistent
# (no environment needed), and its buckets must match the live per-directory
# Counts (same measurable envelope as the badge check).

BUDGET_BUCKETS = ("unit", "integration", "benchmark")
BROKEN = re.compile(r"tests-(\d+)%20collected")
BUCKET_TOKEN = re.compile(r"tests/(unit|integration|benchmark)/?")


def breakdown_claims(text: str) -> dict[str, int]:
    """Per-directory counts as advertised, reading either word order.

    Handles ``其中 `tests/unit/` 831 个`` (number follows the path) and
    ``831 in `tests/unit/``` (number precedes it). The following position is
    tried first: a number *before* a path token is the previous bucket's count
    half the time, and reading it as this bucket's would make the check lie.
    """
    claims: dict[str, int] = {}
    for match in BUCKET_TOKEN.finditer(text):
        tail = text[match.end():]
        after = re.match(r"`?\s*(\d{2,4})", tail)
        if after:
            value = int(after.group(1))
        else:
            head = text[max(0, match.start() - 16):match.start()]
            before = re.search(r"(\d{2,4})\s*(?:in\s*)?[`'\"]?$", head)
            if not before:
                continue
            value = int(before.group(1))
        claims.setdefault(match.group(1), value)
    return claims


def test_control_a_badge_that_contradicts_its_own_breakdown_is_caught():
    """The exact shipped string, verbatim, must be rejected -- not reinterpreted."""
    historical = (
        "[![CI](https://img.shields.io/badge/tests-897%20collected-brightgreen.svg)]\n"
        "**897** cases — 831 in `tests/unit/`, 14 in `tests/integration/` and "
        "51 in `tests/benchmark/`."
    )
    badge = int(BROKEN.search(historical).group(1))
    claims = breakdown_claims(historical)
    assert claims == {"unit": 831, "integration": 14, "benchmark": 51}, claims
    assert sum(claims.values()) == 896 != badge, (
        "the control is toothless: 897 vs 831+14+51=896 must be detectable"
    )


def test_control_breakdown_reading_handles_the_chinese_word_order():
    text = "其中 `tests/unit/` 831 个、`tests/integration/` 14 个、`tests/benchmark/` 51 个。"
    assert breakdown_claims(text) == {
        "unit": 831, "integration": 14, "benchmark": 51,
    }


def test_every_readme_badge_equals_its_own_breakdown():
    for name in READMES:
        text = (REPO / name).read_text(encoding="utf-8")
        badges = {int(b) for b in BROKEN.findall(text)}
        claims = breakdown_claims(text)
        if not claims:
            continue  # A README with no per-directory sentence is not self-contradicting
        assert len(badges) == 1, f"{name} carries more than one badge total: {badges}"
        badge = badges.pop()
        assert set(claims) == set(BUDGET_BUCKETS), f"{name}: missing bucket(s) {claims}"
        assert sum(claims.values()) == badge, (
            f"{name}: badge says {badge} but its own breakdown sums to "
            f"{sum(claims.values())} ({claims})"
        )


def test_readme_breakdown_matches_the_live_per_directory_counts():
    skip_reason = count_guard_reason()
    if skip_reason:
        pytest.skip(skip_reason)
    counts = measured_counts()
    for name in READMES[:2]:  # The localised READMEs carry the prose
        claims = breakdown_claims((REPO / name).read_text(encoding="utf-8"))
        assert claims, f"{name} states no per-directory breakdown at all"
        assert claims == {b: counts[b] for b in BUDGET_BUCKETS}, (
            f"{name} breakdown {claims} != live {dict((b, counts[b]) for b in BUDGET_BUCKETS)}"
        )
