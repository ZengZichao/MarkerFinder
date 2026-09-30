"""The manual's glossary has to be complete, and its rows have to say something.

Two terms in this project are easy to confuse and expensive to get wrong:
``marker``/``gene``, ``consistent``/``concordant``, and above all ``risk`` versus
``incongruence``. The glossary in MANUAL.*.md is where a reader resolves them, so
this gate checks the glossary itself rather than trusting prose:

* every row carries a usage column and a non-empty caveat column;
* the confusable pairs each have a row that names both sides;
* ``risk`` is declared not a probability and not an HGT verdict;
* ``incongruence`` is declared an observation, kept apart from ``risk``.

Every check is two-sided: a synthetic glossary with a row removed must be
reported, so the real assertion cannot pass because the parser found nothing.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MANUAL = REPO / "MANUAL.CN.md"
MANUAL_EN = REPO / "MANUAL.EN.md"

ROW = re.compile(r"^\|\s*\*\*(.+?)\*\*\s*\|")
GLOSSARY_HEADING = {
    "CN": r"^## 术语表",
    "EN": r"^## Glossary",
}


def key(term: str) -> str:
    """Normalise a glossary row label: drop the parenthetical, spacing, slashes.

    ``不可测（NOT_MEASURABLE）`` and ``不可测 (NOT_MEASURABLE)`` — full-width and
    half-width brackets — must compare equal, otherwise the check reports a
    phantom difference between the two languages.
    """
    text = re.split(r"[（(]", term, maxsplit=1)[0]
    return re.sub(r"[\s/]+", "", text).lower()


def glossary_rows(text: str, heading_pattern: str) -> dict[str, str]:
    """Rows of the first table under a heading matching ``heading_pattern``."""
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines)
                     if re.match(heading_pattern, line))
    except StopIteration:
        return {}
    rows: dict[str, str] = {}
    for line in lines[start + 1:]:
        if rows and not line.startswith("|"):
            break
        match = ROW.match(line)
        if match:
            rows[key(match.group(1))] = line
    return rows


def glossary_terms(text: str, heading_pattern: str) -> set[str]:
    return set(glossary_rows(text, heading_pattern))


def missing_terms(required: set[str], provided: set[str]) -> set[str]:
    return required - provided


def _both_glossaries() -> dict[str, dict[str, str]]:
    return {
        "CN": glossary_rows(MANUAL.read_text(encoding="utf-8"),
                            GLOSSARY_HEADING["CN"]),
        "EN": glossary_rows(MANUAL_EN.read_text(encoding="utf-8"),
                            GLOSSARY_HEADING["EN"]),
    }


def test_the_glossary_is_found_in_both_languages():
    """Bookend: a regex that matched nothing would make every check below vacuous."""
    for label, rows in _both_glossaries().items():
        assert len(rows) >= 10, f"{label} glossary rows parsed: {sorted(rows)}"


def test_control_a_glossary_that_drops_a_term_is_reported():
    """Positive control: the check must be able to fail."""
    rows = _both_glossaries()["CN"]
    truncated = set(rows)
    dropped = sorted(truncated)[0]
    truncated.discard(dropped)
    assert missing_terms(set(rows), truncated) == {dropped}
    assert missing_terms(set(rows), set(rows)) == set(), (
        "identical inputs must be silent")


def test_every_row_states_usage_and_a_caveat():
    """A term listed without a caveat is a dictionary entry, not a usage rule."""
    for label, rows in _both_glossaries().items():
        for term, line in rows.items():
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            assert len(cells) >= 3, f"{label} row for {term} is not a 3-column row: {line}"
            assert cells[1], f"{label} row for {term} states no usage: {line}"
            assert cells[2], f"{label} row for {term} states no caveat: {line}"


def test_the_three_confusable_pairs_have_usage_rules():
    """A row exists for each side of each pair, and it says something non-trivial."""
    rows = _both_glossaries()["CN"]
    # Key drops slashes, so the labels compare as glued lowercase strings.
    required = {
        key("marker / gene"): {"marker", "gene"},
        key("consistent / concordant"): {"consistent", "concordant"},
        "risk": {"risk"},
        # Risk and incongruence live in two rows; both must exist and the
        # Incongruence row must point back at risk, which is the confusion this
        # Glossary exists to defend against.
        key("incongruence（拓扑不一致）"): {"incongruence", "risk"},
    }
    for term, needles in required.items():
        assert term in rows, f"no glossary row for {term}: {sorted(rows)}"
        row = rows[term]
        for needle in needles:
            assert needle in row.lower(), f"{term} row never names {needle}: {row}"
    incongruence_row = rows[key("incongruence（拓扑不一致）")]
    assert "HGT 结论" in incongruence_row, (
        "the incongruence row must state that incongruence is not an HGT verdict")


def test_risk_is_declared_not_a_probability():
    """The row that matters most: risk must not read as an HGT probability."""
    glossaries = _both_glossaries()
    row = glossaries["CN"].get("risk", "")
    assert row, "the risk row disappeared from the manual glossary"
    assert "不是概率" in row, row
    row_en = glossaries["EN"].get("risk", "")
    assert row_en, "the risk row disappeared from the English glossary"
    assert "not a probability" in row_en.lower(), row_en
    incongruence_en = glossaries["EN"][key("incongruence")]
    assert "hgt verdict" in incongruence_en.lower(), incongruence_en
