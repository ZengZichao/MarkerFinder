"""Every document in this repository must exist in both languages, English first.

The rule this enforces is a published one: documentation ships in an English and
a Chinese version, and the English version is the primary one. That is easy to
state and easy to break — a new ``docs/`` page written only in Chinese, a
``README.EN.md`` whose Chinese twin was never updated, a link that lands on the
Chinese page from the English index. None of those are visible to a reviewer
reading one file.

So the check is structural, over the filesystem:

1. every ``*.md`` under the repository root and ``docs/`` that carries a
   ``.EN.md`` / ``.CN.md`` suffix must have both halves;
2. the language-neutral entry point (``README.md``) must be English prose and
   must link to both language versions;
3. no document may be English-only or Chinese-only;
4. the pair must be *parallel in shape* -- same count of level-1/2 headings -- so
   a translation that quietly drops a section fails rather than looking fine;
5. the parts of a document that have no translation at all -- command-line
   flags, identifiers, measured numbers -- must be the same in both halves.

Check 4 is what a reviewer can eyeball, and it is not enough: a dropped
"only when", a threshold mistyped as 0.56, and a flag transliterated into
Chinese characters each leave every heading, line, and fence count untouched.
Check 5 is therefore mechanical.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

# Directories that are documentation, not code comments or generated output.
DOC_DIRS = [REPO, REPO / "docs", REPO / "validation", REPO / "tests" / "benchmark",
            REPO / "scripts"]
SKIP_PARTS = {".git", "__pycache__", ".work", "node_modules", "outputs",
              "results", ".pytest_cache", ".mimosa"}

CJK = re.compile(r"[\u4e00-\u9fff]")
PAIR = re.compile(r"^(?P<stem>.+)\.(?P<lang>EN|CN)\.md$")

# ── what the semantic check compares ─────────────────────────────────────
#
# Numbers that carry a decimal point, a percentage, or a data-size unit. Time
# units are deliberately absent: "23 min" and "23 分钟" are the same fact
# written in two languages, and teaching this regex to translate unit words
# would make it the very thing it exists to catch.
FIGURE = re.compile(r"\d+\.\d+|\d+(?:\.\d+)?\s*(?:%|GB|MB|TB|kb|bp)")
RANGE = re.compile(r"(?<![\d.])(\d+(?:\.\d+)?)\s*(%|GB|MB|TB|kb|bp)?"
                   r"\s*(?:[-–—~～]|\.\.\.)\s*"
                   r"(\d+(?:\.\d+)?)\s*(%|GB|MB|TB|kb|bp)")
# A flag is the one thing in a command line that cannot be paraphrased.
FLAG = re.compile(r"(?<![\w/])-{1,2}[A-Za-z][\w-]*")
CODE_SPAN = re.compile(r"`([^`\n]+)`")
ASCII_ONLY = re.compile(r"^[\x20-\x7e]+$")


def _markdown_files() -> list[Path]:
    found: list[Path] = []
    for base in DOC_DIRS:
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*.md")):
            if any(part in SKIP_PARTS for part in path.parts):
                continue
            found.append(path)
    return found


def _split_fenced(text: str) -> tuple[str, str]:
    """Return (prose, fenced bodies).

    Without this, a shell or Python comment inside a code fence (``# install``,
    ``# nexus``) counts as a level-1 heading, and two faithful translations
    would differ in "heading count" merely because one wrote more comments.
    """
    prose: list[str] = []
    blocks: list[str] = []
    fence = None
    buf: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if fence is None and stripped.startswith("```"):
            fence = stripped[:3]
            buf = []
            continue
        if fence is not None:
            if stripped.startswith(fence):
                blocks.append("\n".join(buf))
                fence = None
            else:
                buf.append(line)
            continue
        prose.append(line)
    return "\n".join(prose), "\n".join(blocks)


def _prose_lines(text: str) -> list[str]:
    return _split_fenced(text)[0].splitlines()


def _expand_ranges(text: str) -> str:
    """Write every numeric range as its two endpoints, unit attached to each.

    English writes ``1-12 GB`` and Chinese writes ``1 GB～12 GB``. Both name
    the same interval, so the normalising belongs here rather than in the
    documents; without it the check would flag a style difference it cannot
    tell from a wrong number.
    """
    return RANGE.sub(
        lambda m: f"{m.group(1)}{m.group(2) or m.group(4)} {m.group(3)}{m.group(4)}",
        text,
    )


def _figures(text: str) -> set:
    """The measured numbers a translation is not free to change.

    Spaces inside a token are dropped, so "12 GB", "12GB" and "12  GB" are one
    value: whether the unit is written with a gap is a typographic habit of
    the language, not a fact about the pipeline.
    """
    prose, _fenced = _split_fenced(text)
    return {found.replace(" ", "") for found in FIGURE.findall(_expand_ranges(prose))}


def _flags(text: str) -> Counter:
    """Every command-line flag inside a fenced block.

    Only fences are scanned for bare flags: a flag mentioned in prose arrives
    as an inline code span, and _identifiers compares those.
    """
    _prose, fenced = _split_fenced(text)
    return Counter(FLAG.findall(fenced))


def _identifiers(text: str) -> Counter:
    """ASCII-only inline code spans in prose: names, paths, flags, columns.

    A span with no non-Latin character in it is not prose, so it must survive
    translation unchanged. Spans that do carry Chinese (a quoted program
    output line, which this repository translates on purpose) are skipped.
    """
    prose, _fenced = _split_fenced(text)
    return Counter(span for span in CODE_SPAN.findall(prose) if ASCII_ONLY.match(span))


def _headings(text: str, level: int) -> list[str]:
    prefix = "#" * level + " "
    return [ln for ln in _prose_lines(text) if ln.startswith(prefix)]


def _is_mostly_english(text: str) -> bool:
    """English-primary means the prose backbone is Latin script.

    A Chinese document is allowed to quote English identifiers; an English one
    is not allowed to carry Chinese sentences, because the English file is the
    one users are pointed at first. Code fences are ignored: a translation may
    keep a shell comment in the original language.
    """
    body = [ln for ln in _prose_lines(text) if ln.strip()]
    if not body:
        return True
    cjk_lines = sum(1 for ln in body if len(CJK.findall(ln)) >= 5)
    return cjk_lines / len(body) < 0.05


def _pairs() -> list[tuple[str, Path, Path]]:
    out = []
    for path in _markdown_files():
        match = PAIR.match(path.name)
        if not match:
            continue
        stem = match.group("stem")
        other = "CN" if match.group("lang") == "EN" else "EN"
        out.append((stem, path, path.with_name(f"{stem}.{other}.md")))
    return out


@pytest.mark.parametrize("stem, this, sibling",
                         _pairs(),
                         ids=[f"{s[0]}" for s in _pairs()])
def test_every_bilingual_document_has_both_halves(stem, this, sibling):
    assert sibling.exists(), (
        f"{this.name} exists but its counterpart {sibling.name} does not. "
        "Every document ships in both languages."
    )


@pytest.mark.parametrize("stem, this, sibling",
                         _pairs(),
                         ids=[f"{s[0]}" for s in _pairs()])
def test_the_two_languages_are_parallel_in_shape(stem, this, sibling):
    a = this.read_text(encoding="utf-8")
    b = sibling.read_text(encoding="utf-8")
    for level in (1, 2):
        na, nb = len(_headings(a, level)), len(_headings(b, level))
        # A translation is allowed to fold or add at most one section; a large
        # Divergence means one language is carrying content the other dropped.
        assert abs(na - nb) <= 1, (
            f"{this.name} has {na} level-{level} headings, {sibling.name} has "
            f"{nb}: the two versions are no longer the same document"
        )
    assert abs(len(a.splitlines()) - len(b.splitlines())) < len(a.splitlines()), (
        f"{this.name} and {sibling.name} differ by more than a whole document's "
        "length; one of them is a stub"
    )


def test_english_halves_are_actually_english():
    offenders = []
    for stem, this, _sibling in _pairs():
        if not this.name.endswith(".EN.md"):
            continue
        if not _is_mostly_english(this.read_text(encoding="utf-8")):
            offenders.append(this.name)
    assert not offenders, (
        f"these English documents contain Chinese prose and are therefore not "
        f"the primary version: {offenders}"
    )


def test_no_bilingual_document_is_language_lone():
    """A ``foo.md`` with no language suffix must not be a third, orphan copy."""
    singles = [
        p.name for p in _markdown_files()
        if not PAIR.match(p.name) and p.name != "README.md"
        and not _is_mostly_english(p.read_text(encoding="utf-8"))
    ]
    assert not singles, (
        f"these unsuffixed documents are neither English nor Chinese and cannot "
        f"be paired: {singles}. Name them foo.EN.md / foo.CN.md, or keep them "
        "out of the documentation directories."
    )


def test_root_index_is_english_and_links_both_language_versions():
    index = REPO / "README.md"
    text = index.read_text(encoding="utf-8")
    assert _is_mostly_english(text), (
        "README.md is the entry point and must be the English-primary page"
    )
    assert re.search(r"README\.EN\.md", text), "README.md must link the EN page"
    assert re.search(r"README\.CN\.md", text), "README.md must link the CN page"
    assert re.search(r"MANUAL\.EN\.md", text), "README.md must link the manual"
    # The link order is the display order: English first.
    assert text.index("README.EN.md") < text.index("README.CN.md"), (
        "the English version must be linked before the Chinese one"
    )


# ── semantic parity: the parts of a document that have no translation ────

@pytest.mark.parametrize("stem, this, sibling",
                         _pairs(),
                         ids=[f"{s[0]}" for s in _pairs()])
def test_measured_numbers_are_the_same_in_both_languages(stem, this, sibling):
    """A threshold, a version, and a percentage are not style choices.

    Compared as sets, not counts: "3.10 through 3.14" and
    "3.10、3.11、3.12、3.13、3.14" are the same statement about the same
    interpreters, and demanding equal repetition would flag that as a defect.
    The trade-off is on purpose -- a figure whose value is wrong or absent is
    caught, and a duplicated figure that merely stops being repeated on one
    side is not.
    """
    a = _figures(this.read_text(encoding="utf-8"))
    b = _figures(sibling.read_text(encoding="utf-8"))
    assert a == b, (
        f"{this.name} states {sorted(a - b)} that {sibling.name} does not, and "
        f"{sibling.name} states {sorted(b - a)} that {this.name} does not: one "
        "of the two has a number the other does not"
    )


@pytest.mark.parametrize("stem, this, sibling",
                         _pairs(),
                         ids=[f"{s[0]}" for s in _pairs()])
def test_command_line_flags_are_not_translated(stem, this, sibling):
    a = _flags(this.read_text(encoding="utf-8"))
    b = _flags(sibling.read_text(encoding="utf-8"))
    assert a == b, (
        f"the commands in {this.name} take {sorted(a - b)} and the ones in "
        f"{sibling.name} take {sorted(b - a)}: a flag is not translatable, so "
        "one of the two documents is telling users something the tool does not "
        "accept"
    )


@pytest.mark.parametrize("stem, this, sibling",
                         _pairs(),
                         ids=[f"{s[0]}" for s in _pairs()])
def test_identifiers_in_prose_are_not_localized(stem, this, sibling):
    a = _identifiers(this.read_text(encoding="utf-8"))
    b = _identifiers(sibling.read_text(encoding="utf-8"))
    assert a == b, (
        f"these backticked names exist in {this.name} only: {sorted(a - b)}, "
        f"and in {sibling.name} only: {sorted(b - a)}. A column header, path, "
        "or output label written in one language's characters is a name the "
        "software never emits"
    )


# ── controls: each semantic check has to be able to fail ─────────────────
#
# The three checks above read real documents, so a regex that matches nothing
# would make all three green forever -- the failure mode this repository has
# already had to write two other control tests about.

EN_SAMPLE = """\
# Sample

| item | value |
|---|---|
| threshold | 0.55 |
| free | 30-70% |
| space | 1-12 GB |

Run it with:

```bash
markerfinder run --min-occupancy 0.55 --threads 4
```

The column is `monophyly_proportion` and the file is `Phase5_reports/mag_quality.tsv`.
"""


def _cn_sample(**changes):
    text = EN_SAMPLE.replace("Sample", "样例").replace("threshold", "阈值")
    for old, new in changes.items():
        text = text.replace(old, new)
    return text


def test_control_the_figure_check_fires_on_a_mistyped_threshold():
    en = _figures(EN_SAMPLE)
    assert "0.55" in en and "30%" in en and "1GB" in en, en
    assert _figures(_cn_sample(**{"0.55": "0.56"})) != en
    assert _figures(_cn_sample(**{"30-70%": "30-90%"})) != en
    assert _figures(_cn_sample(**{"1-12 GB": "12 GB"})) != en
    # ... and stays quiet when only the wording of the range differs.
    assert _figures(_cn_sample(**{"1-12 GB": "1 GB～12 GB"})) == en


def test_control_the_flag_check_fires_on_a_translated_or_invented_flag():
    en = _flags(EN_SAMPLE)
    assert en["--min-occupancy"] == 1 and en["--threads"] == 1, en
    assert _flags(_cn_sample(**{"--threads 4": "--线程 4"})) != en
    assert _flags(_cn_sample(**{"--threads": "--threads-all"})) != en
    assert _flags(_cn_sample(**{"--threads 4": "4"})) != en


def test_control_the_identifier_check_fires_on_a_localized_name():
    en = _identifiers(EN_SAMPLE)
    assert "monophyly_proportion" in en and "Phase5_reports/mag_quality.tsv" in en, en
    assert _identifiers(_cn_sample(**{
        "`monophyly_proportion`": "`单系性比例`"})) != en
    assert _identifiers(_cn_sample(**{
        "`Phase5_reports/mag_quality.tsv`": "`Phase5_报告/mag_quality.tsv`"})) != en
    # A Chinese span is prose, not a name, and is left out of the comparison.
    assert _identifiers(_cn_sample(**{
        "The column is": "该列名为"})) == en
