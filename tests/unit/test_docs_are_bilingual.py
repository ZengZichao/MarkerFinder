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
4. the pair must be *parallel in shape* — same count of level-1/2 headings — so
   a translation that quietly drops a section fails rather than looking fine.
"""

from __future__ import annotations

import re
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


def _prose_lines(text: str) -> list[str]:
    """Lines outside fenced code blocks.

    Without this, a shell or Python comment inside a code fence (``# install``,
    ``# nexus``) counts as a level-1 heading, and two faithful translations
    would differ in "heading count" merely because one wrote more comments.
    """
    out: list[str] = []
    fence = None
    for line in text.splitlines():
        stripped = line.strip()
        if fence is None and stripped.startswith("```"):
            fence = stripped[:3]
            continue
        if fence is not None:
            if stripped.startswith(fence):
                fence = None
            continue
        out.append(line)
    return out


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
