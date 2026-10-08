#!/usr/bin/env python3
"""Build a genuinely anonymous peer-review package.

Why
---
A manuscript that tells reviewers it supplied "an anonymised package" is making a
checkable claim. Zipping the repository does not deliver it: author name, e-mail,
ORCID iD, own-repository URLs and BibTeX fields sit in the READMEs, ``pyproject.toml``,
the package's version module and ``LICENSE``, and ``.git`` plus the local tool caches
add commit identities and machine paths on top. This script makes the claim true, and
then *proves* it: the last line it prints is ``identifying-string hits: N``, and a
non-zero N fails the build.

Design rule that shapes the whole file
--------------------------------------
**The identifiers are derived from the repository, never written into this script.**
The names, e-mail, institution token and ORCID iD are read from git metadata, from the
project's own version module, and from the remote URL. That matters twice over: the
tool works on any repository, and the package can ship this script without the script's
own source re-leaking the surname it is meant to remove. Structural patterns (e-mail
shape, ORCID-iD shape, absolute-path shape) are assembled with ``_ci()`` so their
expressions do not match themselves either.

Two hard limits
---------------
* Scrubbing happens **only in the exported copy**. Nothing under the repository root is
  ever modified or deleted here.
* Nothing is invented. What the tool cannot know is written down instead of being filled
  in with a plausible placeholder -- see ``REVIEW_PACKAGE_NOTES.md`` in the package.

Usage
-----
    python scripts/make_review_package.py --out /path/MarkerFinder-review.zip
    python scripts/make_review_package.py --ref v0.1.1 --out ...
    python scripts/make_review_package.py --worktree --out ...   # pre-commit content
    python scripts/make_review_package.py --out ... --hmm-dir db/gtdb_markers
    python scripts/make_review_package.py --out ... --audit-extra 'orcid'

``--ref`` (default ``HEAD``) exports through ``git archive``, so ``.git`` cannot be
included by construction. ``--worktree`` exports tracked + untracked-but-not-ignored
files read from disk; ``.git`` and ignored paths are still dropped explicitly.

``--hmm-dir`` closes the one reproducibility hole a git-derived package cannot otherwise
escape: ``db/gtdb_markers/`` is git-ignored (the TIGRFAM/Pfam profiles are third-party
work assembled by ``scripts/fetch_marker_db.py``), so ``--marker-mode hmm`` -- and with
it the worked example -- is not runnable from the package alone. When the author passes
that directory deliberately, it is copied in and hashed: per file, and in the same
aggregate form ``db/expected_hashes.json`` uses, so the library can be verified rather
than trusted.

Standard library only.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path
from typing import Iterable, Optional

REDACTED = "withheld for anonymous peer review"

#: Tool caches, VCS state and OS droppings that must never reach a reviewer.
DROP_DIR_NAMES = {".git", "__pycache__", ".mypy_cache", ".pytest_cache", ".work",
                  ".mimosa", "node_modules", ".venv", "venv", "build", "dist",
                  "htmlcov", ".tox", ".idea", ".vscode", "smoketest_output"}
DROP_SUFFIXES = {".pyc", ".pyo", ".DS_Store", ".swp", ".swo", ".bak", ".tmp"}
DROP_EXACT = {".DS_Store", "Thumbs.db", "desktop.ini", ".coverage"}
#: Trees dropped wholesale: the git-ignored third-party data, plus run output.
DROP_REL_PREFIXES = ("db/gtdb_markers/", "validation/data/genomes/",
                     "validation/data/hmms/", "tests/benchmark/downloads/",
                     "output/", ".mypy_cache/", ".pytest_cache/")

TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".yaml", ".yml", ".toml", ".tsv",
                 ".csv", ".cfg", ".ini", ".xml", ".html", ".sh", ".bib", ".tex"}


def log(msg: str) -> None:
    print(msg, flush=True)


def _ci(word: str) -> str:
    """A case-insensitive literal, written so the builder cannot match itself."""
    return "".join(f"[{c.lower()}{c.upper()}]" if c.isalpha() else re.escape(c)
                   for c in word)


def _esc(word: str) -> str:
    return re.escape(word)


# --- structural patterns (repo-independent, self-escaping) -------------------


def structural_patterns() -> list[tuple[str, re.Pattern]]:
    return [
        ("orcid iD", re.compile(r"\b\d{4}-\d{4}-\d{4}-\d{3}[\dX]\b")),
        ("orcid profile link",
         re.compile(r"https?://" + _ci("orcid") + r"\.org/\S+", re.I)),
        ("e-mail address",
         re.compile(r"[A-Za-z0-9_.+-]+@" + r"(?:[A-Za-z0-9-]+\.)+[A-Za-z]{2,}")),
        ("mailto link", re.compile(_ci("mailto") + r":\s*\S+", re.I)),
        ("local absolute path",
         re.compile(r"/(?:" + _ci("Users") + "|" + _ci("home") + "|"
                    + _ci("private/var/folders") + "|" + _ci("var/folders") + ")"
                    + r"/[^\s\"',;:)\]}]*")),
        ("windows absolute path",
         re.compile(r"[A-Za-z]:" + r"\\+" + _ci("Users") + r"\\+\w+")),
    ]


# --- identifier derivation ---------------------------------------------------

_ATTR = re.compile(r'^(__\w+)\s*=\s*["\']([^"\']*)["\']', re.M)

#: Technical strings that show up in e-mails and remotes but identify nobody. Without
#: this filter "users" or "github" would become an author-name token and the scrubber
#: would rewrite every ordinary occurrence of the word.
STOP_TOKENS = {"users", "noreply", "github", "gitlab", "com", "org", "net", "edu",
               "cn", "www", "http", "https", "localhost", "example", "outlook",
               "hotmail", "gmail", "qq", "foxmail", "mail", "e-mail", "no-reply"}


def _acceptable(value: str) -> bool:
    if len(value) < 3 or value.lower() in STOP_TOKENS or "@" in value:
        return False
    if value.isdigit():
        return False
    return True


def run_git(repo: Path, *args: str) -> str:
    proc = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if proc.returncode != 0:
        return ""
    return proc.stdout.strip()


def derived_tokens(repo: Path, version_module: Path) -> list[tuple[str, str]]:
    """Return (label, literal) identifier strings taken from the repository itself."""
    tokens: list[tuple[str, str]] = []

    def add(label: str, value: str) -> None:
        value = (value or "").strip()
        if not _acceptable(value) or any(value == v for _, v in tokens):
            return
        tokens.append((label, value))

    # Commit identity: name and e-mail of the author the tag was made under.
    for fmt in ("%an", "%ae", "%cn", "%ce"):
        out = run_git(repo, "log", "-1", "--format=" + fmt)
        add("commit identity", out)
        if "@" in out:
            local, _, domain = out.partition("@")
            add("commit e-mail user", local.split("+")[-1])
            add("institution domain", domain)
            add("institution token", domain.split(".")[0])

    # Remote URL: the owner segment is the author's GitHub namespace.
    remote = run_git(repo, "remote", "get-url", "origin")
    m = re.search(r"github\.[a-z]+/([^/]+)", remote, re.I)
    if m:
        add("repository owner", m.group(1))

    # The project's own version module: whatever it declares as authorship.
    if version_module.is_file():
        for name, value in _ATTR.findall(version_module.read_text(encoding="utf-8")):
            if not name.strip("_").startswith(("author", "email", "orcid", "maintainer")):
                continue
            if name.strip("_") == "orcid":
                add(f"{name} value", value)  # an iD is a whole token, never split
                continue
            add(f"{name} value", value)
            if "@" in value:
                local, _, domain = value.partition("@")
                add(f"{name} user", local)
                add(f"{name} domain", domain)
                add("institution token", domain.split(".")[0])
    return tokens


def identifier_rules(repo: Path, version_module: Path) -> list[tuple[str, re.Pattern]]:
    """Structural rules plus the repository-derived ones."""
    rules = structural_patterns()
    tokens = derived_tokens(repo, version_module)
    for label, literal in tokens:
        rules.append((label, re.compile(re.escape(literal), re.I)))
    # Split compound tokens (a name written with a separator inside it, an address
    # like "user@host") into their parts too, so a
    # README that writes only the surname is still caught. Parts must clear the same
    # acceptability bar as whole tokens, or "com" becomes an author-name rule.
    # Limitation worth knowing: a rule only exists for a spelling the repository
    # actually uses. If git metadata stores one concatenated word, a two-word
    # rendering of the same name is invisible here -- keep name variants out of this
    # file's own comments, or pass them with --audit-extra.
    for label, literal in tokens:
        for part in re.split(r"[\s_@.\-+]+", literal):
            if not _acceptable(part):
                continue
            if re.search(r"[一-鿿]", part):
                rules.append((f"{label} part", re.compile(part)))
            elif len(part) >= 5:
                rules.append((f"{label} part", re.compile(re.escape(part), re.I)))
    remote = run_git(repo, "remote", "get-url", "origin")
    m = re.search(r"github\.[a-z]+/([^/]+)", remote, re.I)
    if m:
        rules.append(("own-repository URL",
                      re.compile(r"github\.[a-z]+/" + re.escape(m.group(1)), re.I)))
    seen: set[str] = set()
    unique: list[tuple[str, re.Pattern]] = []
    for lab, pat in rules:
        if pat.pattern in seen:
            continue
        seen.add(pat.pattern)
        unique.append((lab, pat))
    return unique


# --- source collection --------------------------------------------------------


def files_from_ref(repo: Path, ref: str) -> list[tuple[str, bytes]]:
    tar_bytes = subprocess.run(["git", "-C", str(repo), "archive", "--format=tar", ref],
                               capture_output=True).stdout
    if not tar_bytes:
        raise SystemExit(f"git archive of {ref} produced nothing")
    out: list[tuple[str, bytes]] = []
    with tarfile.open(fileobj=io.BytesIO(tar_bytes)) as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            name = member.name
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise SystemExit(f"unsafe path in git archive: {name}")
            handle = tar.extractfile(member)
            if handle is not None:
                out.append((name.replace("\\", "/"), handle.read()))
    return out


def files_from_worktree(repo: Path) -> list[tuple[str, bytes]]:
    listing = subprocess.run(["git", "-C", str(repo), "ls-files", "--cached",
                              "--others", "--exclude-standard", "-z"],
                             capture_output=True).stdout
    out: list[tuple[str, bytes]] = []
    for raw in listing.split(b"\0"):
        if not raw:
            continue
        name = raw.decode()
        path = repo / name
        if path.is_file():
            out.append((name.replace("\\", "/"), path.read_bytes()))
    return out


def dropped(rel: str) -> bool:
    parts = rel.split("/")
    if any(part in DROP_DIR_NAMES for part in parts[:-1]):
        return True
    name = parts[-1]
    if name in DROP_EXACT or any(name.endswith(s) for s in DROP_SUFFIXES):
        return True
    if ".egg-info" in rel:
        return True
    return any(rel.startswith(p) for p in DROP_REL_PREFIXES)


# --- redaction ----------------------------------------------------------------


def is_text(rel: str, blob: bytes) -> bool:
    if Path(rel).suffix.lower() in TEXT_SUFFIXES:
        return True
    if blob[:4] in {b"%PDF", b"\x89PNG", b"GIF8", b"\xff\xd8"}:
        return False
    return b"\x00" not in blob[:2048] and len(blob) < 400_000


class Scrubber:
    def __init__(self, rules: list[tuple[str, re.Pattern]],
                 tokens: list[tuple[str, str]]) -> None:
        self.rules = rules
        self.tokens = tokens
        # The ORCID profile URL is masked before the bare word, and the author's own
        # repository URL before the owner name, so replacements never nest.
        self.own_repo = [pat for label, pat in rules if label == "own-repository URL"]

    def generic(self, text: str) -> str:
        # Own-repository URLs first, as a unit: masking just the owner would leave a
        # mangled URL that is worse for a reviewer and no less identifying in git terms.
        # The match starts at the host, so any "https://" already in the text stays put
        # -- a replacement carrying its own scheme would double it.
        for pat in self.own_repo:
            text = re.sub(pat.pattern + r"[^\s)\]\"'`]*",
                          "github.com/<" + REDACTED + ">", text)
        text = re.sub(r"https?://" + _ci("orcid") + r"\.org/\S+",
                      "[" + REDACTED + "]", text, flags=re.I)
        text = re.sub(r"\b\d{4}-\d{4}-\d{4}-\d{3}[\dX]\b",
                      "[ORCID iD " + REDACTED + "]", text)
        text = re.sub(_ci("mailto") + r":\s*[\w.+-]+@[\w-]+(?:\.[\w-]+)+",
                      "[e-mail " + REDACTED + "]", text, flags=re.I)
        text = re.sub(r"[A-Za-z0-9_.+-]+@(?:[A-Za-z0-9-]+\.)+[A-Za-z]{2,}",
                      "[e-mail " + REDACTED + "]", text)
        text = re.sub(r"/(?:" + _ci("Users") + "|" + _ci("home") + "|"
                      + _ci("private/var/folders") + "|" + _ci("var/folders") + ")"
                      + r"/[^\s\"',;:)\]}]*", "<local-path>", text)
        text = re.sub(r"[A-Za-z]:\\+" + _ci("Users") + r"\\+\w+\\*", "<local-path>", text)
        for label, literal in self.tokens:
            if "@" in literal:
                continue
            text = re.sub(re.escape(literal), "[author name " + REDACTED + "]",
                          text, flags=re.I)
        return text

    def __call__(self, rel: str, text: str) -> str:
        if rel.endswith("pyproject.toml"):
            text = re.sub(r"\nauthors\s*=\s*\[(?:[^\[\]]|\[[^\]]*\])*?\]", "", text,
                          flags=re.S)
            text = re.sub(r"\n\[project\.urls\]\n(?:[^\[\n]*\n)*", "\n\n", text)
        if Path(rel).name == "LICENSE":
            text = re.sub(r"^Copyright \(c\) (\d{4}).*$",
                          r"Copyright (c) \1 [author name " + REDACTED + "]",
                          text, flags=re.M)
        if rel.endswith("_version.py"):
            for attr in ("__author__", "__email__", "__orcid__"):
                text = re.sub(rf"^{attr}\s*=.*$",
                              f'{attr} = ""  # {REDACTED}', text, flags=re.M)
        if rel.endswith("__init__.py"):
            text = re.sub(r"^Author:.*$",
                          "Author: [" + REDACTED + "] (see LICENSE)", text, flags=re.M)
        if re.fullmatch(r"README(?:\.(?:EN|CN))?\.md", rel):
            text = self._readme(text)
        return self.generic(text)

    _AUTHOR_H = (r"#{1,3}\s*(?:" + _ci("Authors?") + "|" + _esc("作者") + "|"
                 + _esc("通讯作者") + r")")
    _CONTACT_H = (r"#{1,3}\s*(?:" + _ci("Contact") + "|" + _esc("联系方式") + "|"
                  + _ci("Support") + "|" + _esc("反馈") + r")")

    @staticmethod
    def _section(text: str, heading: str, note: str) -> str:
        return re.sub(rf"^({heading}[^\n]*\n)(?:.*?)(?=^\s*#{{1,3}}\s|\Z)",
                      lambda m: f"{m.group(1)}\n{note}\n", text, flags=re.M | re.S)

    def _readme(self, text: str) -> str:
        text = re.sub(r"^(\s*)" + _ci("author") + r"\s*=\s*\{[^}]*\},?\s*$",
                      r"\1author       = {[author name " + REDACTED + "]},",
                      text, flags=re.M | re.I)
        text = re.sub(r"^(\s*)" + _ci("orcid") + r"\s*=\s*\{[^}]*\},?\s*$",
                      r"\1note         = {identifiers " + REDACTED + "},",
                      text, flags=re.M | re.I)
        text = self._section(
            text, self._AUTHOR_H,
            "[Author identity " + REDACTED + ". The corresponding author's name, e-mail "
            "and ORCID iD are held by the editorial office and supplied with the "
            "accepted version.]")
        text = self._section(text, self._CONTACT_H,
                             "[Contact details " + REDACTED + ".]")
        text = re.sub(r"^- \*\*(?:" + _ci("Issues") + r"[^\n]*|"
                      + _ci("Email") + r"[^\n]*|" + _ci("ORCID") + r"[^\n]*|"
                      + _esc("问题与讨论") + r"[^\n]*|" + _esc("邮箱") + r"[^\n]*)\n",
                      "", text, flags=re.M | re.I)
        text = re.sub(r"git clone \S+", "git clone <repository URL " + REDACTED + ">",
                      text)
        return text


# --- the HMM library ----------------------------------------------------------


def aggregate_sha256(root: Path) -> tuple[str, list[tuple[str, str]]]:
    """SHA-256 over the concatenation of per-file digests in sorted path order.

    Identical to how ``markerfinder.utils.db_versioning`` computes the reference hash
    stored in ``db/expected_hashes.json``, so the number recorded in the package is
    directly comparable with what ``--check`` verifies.
    """
    per: list[tuple[str, str]] = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 16), b""):
                digest.update(chunk)
        per.append((str(path.relative_to(root).replace(os.sep, "/")), digest.hexdigest()))
    combined = hashlib.sha256()
    for _, hexdigest in per:
        combined.update(hexdigest.encode("ascii"))
    return combined.hexdigest(), per


def install_hmm_dir(repo: Path, src: Path, staged: dict[str, bytes]) -> dict:
    if not src.is_absolute():
        src = (repo / src).resolve()
    if not src.is_dir():
        raise SystemExit(f"--hmm-dir {src} is not a directory")
    aggregate, per = aggregate_sha256(src)
    for rel, _ in per:
        staged[f"db/gtdb_markers/{rel}"] = (src / rel).read_bytes()
    lines = ["# SHA-256 of each file under db/gtdb_markers/, sorted path order.",
             f"# Aggregate (sha256 over the concatenated per-file digests): {aggregate}",
             "# Compare with db/expected_hashes.json['gtdb_markers']. A mismatch is a",
             "# warning, not a failure -- see _note_gtdb_markers in that file.", ""]
    lines += [f"{d}  {r}" for r, d in per]
    staged["db/gtdb_markers.CHECKSUMS.txt"] = ("\n".join(lines) + "\n").encode()
    return {"provided": True, "source": str(src), "files": len(per),
            "aggregate_sha256": aggregate}


# --- audit --------------------------------------------------------------------


def audit(staged: dict[str, bytes], rules: list[tuple[str, re.Pattern]],
          extra: Optional[Iterable[re.Pattern]] = None,
          quiet: bool = False) -> tuple[dict[str, int], int, int]:
    """Count identifier occurrences across a set of file blobs.

    ``quiet`` suppresses the per-file lines, which is what the before-scrubbing
    positive control uses so a large count does not bury the run log.
    """
    extra = list(extra or ())
    per_rule = {label: 0 for label, _ in rules}
    total = extra_total = 0
    for rel, blob in sorted(staged.items()):
        try:
            text = blob.decode("utf-8")
        except UnicodeDecodeError:
            text = blob.decode("latin-1", errors="replace")
        for label, pattern in rules:
            hits = len(pattern.findall(text))
            if hits:
                per_rule[label] += hits
                total += hits
                if not quiet:
                    log(f"  HIT    {rel}: {hits} x '{label}'")
        for pattern in extra:
            hits = len(pattern.findall(text))
            if hits:
                extra_total += hits
                log(f"  EXTRA  {rel}: {hits} x /{pattern.pattern}/")
    return per_rule, total, extra_total


# --- packaging ----------------------------------------------------------------


def make_zip(staged: dict[str, bytes], out: Path, when: _dt.datetime) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    date_time = when.timetuple()[:6]
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for rel in sorted(staged):
            info = zipfile.ZipInfo(filename=rel, date_time=date_time)
            info.external_attr = 0o644 << 16   # one neutral mode, no owner/group
            info.create_system = 3
            info.create_version = 20
            zf.writestr(info, staged[rel])
    return len(staged)


def notes_text(meta: dict) -> str:
    return f"""# Notes on this peer-review package

Generated: {meta['generated_utc']}
Source: `{meta['source_description']}`
Files in package: {meta['file_count']}

## What was removed

The tree was exported from git and then scrubbed. Only the exported copy was written
to; the repository itself was not modified. Removed or withheld:

* the author name (Latin and Han forms), the e-mail address, the institution token, the
  ORCID iD and the ORCID profile link, wherever they appear -- the README author and
  contact sections, the BibTeX `author=`/`orcid=` fields, the clone line,
  `pyproject.toml` (`authors`, and the `[project.urls]` Homepage/Issues block), the
  version module's `__author__`/`__email__`/`__orcid__` values (attribute names kept,
  because the banner imports them), the package docstring, and the LICENSE copyright
  line;
* every URL pointing at the author's own repository;
* local absolute paths (home-directory, temporary-directory and Windows
  drive-qualified forms), replaced by `<local-path>`;
* `.git` -- and with it commit-author names, e-mails and the remote URL -- plus
  `__pycache__`, `.mypy_cache`, `.pytest_cache`, `.DS_Store`, `*.pyc`, build and
  egg-info output, and the third-party data directories.

The identifying strings were **derived from the repository at build time** (git commit
identity, the remote URL's owner segment, and the version module's own authorship
attributes). They are never hard-coded in `scripts/make_review_package.py`, so shipping
that script cannot re-leak what it removes.

## What was deliberately kept

* Links to third-party repositories (ASTRAL, CheckM, GTDB-Tk, HMMER and so on). They
  name other people, not the authors, and a reviewer needs them; masking them would
  cost information and buy no anonymity.
* The word "ORCID" as a field label, and the `__orcid__` attribute name. The iD is what
  identifies a person, and the iD is gone.

A reviewer-side check that greps for the bare host string `github.com` or for the word
`orcid` will still find those retained references. `--audit-extra REGEX` reports such
counts at build time without gating on them, so the numbers are known before sending.

## What the builder could not know

* **Real DOIs.** No Zenodo or journal DOI for this software exists in the repository;
  the benchmark contract test uses a deliberately synthetic placeholder identifier.
  Only the author knows the deposited record's version and concept DOIs, so nothing was
  substituted here and nothing was invented.
* **Which files were AI-assisted.** `ai-annotation-checklist.md` defines the marker line
  and ships an empty author-completed list; `scripts/add_llm_markers.py` refuses to
  guess which files those are.
* **The worked-example run outputs.** `worked_example/README.md` holds the reproduction
  recipe and no results, because the repository holds no results to report.
* **Whether the archived validation evidence matches this build.** The
  `validation/results/` artefacts were recorded on the previous tag; re-running
  `validation/` on this one is an author action.

## Reproducing HMM marker mode

`db/gtdb_markers/` is git-ignored: the TIGRFAM/Pfam profiles are third-party work
assembled locally by `scripts/fetch_marker_db.py`. {meta['hmm_note']}

## Audit

`identifying-string hits: {meta['audit_hits']}` against the derived rule set. The
builder exits non-zero when that is not 0, so a non-anonymous package cannot be
produced by accident.

The zero is not vacuous: the same rule set found **{meta['before_hits']}** identifier
occurrences in the pre-scrub text of the files being packaged. If that number were 0,
the final 0 would prove nothing -- the builder prints a warning in that case and marks
`audit_is_vacuous: true` in the manifest. Per-rule detail is in
`REVIEW_PACKAGE_MANIFEST.json`.
"""


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="make_review_package.py",
        description="Build an anonymous peer-review package (standard library only).",
    )
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1],
                        help="repository root (default: the repo this script lives in)")
    parser.add_argument("--ref", default="HEAD",
                        help="git ref exported with 'git archive' (default: HEAD)")
    parser.add_argument("--worktree", action="store_true",
                        help="export tracked + untracked-not-ignored files from disk")
    parser.add_argument("--out", type=Path, required=True, help="output .zip path")
    parser.add_argument("--hmm-dir", type=Path, default=None,
                        help="copy this marker-HMM directory in as db/gtdb_markers/ "
                             "and record its SHA-256")
    parser.add_argument("--audit-extra", action="append", default=[], metavar="REGEX",
                        help="extra pattern to report (repeatable); does not gate")
    parser.add_argument("--allow-dirty", action="store_true",
                        help="with --ref, do not note uncommitted changes")
    parser.add_argument("--keep-staging", action="store_true",
                        help="leave the scrubbed tree beside the zip for inspection")
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    if not (repo / ".git").exists():
        raise SystemExit(f"{repo} is not a git repository")

    rules = identifier_rules(repo, repo / "markerfinder" / "_version.py")
    tokens = derived_tokens(repo, repo / "markerfinder" / "_version.py")
    log(f"derived  {len(tokens)} identifier token(s) from git metadata and "
        "_version.py; rule set: "
        + ", ".join(sorted({label for label, _ in rules})))

    dirty = run_git(repo, "status", "--porcelain")
    if dirty and not args.worktree and not args.allow_dirty:
        log(f"NOTE     working tree has {len(dirty.splitlines())} uncommitted path(s); "
            "exporting the committed state of --ref. Use --worktree to include them.")

    if args.worktree:
        collected, source_desc = files_from_worktree(repo), f"worktree of {repo.name}"
    else:
        collected, source_desc = files_from_ref(repo, args.ref), f"git archive of {args.ref}"
    commit = run_git(repo, "rev-parse", "HEAD")

    scrub = Scrubber(rules, tokens)
    staged: dict[str, bytes] = {}
    raw: dict[str, bytes] = {}
    dropped_n = 0
    for rel, blob in collected:
        if dropped(rel):
            dropped_n += 1
            continue
        if is_text(rel, blob):
            try:
                raw[rel] = blob
                staged[rel] = scrub(rel, blob.decode("utf-8")).encode("utf-8")
            except UnicodeDecodeError:
                staged[rel] = blob
        else:
            staged[rel] = blob
    log(f"export   {len(staged)} file(s) staged, {dropped_n} dropped by exclusion list")

    # Positive control: a zero after scrubbing only means something if the same rules
    # found identifiers before scrubbing. Without this, "identifying-string hits: 0"
    # would be indistinguishable from a rule set that matches nothing at all.
    _, before_hits, _ = audit(raw, rules, quiet=True)
    log(f"control  {before_hits} identifier occurrence(s) in the pre-scrub text")
    if before_hits == 0:
        log("WARNING  the rule set matched nothing before scrubbing, so the final zero "
            "proves nothing; check that git metadata and the version module expose the "
            "authorship you expect (or pass --worktree).")

    hmm: dict[str, object] = {"provided": False}
    if args.hmm_dir is not None:
        hmm = install_hmm_dir(repo, args.hmm_dir, staged)
        note = (f"Its {hmm['files']} files were copied in as `db/gtdb_markers/`, and "
                f"`db/gtdb_markers.CHECKSUMS.txt` records every per-file SHA-256 plus "
                f"the aggregate `{hmm['aggregate_sha256']}`, computed the way "
                f"`db/expected_hashes.json` is, so the library can be verified rather "
                f"than trusted.")
        log(f"hmm-dir  {hmm['files']} file(s), aggregate sha256 {hmm['aggregate_sha256']}")
    else:
        note = ("**No HMM library was supplied to this build**, so `--marker-mode hmm` "
                "-- and with it the worked example -- cannot be reproduced from the "
                "package alone. Run `python scripts/fetch_marker_db.py --auto`, or "
                "rebuild with `--hmm-dir db/gtdb_markers`.")

    extras = [re.compile(p) for p in args.audit_extra]
    per_rule, hits, extra_hits = audit(staged, rules, extras)
    log(f"audit    {len(rules)} rule(s) over {len(staged)} file(s)")
    log(f"identifying-string hits (sources): {hits}")
    if extra_hits:
        log(f"extra-pattern hits: {extra_hits} (reported, not gating)")

    now = _dt.datetime.now(_dt.timezone.utc)
    stamp = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    staged["REVIEW_PACKAGE_NOTES.md"] = notes_text({
        "generated_utc": stamp, "source_description": source_desc,
        "file_count": len(staged) + 2, "hmm_note": note,
        "audit_hits": hits, "before_hits": before_hits,
    }).encode("utf-8")

    manifest = {
        "generated_utc": stamp,
        "source": source_desc,
        "head_commit": commit,
        "ref": None if args.worktree else args.ref,
        "files_in_package": len(staged),
        "dropped_by_exclusion": dropped_n,
        "rules_checked": sorted({label for label, _ in rules}),
        "identifier_occurrences_before_scrubbing": before_hits,
        "identifying_string_hits": hits,
        "audit_is_vacuous": before_hits == 0,
        "per_rule_hits": per_rule,
        "extra_pattern_hits": extra_hits,
        "hmm_library": hmm,
        "note": ("Author name/e-mail/institution/ORCID, own-repository URLs and local "
                 "absolute paths were scrubbed in this exported copy only. DOIs and "
                 "AI-annotation targets are not invented; see REVIEW_PACKAGE_NOTES.md."),
    }
    staged["REVIEW_PACKAGE_MANIFEST.json"] = (json.dumps(
        manifest, indent=2, ensure_ascii=False) + "\n").encode("utf-8")

    # Final pass: the two generated files are part of the package too.
    per_rule, hits, _ = audit(staged, rules, [])
    log(f"identifying-string hits: {hits}")

    count = make_zip(staged, args.out, now)
    log(f"zip      {args.out} :: {count} file(s), {args.out.stat().st_size:,} bytes")

    if args.keep_staging:
        target = args.out.with_suffix(".staged")
        for rel, blob in staged.items():
            path = target / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(blob)
        log(f"staging  kept at {target}")

    if hits:
        log("FAIL     the package still carries author identifiers; do not send it")
        return 3
    log("ok       package is anonymous against the derived rule set")
    return 0


if __name__ == "__main__":
    sys.exit(main())
