"""The anonymity packager (``scripts/make_review_package.py``) must scrub, prove and fail closed.

The manuscript tells reviewers the code reaches them as "an anonymised package
(version-control history and identifying metadata removed)". That claim is executed by
this script, which until now had no test: a green ``identifying-string hits: 0`` could
have come from a rule set that matched nothing at all.

Cases run against a scratch git repository in ``tmp_path`` whose authorship is
fictitious, so the real end-to-end path (collect -> scrub -> audit -> zip -> exit code)
is exercised without touching this repository. They pin down:

(a) every identifier the rule set derives is gone from the shipped zip, including a
    surname written on its own (the audit matches parts, so removal must too);
(b) the positive control is real -- the same rules found those identifiers *before*
    scrubbing, and the manifest records it;
(c) the tool fails closed (exit 3) when a leak survives;
(d) the documented blind spot (X-045): when the repository stores the name as one
    concatenated word, a two-word rendering of it matches no derived rule, the gate
    stays green while the variant ships, and only ``--audit-extra`` surfaces it.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "make_review_package.py"

NAME_SPACED = "Testperson Testfamily"
NAME_ONEWORD = "TestpersonTestfamily"
# Assembled at runtime, never written out whole. This file is itself inside the tree the
# packager scrubs, so a literal e-mail address, ORCID iD or absolute path here would be
# rewritten in the shipped copy — the fixture would then build a fake repo carrying no
# identity, the leak gate would correctly fire, and `pytest tests/unit` would fail for the
# one person who is supposed to be able to run it: the reviewer.
EMAIL = "testperson" + "@" + "institute.example"
ORCID = "-".join(["0000", "0001", "2345", "6789"])
REPO_URL = "https://github.com/testperson/markerfinder"
ABS_PATH = "/home" + "/testperson/work/markerfinder/run.py"

GIT_ENV = dict(GIT_CONFIG_GLOBAL="/dev/null", GIT_CONFIG_SYSTEM="/dev/null",
               HOME="/tmp", PATH=os.environ.get("PATH", ""))


def _load():
    spec = importlib.util.spec_from_file_location("make_review_package_under_test", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _make_repo(tmp_path: Path, author_name: str) -> Path:
    """A committed scratch repository carrying a full set of author identifiers."""
    repo = tmp_path / "repo"
    (repo / "markerfinder").mkdir(parents=True)
    (repo / "markerfinder" / "_version.py").write_text(
        '"""Version and authorship."""\n'
        '__version__ = "9.9.9"\n'
        f'__author__ = "{author_name}"\n'
        f'__email__ = "{EMAIL}"\n'
        f'__orcid__ = "{ORCID}"\n', encoding="utf-8")
    (repo / "markerfinder" / "__init__.py").write_text(
        f'"""MarkerFinder by {author_name} <{EMAIL}>, see {REPO_URL}."""\n',
        encoding="utf-8")
    (repo / "README.md").write_text(
        f"# MarkerFinder\n\nAuthor: {author_name} ({EMAIL}, ORCID {ORCID})\n\n"
        f"Clone: git clone {REPO_URL}\nDebug with: {ABS_PATH}\n", encoding="utf-8")
    (repo / ".gitignore").write_text("__pycache__/\n*.pyc\nbuild/\n", encoding="utf-8")
    (repo / "build").mkdir()
    (repo / "build" / "leftover.py").write_text(f"# {author_name}\n", encoding="utf-8")
    (repo / "markerfinder" / "__pycache__").mkdir()
    (repo / "markerfinder" / "__pycache__" / "x.pyc").write_bytes(b"\x00\x01fake")

    def git(*args):
        return subprocess.run(["git", "-C", str(repo), *args], capture_output=True,
                              text=True, env=GIT_ENV)

    ident = ["-c", f"user.name={author_name}", "-c", f"user.email={EMAIL}"]
    assert git("init", "-q").returncode == 0
    assert git(*ident, "add", "-A").returncode == 0
    committed = git(*ident, "-c", "commit.gpgsign=false", "commit", "-q", "-m", "seed")
    assert committed.returncode == 0, committed.stderr
    git("remote", "add", "origin", REPO_URL)
    return repo


@pytest.fixture()
def fake_repo(tmp_path):
    """Metadata stores the name as two words, so both parts become rules."""
    return _make_repo(tmp_path, NAME_SPACED)


@pytest.fixture()
def oneword_repo(tmp_path):
    """Metadata stores the name as one concatenated word (this repository's shape)."""
    return _make_repo(tmp_path, NAME_ONEWORD)


def _run(module, repo, out, extra=()):
    return module.main(["--repo", str(repo), "--worktree", "--out", str(out), *extra])


def _members(out):
    with zipfile.ZipFile(out) as z:
        return {n: z.read(n) for n in z.namelist()}


def test_every_derived_identifier_is_scrubbed(fake_repo, tmp_path):
    out = tmp_path / "pkg.zip"
    assert _run(_load(), fake_repo, out) == 0
    text = "\n".join(b.decode("utf-8", "replace") for n, b in _members(out).items())
    for secret in (NAME_SPACED, EMAIL, ORCID, REPO_URL, ABS_PATH):
        assert secret not in text, f"{secret} survived into the shipped package"
    assert "markerfinder/_version.py" in _members(out)
    assert '__author__ = ""' in text  # attribute kept for the banner, value withheld


def test_surname_written_alone_is_scrubbed_too(fake_repo, tmp_path):
    """Detection matches token parts, so removal has to as well."""
    (fake_repo / "CITATION.md").write_text(
        "Cite as: Testfamily et al., MarkerFinder.\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(fake_repo), "add", "-A"], env=GIT_ENV,
                   capture_output=True)
    out = tmp_path / "pkg.zip"
    module = _load()
    assert _run(module, fake_repo, out) == 0
    assert "Testfamily" not in _members(out)["CITATION.md"].decode()


def test_positive_control_proves_the_final_zero_is_meaningful(fake_repo, tmp_path):
    out = tmp_path / "pkg.zip"
    assert _run(_load(), fake_repo, out) == 0
    manifest = json.loads(_members(out)["REVIEW_PACKAGE_MANIFEST.json"].decode())
    assert manifest["identifier_occurrences_before_scrubbing"] > 0
    assert manifest["audit_is_vacuous"] is False
    assert manifest["identifying_string_hits"] == 0


def test_vcs_and_build_artifacts_never_ship(fake_repo, tmp_path):
    out = tmp_path / "pkg.zip"
    assert _run(_load(), fake_repo, out) == 0
    names = _members(out)
    assert not [n for n in names if n.startswith(".git/") or "__pycache__" in n
                or n.startswith("build/") or n.endswith(".pyc")]


def test_tool_fails_closed_when_a_leak_survives(fake_repo, tmp_path, monkeypatch):
    module = _load()

    class NoScrub:  # every rule matches, nothing is removed
        def __init__(self, rules, tokens):
            pass

        def __call__(self, rel, text):
            return text

    monkeypatch.setattr(module, "Scrubber", NoScrub)
    assert module.main(["--repo", str(fake_repo), "--worktree",
                        "--out", str(tmp_path / "pkg.zip")]) == 3


def test_concatenated_metadata_cannot_see_a_spaced_variant_until_audit_extra(
        oneword_repo, tmp_path, capsys):
    """X-045: a rule exists only for a spelling the repository actually uses.

    With the name stored as ``TestpersonTestfamily``, no rule covers the surname on its
    own, so a two-word rendering is only half-masked -- the e-mail user token happens to
    match the given name, and ``Testfamily`` ships untouched -- while the gate still
    reports zero. ``--audit-extra`` is the escape hatch that surfaces it, and it is why
    the shipped package is re-scanned independently of this tool's own verdict.
    """
    (oneword_repo / "NOTES.md").write_text(
        "The author signs as Testperson Testfamily.\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(oneword_repo), "add", "-A"], env=GIT_ENV,
                   capture_output=True)
    out = tmp_path / "pkg.zip"
    module = _load()
    assert _run(module, oneword_repo, out) == 0, "the gate should be blind here"
    body = _members(out)["NOTES.md"].decode()
    assert "Testfamily" in body, "expected the surname to escape the derived rules"
    assert "Testperson Testfamily" not in body, "given name is masked via the e-mail user"

    capsys.readouterr()
    assert _run(module, oneword_repo, out, extra=("--audit-extra", "Testfamily")) == 0
    assert "extra-pattern hits: 1" in capsys.readouterr().out
