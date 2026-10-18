"""The LLM-annotation mechanism (``scripts/add_llm_markers.py``) must mark, hold, and refuse.

``AI_ASSISTANCE.md`` states a convention -- one marker line, placed after the module
docstring or above an affected function -- and pointed at two artefacts that did not
exist. This test is the mechanism half of that gap. It proves the script:
(a) puts the marker where the convention says;
(b) is idempotent, so re-running is safe;
(c) gates CI in *both* directions -- a declared marker missing, and a marker nobody
    declared, are each a failure;
(d) refuses to guess: an absent list, a stale path, or a missing symbol is an error,
    never a fallback to marking the whole module.

Every case runs against a scratch tree in ``tmp_path`` via ``--repo``, never against
this repository's own sources. The autouse fixture restores the module's root so no
other test inherits a scratch tree, and keeps ``__pycache__`` out of the repository.
"""
from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "add_llm_markers.py"
DOCSUB = '"""Docstring line.\n\nMore docstring.\n"""\n\nimport os\n\n\ndef f():\n    return os\n'
NODOCSUB = "import os\n\n\ndef g():\n    return os\n"


def _load():
    spec = importlib.util.spec_from_file_location("add_llm_markers_under_test", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


mod = _load()
MARKER = mod.MARKER


@pytest.fixture(autouse=True)
def _restore_root():
    original = mod.REPO
    yield
    mod.REPO = original


@pytest.fixture
def scratch(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "dont_write_bytecode", True)
    mod.set_root(tmp_path)
    return tmp_path


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def run(root: Path, *args: str) -> int:
    return mod.main(["--repo", str(root), *args])


# -- apply: placement follows the convention ---------------------------------


def test_apply_places_marker_after_module_docstring(scratch):
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    assert run(scratch, "apply", "--list", str(lst)) == 0
    lines = target.read_text(encoding="utf-8").splitlines()
    # the docstring spans lines 1-4 (its closing quotes are line 4), so the marker
    # goes on line 5 -- immediately after it, and before the first import
    assert lines[3] == '"""'
    assert lines[4].strip() == MARKER
    assert lines[5] == ""
    assert lines[6] == "import os"


def test_apply_places_marker_at_top_when_there_is_no_docstring(scratch):
    target = write(scratch / "pkg" / "b.py", NODOCSUB)
    lst = write(scratch / "list.txt", "pkg/b.py\n")
    run(scratch, "apply", "--list", str(lst))
    assert target.read_text(encoding="utf-8").splitlines()[0].strip() == MARKER


def test_apply_places_marker_above_named_function_ahead_of_decorators(scratch):
    src = '"""Doc."""\n\n\nclass C:\n\n    @property\n    def meth(self):\n        return 1\n'
    target = write(scratch / "pkg" / "c.py", src)
    lst = write(scratch / "list.txt", "pkg/c.py::C.meth\n")
    assert run(scratch, "apply", "--list", str(lst)) == 0
    lines = target.read_text(encoding="utf-8").splitlines()
    at = [i for i, ln in enumerate(lines) if ln.strip() == MARKER]
    assert len(at) == 1, lines
    assert lines[at[0]] == "    " + MARKER, "marker must carry the definition's indent"
    assert lines[at[0] + 1].strip() == "@property"
    assert lines[at[0] + 2].strip() == "def meth(self):"


def test_apply_is_idempotent_for_module_targets(scratch):
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    run(scratch, "apply", "--list", str(lst))
    first = target.read_text(encoding="utf-8")
    assert run(scratch, "apply", "--list", str(lst)) == 0
    assert target.read_text(encoding="utf-8") == first, "second run changed the file"
    assert first.count(MARKER) == 1


def test_apply_is_idempotent_for_symbol_targets(scratch):
    """Line numbers shift once inserted; the marker must not be added twice."""
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py::f\n")
    run(scratch, "apply", "--list", str(lst))
    first = target.read_text(encoding="utf-8")
    run(scratch, "apply", "--list", str(lst))
    assert target.read_text(encoding="utf-8") == first
    assert first.count(MARKER) == 1


def test_apply_touches_no_other_file(scratch):
    write(scratch / "pkg" / "a.py", DOCSUB)
    other = write(scratch / "pkg" / "z.py", NODOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    run(scratch, "apply", "--list", str(lst))
    assert other.read_text(encoding="utf-8") == NODOCSUB


def test_dry_run_writes_nothing(scratch):
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    assert run(scratch, "apply", "--dry-run", "--list", str(lst)) == 0
    assert target.read_text(encoding="utf-8") == DOCSUB


def test_apply_preserves_trailing_newline(scratch):
    target = write(scratch / "pkg" / "n.py", '"""Doc."""\nx = 1')  # no final newline
    lst = write(scratch / "list.txt", "pkg/n.py\n")
    run(scratch, "apply", "--list", str(lst))
    assert target.read_text(encoding="utf-8").endswith("\n")


def test_result_still_compiles(scratch):
    target = write(scratch / "pkg" / "c.py",
                   '"""Doc."""\n\n\nclass C:\n    def m(self):\n        return 1\n')
    lst = write(scratch / "list.txt", "pkg/c.py::C.m\n")
    run(scratch, "apply", "--list", str(lst))
    compile(target.read_text(encoding="utf-8"), str(target), "exec")


# -- refusing to guess --------------------------------------------------------


def test_missing_list_is_an_error_not_a_guess(scratch):
    assert run(scratch, "check", "--list", str(scratch / "nope.txt")) == 2


def test_unknown_path_is_refused(scratch):
    write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/ghost.py\n")
    assert run(scratch, "apply", "--list", str(lst)) == 2


def test_unknown_symbol_is_refused_without_falling_back_to_module(scratch):
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py::not_there\n")
    assert run(scratch, "apply", "--list", str(lst)) == 2
    assert target.read_text(encoding="utf-8") == DOCSUB, \
        "the tool marked the whole module after failing to locate the symbol"


def test_path_resolving_outside_the_repo_is_refused(scratch):
    lst = write(scratch / "list.txt", "../escape.py\n")
    assert run(scratch, "apply", "--list", str(lst)) == 2


def test_markdown_prose_and_table_lines_are_not_targets(scratch):
    """The list lives inside a readable document; only path lines may ever mark."""
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt",
                "# heading\n\n"
                "| File | Assisted? |\n|---|---|\n| pkg/a.py | yes |\n\n"
                "- a bullet\n"
                "> a quote naming pkg/a.py\n"
                "`pkg/a.py` inline code\n")
    assert run(scratch, "check", "--list", str(lst)) == 0
    assert target.read_text(encoding="utf-8") == DOCSUB


def test_trailing_comment_on_a_target_line_is_ignored(scratch):
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py   # whole module, reviewed 2026-10\n")
    run(scratch, "apply", "--list", str(lst))
    assert MARKER in target.read_text(encoding="utf-8")


def test_duplicate_target_is_applied_once(scratch):
    target = write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\npkg/a.py\n")
    run(scratch, "apply", "--list", str(lst))
    assert target.read_text(encoding="utf-8").count(MARKER) == 1


# -- check: the CI gate, in both directions ----------------------------------


def test_check_passes_after_apply(scratch):
    write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    run(scratch, "apply", "--list", str(lst))
    assert run(scratch, "check", "--list", str(lst)) == 0


def test_check_fails_when_a_declared_marker_is_absent(scratch):
    write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    assert run(scratch, "check", "--list", str(lst)) == 1


def test_check_fails_on_an_unsanctioned_marker(scratch):
    """A marker nobody declared is as much a defect as a missing one."""
    write(scratch / "pkg" / "a.py", DOCSUB)
    write(scratch / "pkg" / "rogue.py", DOCSUB.replace("import os", MARKER + "\nimport os"))
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    run(scratch, "apply", "--list", str(lst))
    assert run(scratch, "check", "--list", str(lst)) == 1


def test_check_passes_when_symbol_target_marker_is_present(scratch):
    write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py::f\n")
    run(scratch, "apply", "--list", str(lst))
    assert run(scratch, "check", "--list", str(lst)) == 0


def test_empty_list_passes_and_says_so(scratch, capsys):
    lst = write(scratch / "list.txt", "# nothing filled in yet\n")
    assert run(scratch, "check", "--list", str(lst)) == 0
    assert "0 targets" in capsys.readouterr().out


def test_require_nonempty_turns_an_empty_list_into_a_failure(scratch):
    lst = write(scratch / "list.txt", "# nothing filled in yet\n")
    assert run(scratch, "check", "--require-nonempty", "--list", str(lst)) == 1


# -- inventory, and the tool's freedom from its own string --------------------


def test_inventory_reports_marked_files(scratch, capsys):
    write(scratch / "pkg" / "a.py", DOCSUB)
    lst = write(scratch / "list.txt", "pkg/a.py\n")
    run(scratch, "apply", "--list", str(lst))
    capsys.readouterr()
    assert run(scratch, "inventory") == 0
    out = capsys.readouterr().out
    assert "pkg/a.py" in out
    assert "1 marked file(s)" in out


def test_inventory_is_empty_on_this_repository():
    """The tool must not read as annotated merely because it names the marker.

    Its docstring describes the line in prose. If that ever becomes a literal quote of
    it, this test fires -- inventory would report a file nobody declared, exactly the
    drift ``check`` exists to catch.
    """
    mod.set_root(REPO)
    assert mod.marked_files() == []


def test_marker_string_is_the_one_the_convention_document_defines():
    doc = (REPO / "AI_ASSISTANCE.md").read_text(encoding="utf-8")
    assert MARKER in doc, "the script's marker is not the line AI_ASSISTANCE.md states"


def test_marker_literal_appears_once_and_never_as_a_standalone_line():
    """The script defines the marker exactly once and must never *carry* it.

    A standalone marker line in this file would make the tool report itself as an
    annotated module, so ``inventory`` and ``check`` would disagree with reality.
    """
    body = SCRIPT.read_text(encoding="utf-8")
    assert body.count(MARKER) == 1, "the marker string must be defined once, not quoted"
    assert not any(ln.strip() == MARKER for ln in body.splitlines()), \
        "the script carries the marker as its own line and would mark itself"


def test_script_imports_only_the_standard_library():
    """Review tooling must run in a bare environment."""
    allowed = {"__future__", "argparse", "ast", "re", "sys", "pathlib", "typing"}
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert imported <= allowed, f"unexpected dependency: {sorted(imported - allowed)}"
