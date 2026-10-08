#!/usr/bin/env python3
"""Apply and verify the repository's LLM-assistance annotation convention.

Why this script exists
----------------------
``AI_ASSISTANCE.md`` defines one annotation line: the comment reading
"LLM-assisted revision; reviewed against tests", pointing back to that file. Its
exact text lives in the ``MARKER`` constant below and is quoted nowhere else in this
file, so the tool cannot mistake its own source for an annotated module.

``AI_ASSISTANCE.md`` also says that *which files the assistant touched is a factual
question about the project history, and only the author can answer it*. That sentence
is the design constraint here: this tool never inspects code to decide what deserves a
marker, never infers a target from git history, and has no "auto" mode. It takes an
explicit, author-authored list of targets and does two things with it -- put the
marker on those targets, and check that the marker is on those targets.

The list format
---------------
The list lives in ``ai-annotation-checklist.md`` by default, or in any file named by
``--list``. One target per line::

    markerfinder/modules/hgt_filter.py            # whole module
    markerfinder/pipeline.py::run_pipeline        # one function only
    # a comment line, ignored
    markerfinder/models/marker.py::ResolutionPreset

Bare paths are placed immediately after the module docstring (at the top of the file
when there is none). ``path::name`` targets are placed immediately above the named
function or class, ahead of its decorators, at the definition's own indentation.

Commands
--------
    apply     write the marker onto every listed target (default action)
    check     verify every listed target carries the marker, and that no file outside
              the list carries it; exit non-zero on either
    inventory print which files currently carry the marker

``apply --dry-run`` shows what would change and writes nothing.
``check --require-nonempty`` turns an empty list into a failure, so a CI job cannot
go green while the author still has not filled the list in.

Exit codes: 0 ok; 1 a marker is missing or an unsanctioned marker was found; 2 the
list itself is unusable (missing file, a path outside the repository, a missing
target file, or a named symbol that does not exist).
"""

from __future__ import annotations

import argparse
import ast
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator, NamedTuple, Optional

#: The project root every path is resolved against. ``--repo`` overrides it so the
#: tool (and its test) can operate on a scratch copy without touching a live tree.
REPO = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()
DEFAULT_LIST_NAME = "ai-annotation-checklist.md"

#: The one line the convention defines. It appears exactly once in this file, which
#: is why the docstring above describes it in prose instead of quoting it.
MARKER = "# LLM-assisted revision; reviewed against tests (see AI_ASSISTANCE.md)"

#: Directories that are not part of a release and so are not annotation candidates.
SKIP_DIR_PARTS = {".git", "__pycache__", ".mypy_cache", ".pytest_cache", ".work",
                  "node_modules", ".venv", "venv", "build", "dist"}

_TARGET = re.compile(r"^(?P<path>[^\s#]+?\.py)(?:::(?P<symbol>[A-Za-z_][\w.]*))?"
                     r"\s*(?:#.*)?$")


class ListError(RuntimeError):
    """The author's list is unusable -- refuse to proceed rather than guess."""


class Target(NamedTuple):
    path: Path
    symbol: Optional[str]
    line_no: int  # line in the list file, so an error can be pointed at it

    @property
    def key(self) -> str:
        rel = str(self.path.relative_to(REPO))
        return rel + (f"::{self.symbol}" if self.symbol else "")


class Placement(NamedTuple):
    index: int           # line index (0-based) before which the marker goes
    indent: str          # indentation keeping the marker in the right block
    reason: str          # printed justification
    module_level: bool


def set_root(root: Path) -> None:
    global REPO
    REPO = root.resolve()


def _is_target_line(text: str) -> bool:
    """Skip document prose, tables and headings; accept only real target lines.

    The list is allowed to live inside a readable Markdown checklist, so only lines
    that really look like ``path.py`` or ``path.py::symbol`` count as targets.
    """
    if "|" in text or text.startswith(("#", ">", "-", "*", "`", "[")):
        return False
    return bool(_TARGET.match(text))


def parse_list(list_path: Path) -> list[Target]:
    """Read the author-supplied target list. Never invents a target."""
    if not list_path.exists():
        raise ListError(
            f"annotation list not found: {list_path}\n"
            "This tool refuses to guess which files an AI assistant touched. Fill the "
            "list in (one path per line, optionally path::symbol) and re-run."
        )
    targets: list[Target] = []
    seen: set[str] = set()
    for line_no, raw in enumerate(list_path.read_text(encoding="utf-8").splitlines(), 1):
        text = raw.strip()
        if not text or not _is_target_line(text):
            continue
        match = _TARGET.match(text)
        assert match is not None
        rel = match.group("path")
        candidate = (REPO / rel).resolve()
        try:
            candidate.relative_to(REPO)
        except ValueError:
            raise ListError(
                f"list line {line_no}: '{rel}' resolves outside the repository "
                f"({candidate}). Refusing to write outside the project tree."
            ) from None
        if not candidate.is_file():
            raise ListError(
                f"list line {line_no}: '{rel}' does not exist. The list is a factual "
                "claim about the project history, so a stale or mistyped path has to "
                "be corrected, not skipped."
            )
        target = Target(candidate, match.group("symbol"), line_no)
        if target.key in seen:
            continue
        seen.add(target.key)
        targets.append(target)
    return targets


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _has_marker(lines: Iterable[str]) -> bool:
    return any(line.strip() == MARKER for line in lines)


def _docstring_end(tree: ast.Module) -> int:
    """Line index (0-based) just after the module docstring; 0 when there is none."""
    body = tree.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
            and isinstance(body[0].value.value, str):
        return int(body[0].end_lineno)
    return 0


def find_symbol(tree: ast.Module, dotted: str) -> Optional[ast.AST]:
    """Locate a function or class by (possibly dotted) name, searching top-down."""
    nodes: list[ast.AST] = list(tree.body)
    found: Optional[ast.AST] = None
    for name in dotted.split("."):
        found = None
        for node in nodes:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) \
                    and node.name == name:
                found = node
                break
        if found is None:
            return None
        nodes = list(getattr(found, "body", []))
    return found


def placement_for(path: Path, symbol: Optional[str]) -> Placement:
    """Where the marker line belongs for this target."""
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    tree = _tree(path)
    if symbol is None:
        idx = _docstring_end(tree)
        reason = ("module level, after the module docstring" if idx
                  else "module level, top of file")
        return Placement(idx, "", reason, True)
    node = find_symbol(tree, symbol)
    if node is None:
        raise ListError(
            f"{path.relative_to(REPO)}: symbol '{symbol}' not found. Refusing to mark "
            "a location it cannot locate, and refusing to fall back to the module."
        )
    start = int(node.lineno) - 1
    decorators = list(getattr(node, "decorator_list", []))
    if decorators:
        start = min(start, int(decorators[0].lineno) - 1)
    line = lines[start]
    indent = line[: len(line) - len(line.lstrip())]
    kind = "class" if isinstance(node, ast.ClassDef) else "function"
    return Placement(start, indent, f"above the {kind} '{symbol}'", False)


def marker_present(path: Path, placement: Placement) -> bool:
    """Is this target already annotated? Idempotency depends on this being exact.

    After a successful insertion the parsed line numbers shift by one, so the marker
    ends up at ``index - 1`` rather than at ``index``. Both positions are checked,
    which is what makes a second run a no-op instead of a duplicate.
    """
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    window = lines[max(0, placement.index - 1): placement.index + 1]
    if _has_marker(window):
        return True
    if placement.module_level:
        # A module marker may also sit above the docstring, or at line 0.
        return _has_marker(lines[: placement.index + 1])
    return False


def insert(path: Path, placement: Placement) -> None:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    if lines and not lines[-1].endswith("\n"):
        lines.append("\n")
    lines.insert(placement.index, f"{placement.indent}{MARKER}\n")
    path.write_text("".join(lines), encoding="utf-8")


def tracked_python_files() -> Iterator[Path]:
    """Python files the convention could apply to."""
    for path in sorted(REPO.rglob("*.py")):
        if any(part in SKIP_DIR_PARTS for part in path.parts):
            continue
        if path.resolve() == SELF:
            # This tool documents the marker; it is not itself an assisted module.
            continue
        yield path


def marked_files() -> list[Path]:
    """Files that already carry the marker, for ``inventory`` and the drift check."""
    out: list[Path] = []
    for path in tracked_python_files():
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeDecodeError):
            continue
        if _has_marker(lines):
            out.append(path)
    return out


def _syntax_err(exc: SyntaxError) -> str:
    return f"{exc.filename}: {exc.msg} (line {exc.lineno})"


def cmd_apply(targets: list[Target], dry_run: bool, verbose: bool) -> int:
    changed = kept = errors = 0
    for target in targets:
        try:
            placement = placement_for(target.path, target.symbol)
        except ListError as exc:
            print(f"ERROR  {exc}", file=sys.stderr)
            errors += 1
            continue
        except SyntaxError as exc:
            print(f"ERROR  {target.key}: {_syntax_err(exc)}", file=sys.stderr)
            errors += 1
            continue
        if marker_present(target.path, placement):
            kept += 1
            if verbose:
                print(f"ok     {target.key} (marker already present)")
            continue
        changed += 1
        if dry_run:
            print(f"would  mark {target.key} -- {placement.reason}")
            continue
        insert(target.path, placement)
        print(f"marked {target.key} -- {placement.reason}")
    print(f"{'[dry-run] ' if dry_run else ''}{changed} to add, {kept} already marked, "
          f"{errors} error(s), {len(targets)} listed target(s)")
    return 2 if errors else 0


def cmd_check(targets: list[Target], require_nonempty: bool, verbose: bool) -> int:
    """The CI half of the convention: the markers and the list must agree."""
    if not targets:
        msg = ("annotation list declares 0 targets -- nothing is claimed, so nothing "
               "can be verified. Only the author knows which files were AI-assisted; "
               f"fill in {DEFAULT_LIST_NAME}.")
        if require_nonempty:
            print(f"FAIL   {msg}", file=sys.stderr)
            return 1
        print(f"check  {msg}")
        return 0

    missing = 0
    for target in targets:
        try:
            placement = placement_for(target.path, target.symbol)
        except ListError as exc:
            print(f"ERROR  {exc}", file=sys.stderr)
            return 2
        if marker_present(target.path, placement):
            if verbose:
                print(f"ok     {target.key}")
            continue
        missing += 1
        print(f"MISSING {target.key} -- no marker line where the list says one belongs "
              f"({placement.reason})", file=sys.stderr)

    sanctioned = {t.path for t in targets}
    stray = [p for p in marked_files() if p not in sanctioned]
    for path in stray:
        print(f"UNSANCTIONED {path.relative_to(REPO)} carries the marker, but the "
              "annotation list does not declare it", file=sys.stderr)

    if missing or stray:
        print(f"check    {missing} missing, {len(stray)} unsanctioned of "
              f"{len(targets)} listed target(s) -- FAIL", file=sys.stderr)
        return 1
    print(f"check    {len(targets)} listed target(s) all marked, 0 unsanctioned -- OK")
    return 0


def cmd_inventory() -> int:
    found = marked_files()
    for path in found:
        print(path.relative_to(REPO))
    print(f"{len(found)} marked file(s)")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="add_llm_markers.py",
        description="Apply/verify the LLM-assistance marker from an author-supplied "
                    "list. Refuses to guess which files were assisted.",
    )
    parser.add_argument("command", nargs="?", default="apply",
                        choices=["apply", "check", "inventory"],
                        help="apply (default), check (CI gate), or inventory")
    parser.add_argument("--repo", type=Path, default=None,
                        help="project root targets resolve against (default: the "
                             "repository this script lives in)")
    parser.add_argument("--list", dest="list_path", type=Path, default=None,
                        help=f"author-supplied target list (default: "
                             f"<repo>/{DEFAULT_LIST_NAME})")
    parser.add_argument("--dry-run", action="store_true",
                        help="with apply: print what would change, write nothing")
    parser.add_argument("--check", action="store_true",
                        help="run the check command (same as the 'check' positional)")
    parser.add_argument("--require-nonempty", action="store_true",
                        help="with check: fail when the list declares no target")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    if args.repo is not None:
        if not args.repo.is_dir():
            print(f"ERROR  --repo {args.repo} is not a directory", file=sys.stderr)
            return 2
        set_root(args.repo)

    list_path = args.list_path if args.list_path is not None else REPO / DEFAULT_LIST_NAME

    command = "check" if args.check else args.command
    if command == "inventory":
        return cmd_inventory()
    try:
        targets = parse_list(list_path)
    except ListError as exc:
        print(f"ERROR  {exc}", file=sys.stderr)
        return 2
    if command == "check":
        return cmd_check(targets, args.require_nonempty, args.verbose)
    return cmd_apply(targets, args.dry_run, args.verbose)


if __name__ == "__main__":
    sys.exit(main())
