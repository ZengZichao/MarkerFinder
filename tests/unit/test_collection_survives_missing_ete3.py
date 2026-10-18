"""The suite must stay collectable when ete3 cannot be imported.

Finding: ``tests/unit/test_mad_root_regression.py`` and
``tests/unit/test_taxonomy_regression.py`` once did a bare module-level
``from ete3 import Tree``. ete3 imports the stdlib ``cgi`` module at
package-import time, so on any interpreter where that import fails, collection
aborted with

    Interrupted: 2 errors during collection
    ============================== 2 errors in 2.24s ==============================

and ``pytest tests`` exited 2 having reported **zero** test results. A suite
that cannot be collected is not a suite with skips: it is a suite that tells you
nothing, which is precisely the failure mode this file was written to
eliminate. The sanctioned pattern already used by ``test_etree_splits.py`` /
``test_static_check_contract.py`` is a labelled skip whose reason contains
``NOT EXECUTED`` so it can never be mistaken for a pass.

Three gates here, because they fail in different environments:

* the static import ratchet fires wherever it runs, so a re-introduced bare
  import is caught even on 3.10 where ete3 imports happily;
* the skip-labelling ratchet keeps every ete3-driven skip readable as
  NOT EXECUTED in a CI summary;
* the collect-only subprocess proves the property end to end on the interpreter
  that actually matters — the one where ete3 is unimportable. (On Python 3.13+
  ``markerfinder._cgi_compat`` normally makes ete3 importable again, so this
  guard is reached by simulating a missing ete3 rather than by waiting for a
  particular interpreter version.)
"""

from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

TESTS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = TESTS_DIR.parent

REQUIRE_TOKEN = "NOT EXECUTED"


# ---------------------------------------------------------------------------
# Static analysis helpers
# ---------------------------------------------------------------------------

def _module_level_ete3_imports(tree: ast.Module) -> List[int]:
    """Line numbers of module-level (not function/class-local) ete3 imports."""
    lines: List[int] = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and (node.module or "") == "ete3":
            lines.append(node.lineno)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] == "ete3":
                    lines.append(node.lineno)
    return lines


def _module_level_skip_lines(tree: ast.Module) -> List[Tuple[int, bool]]:
    """Return (lineno, reason_is_labelled) for module-level importorskip("ete3")."""
    found: List[Tuple[int, bool]] = []
    for node in tree.body:
        call = node.value if isinstance(node, ast.Expr) else None
        if not isinstance(call, ast.Call):
            continue
        func = call.func
        name = (
            func.attr
            if isinstance(func, ast.Attribute)
            else getattr(func, "id", "") if isinstance(func, ast.Name) else ""
        )
        if name != "importorskip":
            continue
        if not call.args or not isinstance(call.args[0], ast.Constant):
            continue
        if call.args[0].value != "ete3":
            continue
        strings = [
            c.value
            for c in ast.walk(call)
            if isinstance(c, ast.Constant) and isinstance(c.value, str)
        ]
        found.append((call.lineno, any(REQUIRE_TOKEN in s for s in strings)))
    return found


def find_unguarded_phylo_imports(source: str) -> List[str]:
    """Problems in ``source``: bare or unlabelled module-level ete3 imports."""
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:  # Pragma: no cover - a broken file is its own bug
        return [f"unparseable: {exc}"]
    problems: List[str] = []
    skips = _module_level_skip_lines(tree)
    for lineno in _module_level_ete3_imports(tree):
        guarding = [s for s in skips if s[0] < lineno]
        if not guarding:
            problems.append(
                f"line {lineno}: module-level ete3 import with no preceding "
                f"pytest.importorskip('ete3', ...)"
            )
        elif not any(labelled for _, labelled in guarding):
            problems.append(
                f"line {lineno}: importorskip('ete3') reason must contain "
                f"'{REQUIRE_TOKEN}' so a skip is never read as a pass"
            )
    return problems


# ---------------------------------------------------------------------------
# Gate 1 — static ratchet (interpreter independent)
# ---------------------------------------------------------------------------

def test_control_the_scanner_both_fires_and_stays_silent():
    """Self-contained control: the scanner must catch a planted defect.

    Without this, a scanner that silently returns [] for everything — because of
    a wrong attribute name, say — would make the real gate vacuously green.
    """
    bare = "import pytest\nfrom ete3 import Tree as EteTree\n"
    unlabelled = (
        "import pytest\n"
        'pytest.importorskip("ete3", reason="no ete3 here")\n'
        "from ete3 import Tree as EteTree\n"
    )
    guarded = (
        "import pytest\n"
        'pytest.importorskip("ete3", reason="ete3 missing — NOT EXECUTED")\n'
        "from ete3 import Tree as EteTree\n"
    )
    deferred = "import pytest\n\n\ndef test_x():\n    from ete3 import Tree\n"

    assert find_unguarded_phylo_imports(bare), "bare import must be flagged"
    assert find_unguarded_phylo_imports(unlabelled), "unlabelled skip must be flagged"
    assert find_unguarded_phylo_imports(guarded) == []
    assert find_unguarded_phylo_imports(deferred) == []


def test_no_test_module_imports_ete3_at_module_scope_unguarded():
    offenders = {}
    for path in sorted(TESTS_DIR.rglob("*.py")):
        problems = find_unguarded_phylo_imports(path.read_text(encoding="utf-8"))
        if problems:
            offenders[str(path.relative_to(REPO_ROOT))] = problems
    assert not offenders, (
        "Test modules must guard module-level ete3 imports with a labelled "
        "pytest.importorskip('ete3', reason=... NOT EXECUTED ...), otherwise "
        "'pytest tests' aborts during collection on Python >= 3.13: "
        f"{offenders}"
    )


# ---------------------------------------------------------------------------
# Gate 2 — an ete3-driven skip must say so in words the CI summary can grep
# ---------------------------------------------------------------------------

# A reason built from a name or attribute (reason=_NOT_EXECUTED) is not
# Resolvable statically, so it is not flagged; only literal strings are checked.
ETE3_UNAVAILABLE = re.compile(
    r"ete3 (?:is )?(?:not (?:importable|available)|unavailable)|without ete3",
    re.IGNORECASE,
)


def _skip_call_name(call: ast.Call) -> str:
    func = call.func
    if isinstance(func, ast.Attribute):
        return func.attr
    if isinstance(func, ast.Name):
        return func.id
    return ""


def find_unlabelled_ete3_skips(source: str) -> List[str]:
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:  # Pragma: no cover - a broken file is its own bug
        return [f"unparseable: {exc}"]
    problems: List[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or _skip_call_name(node) not in {
            "skip", "importorskip", "xfail",
        }:
            continue
        joined = " ".join(
            c.value
            for c in ast.walk(node)
            if isinstance(c, ast.Constant) and isinstance(c.value, str)
        )
        if ETE3_UNAVAILABLE.search(joined) and REQUIRE_TOKEN not in joined:
            problems.append(
                f"line {node.lineno}: skip caused by a missing ete3 must carry "
                f"'{REQUIRE_TOKEN}' so it can never be read as a pass"
            )
    return problems


def test_control_the_skip_labelling_scanner():
    unlabelled = (
        "import pytest\n\n\ndef test_x():\n"
        '    pytest.skip("ete3 unavailable: cannot measure")\n'
    )
    labelled = (
        "import pytest\n\n\ndef test_x():\n"
        '    pytest.skip("ete3 unavailable — NOT EXECUTED in this interpreter")\n'
    )
    unrelated = (
        "import pytest\n\n\ndef test_x():\n"
        '    pytest.skip("mafft missing on PATH")\n'
    )
    indirect = (
        "import pytest\n\n_NOT_EXECUTED = 'x'\n\n\ndef test_x():\n"
        "    pytest.importorskip('ete3', reason=_NOT_EXECUTED)\n"
    )
    assert find_unlabelled_ete3_skips(unlabelled), "unlabelled must be flagged"
    assert find_unlabelled_ete3_skips(labelled) == []
    assert find_unlabelled_ete3_skips(unrelated) == []
    assert find_unlabelled_ete3_skips(indirect) == []  # Documented limitation


def test_every_ete3_skip_is_labelled_not_executed():
    offenders = {}
    for path in sorted(TESTS_DIR.rglob("*.py")):
        problems = find_unlabelled_ete3_skips(path.read_text(encoding="utf-8"))
        if problems:
            offenders[str(path.relative_to(REPO_ROOT))] = problems
    assert not offenders, (
        "A test that could not run must be visibly NOT EXECUTED. "
        f"Fix these skip reasons: {offenders}"
    )


# ---------------------------------------------------------------------------
# Gate 3 — end-to-end collectability in *this* interpreter
# ---------------------------------------------------------------------------

def _collection_problems(returncode: int, out: str) -> List[str]:
    """Reasons this ``--collect-only`` run is untrustworthy; [] means healthy.

    Markers are deliberately precise. Matching the bare word "error" over
    collect-only output is useless, because test class names
    (``TestErrorScenario`` …) contain it — an earlier version of this gate failed for
    exactly that reason. Collection failures show up as ``ERROR`` lines and an
    "errors during collection" interruption; a successful run prints
    "N tests collected".
    """
    problems: List[str] = []
    if returncode != 0:
        problems.append(f"exit code {returncode}")
    if "errors during collection" in out:
        problems.append("pytest aborted during collection")
    if re.search(r"(?m)^ERROR\b", out):
        problems.append("ERROR lines present")
    # Unanchored on purpose: pytest pads the summary line with '=' signs, so a
    # "^" here would never match and the "no summary" branch would fire forever.
    collected = re.search(r"(\d+) tests? collected\b", out)
    if collected is None:
        problems.append("no 'N tests collected' summary line")
    elif int(collected.group(1)) <= 0:
        problems.append("zero tests collected")
    return problems


# Captured from `python3 -m pytest tests` (Python 3.14, ete3 unimportable)
# while the two bare module-level imports were still present.
# Kept here so the gate's alarm can be proven to fire, not just proven silent.
CAPTURED_COLLECTION_ABORT = """\
==================================== ERRORS ====================================
_______________ ERROR collecting tests/unit/test_cr_mad_root.py ________________
E   ModuleNotFoundError: No module named 'cgi'
_______________ ERROR collecting tests/unit/test_cr_taxonomy.py ________________
E   ModuleNotFoundError: No module named 'cgi'
ERROR tests/unit/test_cr_mad_root.py
ERROR tests/unit/test_cr_taxonomy.py
!!!!!!!!!!!!!!!!!!! Interrupted: 2 errors during collection !!!!!!!!!!!!!!!!!!!!
============================== 2 errors in 2.24s ===============================
"""


def test_control_the_alarm_fires_on_the_real_failure_transcript():
    problems = _collection_problems(2, CAPTURED_COLLECTION_ABORT)
    assert problems, "captured collection abort must be reported as unhealthy"
    # The healthy sample is the live run below; this control only proves the
    # Detector is not a no-op.
    assert _collection_problems(
        0, "========================= 721 tests collected in 1.95s ==========\n"
    ) == []


def test_collection_completes_in_this_interpreter():
    """``pytest --collect-only`` must exit 0 whether or not ete3 is importable.

    Collection never executes test bodies, so this cannot recurse into itself.
    """
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", "tests"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=600,
    )
    out = proc.stdout + proc.stderr
    tail = "\n".join(out.splitlines()[-12:])
    problems = _collection_problems(proc.returncode, out)
    assert not problems, (
        "The suite must be collectable in every supported interpreter; "
        f"{problems}:\n{tail}"
    )
