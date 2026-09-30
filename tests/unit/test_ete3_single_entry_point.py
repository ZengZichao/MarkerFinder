"""One sanctioned ete3 entry point — and a ratchet on the remaining imports.

``etree.require_ete3`` is the single sanctioned way to reach ete3. A bare claim
of "one entry point" is worth nothing on its own: the entry point can exist
while direct ``from ete3 import ...`` statements are still scattered across
production modules. This suite defends against exactly that gap.

The two gateway measurements named in the design notes
(``calculate_rf_distance`` and ``get_quartet_topology`` in ``tree_utils.py``) now
go through ``require_ete3``, and so does every other site. The ledger below is
empty, and the ratchet stays (an empty ledger that any new scattered import makes
non-empty), because a gate that has passed once is not a gate that keeps passing.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

PACKAGE = Path("markerfinder")
IMPORT_LINE = re.compile(r"^\s*(?:from ete3 import|import ete3\b)")

# The sanctioned entry point itself.
EXEMPT = {"utils/etree.py"}

# Debt is paid in full: no module imports ete3 anywhere except utils/etree.py.
# The dict must stay empty — an entry here is a promise that someone will pay,
# and such promises rot if left unguarded (this suite started from an entry
# point that claimed to be the only one but had zero callers).
# Migrated: tree_utils, hgt_filter.detect_monophyly, taxonomy's three
# availability probes, models/tree.py, mad_root, self_test,
# phylogenetic_inference, dependency_check.
DEBT_LEDGER: dict[str, int] = {}


def _direct_imports_by_file() -> dict:
    counts: dict = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        rel = str(path.relative_to(PACKAGE))
        if rel in EXEMPT:
            continue
        n = sum(
            1 for line in path.read_text(encoding="utf-8", errors="replace").split("\n")
            if IMPORT_LINE.match(line)
        )
        if n:
            counts[rel] = n
    return counts


def test_the_two_named_gateways_now_use_the_entry_point():
    """'s cited evidence files must no longer import ete3 directly."""
    import inspect

    from markerfinder.utils import tree_utils

    source = inspect.getsource(tree_utils)
    assert "from ete3 import" not in source, (
        "tree_utils is the historical defect site; it must reach "
        "ete3 only through etree.require_ete3()"
    )
    assert source.count("require_ete3()") >= 2


def test_require_ete3_actually_has_production_callers():
    """An entry point nobody calls is not an entry point (baseline )."""
    callers = []
    for path in sorted(PACKAGE.rglob("*.py")):
        if path.name == "etree.py":
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if "require_ete3()" in text:
            callers.append(str(path.relative_to(PACKAGE)))
    assert callers, (
        "require_ete3() is defined but called nowhere: the single-failure-"
        "type guarantee would be fiction"
    )


def test_direct_import_ledger_matches_reality_exactly():
    """Frozen per-file counts: new scattered imports fail, and so does drift."""
    assert _direct_imports_by_file() == DEBT_LEDGER, (
        "direct ete3 imports changed. If you MIGRATED one, shrink DEBT_LEDGER; "
        "if you ADDED one, use markerfinder.utils.etree.require_ete3() instead. "
        f"actual={_direct_imports_by_file()}"
    )


def test_direct_imports_are_strictly_below_the_baseline():
    # The ledger is empty now; the ratchet keeps it below the historical count.
    total = sum(DEBT_LEDGER.values())
    assert total == 0, total
    assert total < 19


def test_fallback_semantics_survive_the_entry_point(monkeypatch):
    """Require_ete3 must raise the one declared failure type when ete3 is gone."""
    import builtins

    from markerfinder.exceptions import PhyloToolUnavailable
    from markerfinder.utils import tree_utils

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "ete3":
            raise ImportError("simulated: ete3 not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    rf, norm = tree_utils.calculate_rf_distance(
        "((A,B),(C,D));", "((A,C),(B,D));"
    )
    # Must fall back to the pure-Python measurement, NOT return (None, None)
    # And NOT a placeholder: 2 of 2 possible internal splits differ here.
    assert rf is not None and norm is not None
    assert abs(norm - 1.0) < 1e-9, (rf, norm)


def test_mad_root_refuses_to_echo_an_unrooted_tree(monkeypatch):
    """"cannot root" must be an exception, never a plausible tree.

    mad_root used to return its input unchanged (one warning line) whenever ete3
    was unusable, so a caller measuring clades could not tell that it had been
    handed an unrooted tree. Runs on every interpreter: ete3 is made unimportable
    by the patch, so this does not depend on the environment.
    """
    import markerfinder.utils.etree as etree_mod
    from markerfinder.exceptions import PhyloToolUnavailable
    from markerfinder.utils.mad_root import mad_root

    def _unavailable():
        raise PhyloToolUnavailable("simulated: ete3 not importable")

    monkeypatch.setattr(etree_mod, "require_ete3", _unavailable)
    with pytest.raises(PhyloToolUnavailable) as excinfo:
        mad_root("((A:1,B:1):1,(C:1,D:1):1);")
    assert "no pure-Python rooting fallback" in str(excinfo.value)
