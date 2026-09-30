"""Code-review regression tests for ``taxonomy.py`` (#4) monophyly_proportion.

Locks in the performance optimization where ``monophyly_proportion`` pre-builds
a ``node_for_clade`` map for O(1) lookup instead of re-walking the tree per
taxon with ``get_common_ancestor``. The *returned proportion must be unchanged*
(behavioral parity) — verified here against an independent reference that uses
the old per-taxon ``is_monophyletic`` (get_common_ancestor) approach. A
multi-taxon query is exercised as well.
"""

from collections import defaultdict

import pytest

# Differential test: the independent reference reimplements the old per-taxon
# Approach with ete3, so it needs ete3 importable. Guarded at module scope —
# A bare ``from ete3 import...`` made ``pytest tests`` abort during collection
# On Python >= 3.13 (ete3 imports the stdlib ``cgi`` module, removed in 3.13),
# Which reported zero results instead of a labelled NOT EXECUTED skip.
pytest.importorskip(
    "ete3",
    reason=(
        "ete3 unavailable in this interpreter — DIFFERENTIAL TEST NOT "
        "EXECUTED (must run in a 3.10-3.12 environment)"
    ),
)
from ete3 import Tree as EteTree  # Noqa: E402

from markerfinder.taxonomy import is_monophyletic, monophyly_proportion


def _ref_proportion(tree, tax, level):
    """Independent reference: old per-taxon ``is_monophyletic`` approach."""
    t = EteTree(tree, format=1)
    all_tips = {n.name for n in t.get_leaves() if n.name}
    taxon_to_tips = defaultdict(set)
    for tip, d in tax.items():
        v = d.get(level)
        if v:
            taxon_to_tips[v].add(tip)
    total = 0
    mono = 0
    for label, tips in taxon_to_tips.items():
        tips_on = tips & all_tips
        if len(tips_on) < 2:
            continue
        total += 1
        if is_monophyletic(tree, label, tax, level):
            mono += 1
    return (mono / total, total, mono) if total else (None, 0, 0)


class TestMonophylyProportionParity:
    def test_perfect_matches_reference(self):
        # A clean trifurcation: each genus is its own clade, so all three are
        # Monophyletic (proportion == 1.0).
        tree = "((B1,B2),(A1,A2),(C1,C2));"
        tax = {
            "B1": {"domain": "B", "genus": "B"},
            "B2": {"domain": "B", "genus": "B"},
            "A1": {"domain": "B", "genus": "A"},
            "A2": {"domain": "B", "genus": "A"},
            "C1": {"domain": "B", "genus": "C"},
            "C2": {"domain": "B", "genus": "C"},
        }
        res = monophyly_proportion(tree, tax, "genus")
        ref = _ref_proportion(tree, tax, "genus")
        assert res == ref
        assert res[0] == 1.0

    def test_broken_matches_reference(self):
        tree = "((B1,A1),(B2,A2));"
        tax = {
            "B1": {"domain": "B", "genus": "B"},
            "B2": {"domain": "B", "genus": "B"},
            "A1": {"domain": "B", "genus": "A"},
            "A2": {"domain": "B", "genus": "A"},
        }
        res = monophyly_proportion(tree, tax, "genus")
        ref = _ref_proportion(tree, tax, "genus")
        assert res == ref
        assert res[0] == 0.0

    def test_mixed_multi_taxon_query(self):
        # Three genera; only B is monophyletic -> 1/3.
        tree = "((B1,B2),(A1,C1,(A2,C2)));"
        tax = {
            "B1": {"domain": "B", "genus": "B"},
            "B2": {"domain": "B", "genus": "B"},
            "A1": {"domain": "B", "genus": "A"},
            "A2": {"domain": "B", "genus": "A"},
            "C1": {"domain": "B", "genus": "C"},
            "C2": {"domain": "B", "genus": "C"},
        }
        res = monophyly_proportion(tree, tax, "genus")
        ref = _ref_proportion(tree, tax, "genus")
        assert res == ref
        prop, total, mono = res
        assert total == 3 and mono == 1
        assert abs(prop - 1 / 3) < 1e-9
