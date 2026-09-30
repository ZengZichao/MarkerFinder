"""``models/tree.py``: an unreadable tree must not answer the question.

Found while consolidating the ete3 entry point. Every
method in this dataclass carried the same shape::

    try:
        from ete3 import Tree as EteTree
...
    except Exception:
        return <a value that also means a real biological answer>

``get_induced_tree`` handed back the **full tree** as if it had been pruned;
``contains_bipartition`` and ``are_sister_groups`` answered ``False`` ("absent");
``internal_branches`` answered ``[]`` ("no internal branches");
``get_quartet_topology`` answered ``""``; and ``get_average_support`` answered
``0.0`` for a tree that simply carries no support labels ("every branch has zero
support"). Each of those is a claim, indistinguishable from a measurement —
exactly the G1/ placeholder class that forbids and that blames
for the retracted sponge paper's failure mode (Steenwyk & King 2025,
Science, doi:10.1126/science.adw9456; retracted 2026-02-05).

The placeholder was not inert. ``CoalescentInference._filter_gene_trees`` (the
``--min-gene-tree-support`` gate, a documented public CLI flag) compared the
returned ``None`` against a float and raised::

    TypeError: '<' not supported between instances of 'NoneType' and 'float'

so a run that hit one unparseable gene tree aborted in Phase 2. That is reproducible
on Python 3.10 *with* ete3 installed — ``n_tips`` uses the regex fallback and
counts 4 tips, while ``get_average_support`` uses ete3 and gives up. The two
parsers disagree inside one method call.

``mypy --strict`` would have caught the comparison, but ``models/tree.py`` and
``modules/phylogenetic_inference.py`` are not in 's nine-module list, so the
type hole sat outside the gate. Recorded in the closure report as a candidate for
widening that list.

Reading supports was worse than a type mismatch, though. Two further defects were
measured here on ete3 3.1.3:

* ``pipeline._gene_tree_mean_support`` walked ``tree.traverse`` while its only
  caller passes a ``models.tree.Tree`` (no such method), so every marker landed in
  the blanket ``except`` and the helper returned ``None`` for an entire run — the
  PIS / effective-site back-fill never measured anything. Its two unit tests
  passed a ``models.Tree`` too and asserted ``is None``, i.e. they locked in the
  inertness.
* Had it received a real ete3 tree, it read ``node.support`` — which ete3
  defaults to ``1.0`` for unlabelled nodes — and clamped with
  ``min(max(s, 0.0), 1.0)``, mapping IQ-TREE UFBOOT percentages onto 1.0 as well.
  Measured: UFBOOT 98/76, UFBOOT 5/3 and a tree with no labels at all all
  returned ``1.0``. The proxy meant to rank markers by informativeness could not
  separate a 3% branch from a 98% one, and reported maximum support for a tree
  carrying no support data (the zero-discrimination failure mode).
"""

from __future__ import annotations

import logging

import pytest

from markerfinder.config import PhylogeneticConfig
from markerfinder.models.tree import Tree
from markerfinder.modules.phylogenetic_inference import CoalescentInference

# Four tips for the regex counter, but syntactically broken for ete3.
BROKEN = "((A,B),(C,D)"
VALID = "((A,B),(C,D));"


def _ete3_available() -> bool:
    from markerfinder.utils.etree import require_ete3

    try:
        require_ete3()
        return True
    except Exception:
        return False


# ── the probe itself ──────────────────────────────────────────────────────

def test_control_broken_tree_is_visible_to_tips_but_not_to_support():
    """Proves the fixture reaches the buggy path, so the gates below are not
    testing an input that would have been skipped anyway."""
    tree = Tree(newick=BROKEN)
    assert tree.n_tips == 4, "regex fallback must still count the four tips"
    assert tree.get_average_support() is None


def test_control_the_unguarded_comparison_genuinely_crashes():
    """Positive control: the TypeError this package fixes is real."""
    with pytest.raises(TypeError):
        Tree(newick=BROKEN).get_average_support() < 0.9


# ── the crash fix ─────────────────────────────────────────────────────────

def _gate(min_support):
    config = PhylogeneticConfig()
    config.min_gene_tree_support = min_support
    return CoalescentInference(config)


def test_support_gate_drops_unmeasurable_tree_and_says_so(caplog):
    inference = _gate(0.9)
    with caplog.at_level(logging.WARNING, logger="markerfinder"):
        kept = inference._filter_gene_trees({"m_bad": Tree(newick=BROKEN)})
    assert kept == {}, "a tree whose support cannot be read must not pass a " \
                       "quality gate silently"
    text = caplog.text
    assert "NOT MEASURABLE" in text, text
    assert "m_bad" in text, text


def test_support_gate_still_keeps_a_good_tree(caplog):
    """Two-sided control: the new branch must not swallow healthy trees."""
    inference = _gate(0.9)
    good = Tree(newick=VALID, support_values={"n1": 100.0, "n2": 100.0})
    with caplog.at_level(logging.WARNING, logger="markerfinder"):
        kept = inference._filter_gene_trees({"m_good": good})
    assert list(kept) == ["m_good"]
    assert "NOT MEASURABLE" not in caplog.text


def test_support_gate_disabled_leaves_real_low_support_alone():
    """Min_support=None keeps the historical n_tips>=4-only behaviour."""
    inference = _gate(None)
    kept = inference._filter_gene_trees({"m_bad": Tree(newick=BROKEN)})
    assert list(kept) == ["m_bad"]


# ── the placeholder exits ─────────────────────────────────────────────────

def test_reserved_apis_answer_none_rather_than_a_claim():
    tree = Tree(newick=BROKEN)
    assert tree.get_induced_tree({"A", "B", "C", "D"}) is None, (
        "handing back the unpruned tree was the worst variant: it looked like a "
        "result and silently changed every downstream measurement"
    )
    assert tree.contains_bipartition({"A", "B"}, {"C", "D"}) is None
    assert tree.are_sister_groups({"A", "B"}, {"C", "D"}) is None
    assert tree.internal_branches is None
    assert tree.get_mrca_depth({"A", "B"}) is None
    assert tree.get_quartet_topology(("A", "B", "C", "D")) is None


def test_valid_tree_still_measures_when_ete3_is_present():
    if not _ete3_available():
        pytest.skip(
            "ete3 unavailable in this interpreter — ASSERTION NOT EXECUTED "
            "(the ete3 branch of these APIs needs a 3.10-3.12 environment)"
        )
    tree = Tree(newick=VALID)
    induced = tree.get_induced_tree({"A", "B", "C"})
    assert induced is not None and set(induced.get_tips()) == {"A", "B", "C"}
    assert tree.contains_bipartition({"A", "B"}, {"C", "D"}) is True
    assert tree.are_sister_groups({"A", "B"}, {"C", "D"}) is True
    assert tree.are_sister_groups({"A", "C"}, {"B", "D"}) is False
    branches = tree.internal_branches
    assert branches is not None and len(branches) == 2
    depth = tree.get_mrca_depth({"A", "B"})
    assert depth is not None and depth >= 1
    assert tree.get_quartet_topology(("A", "B", "C", "D"))

    # Two-sided control: a real absence must still be a False, never None.
    assert tree.contains_bipartition({"A", "C"}, {"B", "D"}) is False


def test_no_placeholder_number_for_a_tree_without_support_labels():
    """A tree with no internal labels has an UNDEFINED mean support, on either
    interpreter: with ete3 it parses and yields no supports, without ete3 the
    read itself is not measurable. Both must be None, never 0.0."""
    assert Tree(newick=VALID).get_average_support() is None
    # Two-sided control: when labels exist the mean really is a number.
    rated = Tree(newick=VALID, support_values={"n1": 80.0, "n2": 90.0})
    assert rated.get_average_support() == pytest.approx(85.0)


# ── the support-reading contract (PIS back-fill input) ───────────────────

UFBOOT_STRONG = "((A:1,B:1)98:1,(C:1,D:1)76:1);"
UFBOOT_WEAK = "((A:1,B:1)5:1,(C:1,D:1)3:1);"
FASTTREE = "((A:1,B:1)0.98:1,(C:1,D:1)0.76:1);"
NAMED_ONLY = "((A:1,B:1)Node1:1,(C:1,D:1)Node2:1);"


def _mean(newick: str):
    from markerfinder.pipeline import _gene_tree_mean_support

    # The pipeline passes a models.Tree, NOT an ete3 tree. Feeding anything else
    # Here would re-create the bug this package fixes.
    return _gene_tree_mean_support(Tree(newick=newick))


def test_control_unlabelled_tree_contributes_nothing():
    """Ete3's TreeNode.support defaults to 1.0; the label text is the evidence."""
    assert Tree(newick=VALID).support_labels() == []
    assert _mean(VALID) is None


def test_percent_scale_is_not_collapsed_by_a_zero_one_clamp():
    """Before the fix all three of these returned the identical 1.0."""
    if not _ete3_available():
        pytest.skip(
            "ete3 unavailable in this interpreter — ASSERTION NOT EXECUTED "
            "(needs a 3.10-3.12 environment to read labels)"
        )
    assert Tree(newick=UFBOOT_STRONG).support_labels() == [98.0, 76.0]
    assert _mean(UFBOOT_STRONG) == pytest.approx(0.87)
    assert _mean(UFBOOT_WEAK) == pytest.approx(0.04)
    assert _mean(FASTTREE) == pytest.approx(0.87)
    # The discrimination property itself: near-zero support must not tie with
    # Near-perfect support. This is the whole point of the PIS proxy.
    assert _mean(UFBOOT_WEAK) < _mean(UFBOOT_STRONG)
    # A clade name is not a support value.
    assert _mean(NAMED_ONLY) is None


def test_zero_support_counts_and_absurd_scales_refuse():
    if not _ete3_available():
        pytest.skip(
            "ete3 unavailable in this interpreter — ASSERTION NOT EXECUTED "
            "(needs a 3.10-3.12 environment to read labels)"
        )
    # One branch at 0%, one at 100%: the mean is 0.5. Dropping the zero (the old
    # "s > 0" filter) reported 1.0 instead.
    assert _mean("((A:1,B:1)0:1,(C:1,D:1)100:1);") == pytest.approx(0.5)
    # A single label above 1 proves the tree is percent-scaled, so 0.9 is 0.9%.
    assert _mean("((A:1,B:1)0.9:1,(C:1,D:1)95:1);") == pytest.approx(0.4795)
    # Nothing legitimate is above 100.
    assert _mean("((A:1,B:1)980:1,(C:1,D:1)760:1);") is None


def test_internal_branches_report_no_support_rather_than_a_default():
    if not _ete3_available():
        pytest.skip(
            "ete3 unavailable in this interpreter — ASSERTION NOT EXECUTED "
            "(needs a 3.10-3.12 environment to read labels)"
        )
    branches = Tree(newick=VALID).internal_branches
    assert branches is not None and len(branches) == 2
    assert {b["support"] for b in branches} == {None}, (
        "an unlabelled branch must report no support, not ete3's 1.0 default"
    )
    labelled = Tree(newick=UFBOOT_STRONG).internal_branches
    assert sorted(b["support"] for b in labelled) == [76.0, 98.0]
