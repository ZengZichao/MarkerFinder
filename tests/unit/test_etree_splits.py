"""Pure-Python split-set RF / quartet / clade measurement tests.

Differential coverage for:mod:`markerfinder.utils.etree` — the ete3-free
fallback measurement path. The ete3 differential test is
skipped with an explicit "NOT EXECUTED" marker when ete3 is unavailable
(a skip must be visible, never silently counted as a pass).
"""

import math

import pytest

from markerfinder.utils import etree


class TestSplits:
    def test_splits_count_matches_resolution(self):
        # A fully resolved unrooted binary tree has exactly n-3 splits.
        for n in (4, 5, 7, 10):
            # Caterpillar: ((((t1,t2),t3)...,tn)
            inner = "(t1,t2)"
            for i in range(3, n + 1):
                inner = f"({inner},t{i})"
            newick = f"{inner};"
            splits, tips = etree.splits_of(newick)
            assert len(tips) == n
            assert len(splits) == n - 3

    def test_splits_symmetric_bipartition_dedup(self):
        # Root children of ((A,B),(C,D)) encode the SAME bipartition twice;
        # Normalization must dedupe to the single non-trivial split.
        splits, tips = etree.splits_of("((A,B),(C,D));")
        assert len(splits) == 1
        assert splits == {frozenset({"A", "B"})}


class TestRfDistance:
    def test_rf_symmetric_and_zero_for_identical(self):
        a = "(((A,B),C),((D,E),(F,G)));"
        raw, norm = etree.rf_distance(a, a)
        assert raw == 0
        assert norm == 0.0
        # Left-right mirrored rewriting of the same topology (recursive swap).
        mirrored = "(((F,G),(E,D)),(C,(B,A)));"
        raw2, norm2 = etree.rf_distance(a, mirrored)
        assert raw2 == 0
        assert norm2 == 0.0

    def test_rf_symmetric(self):
        t1, t2 = "((A,B),(C,D));", "((A,C),(B,D));"
        r12 = etree.rf_distance(t1, t2)
        r21 = etree.rf_distance(t2, t1)
        assert r12 == r21

    def test_rf_maximum_conflict_on_quartet(self):
        raw, norm = etree.rf_distance("((A,B),(C,D));", "((A,C),(B,D));")
        assert raw == 2
        assert norm == 1.0

    def test_rf_partial_shared_tips(self):
        # Only the shared tip set is compared.
        t1 = "(((A,B),C),D);"
        t2 = "(((A,B),X),Y);"
        raw, norm = etree.rf_distance(t1, t2)
        assert raw == 0
        assert norm == 0.0

    def test_rf_tip_cap_beyond_64(self):
        # >64 shared tips => NOT_MEASURABLE (TipCapExceeded).
        big = "(" + ",".join(f"T{i}" for i in range(65)) + ");"
        with pytest.raises(etree.TipCapExceeded):
            etree.rf_distance(big, big)

    def test_rf_at_tip_cap_boundary_ok(self):
        tips64 = [f"T{i}" for i in range(64)]
        t1 = "(" + ",".join(tips64) + ");"
        # Star tree: no splits, comparison trivially valid at the cap
        raw, norm = etree.rf_distance(t1, t1)
        assert raw == 0

    def test_rf_unparseable_raises(self):
        with pytest.raises(etree.NewickParseError):
            etree.rf_distance("((A,B;", "((A,B),(C,D));")


class TestQuartetTopology:
    def test_quartet_topology_three_resolutions(self):
        base = ("A", "B", "C", "D")
        q1 = etree.quartet_topology("((A,B),(C,D));", base)
        q2 = etree.quartet_topology("((A,C),(B,D));", base)
        q3 = etree.quartet_topology("((A,D),(B,C));", base)
        assert q1 == "((A,B),(C,D))"
        assert q2 == "((A,C),(B,D))"
        assert q3 == "((A,D),(B,C))"
        assert len({q1, q2, q3}) == 3

    def test_quartet_unresolved_star_form(self):
        q = etree.quartet_topology("(A,B,C,D);", ("A", "B", "C", "D"))
        assert q == "(A,B,C,D)"

    def test_quartet_missing_tip_returns_none(self):
        assert etree.quartet_topology("((A,B),(C,D));", ("A", "B", "C", "Z")) is None

    def test_quartet_rotation_invariant(self):
        base = ("A", "B", "C", "D")
        # Same tree written with different rotations / tip orders.
        variants = [
            "((A,B),(C,D));",
            "((B,A),(D,C));",
            "(C,D,(A,B));",
            "((D,C),(B,A));",
        ]
        tops = {etree.quartet_topology(v, base) for v in variants}
        assert len(tops) == 1
        assert tops == {"((A,B),(C,D))"}

    def test_quartet_induced_from_deeper_tree(self):
        # Induced topology on a 4-tip subset of a 6-tip tree.
        q = etree.quartet_topology("(((A,B),C),((D,E),F));", ("A", "B", "C", "D"))
        assert q == "((A,B),(C,D))"


class TestMonophylyFallback:
    def test_clade_monophyletic_true(self):
        nwk = "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));"
        assert etree.is_clade_monophyletic(nwk, {"B1", "B2"}) is True

    def test_clade_monophyletic_false(self):
        nwk = "((B1:0.1,A1:0.2),(B2:0.3,A2:0.4));"
        assert etree.is_clade_monophyletic(nwk, {"B1", "B2"}) is False

    def test_numeric_tip_labels_preserved(self):
        # Numeric labels must survive as tips, not be eaten as support values.
        tips = etree.tip_set("((123,12.5),(x1,y1));")
        assert tips == {"123", "12.5", "x1", "y1"}

    def test_internal_support_labels_not_tips(self):
        tips = etree.tip_set("((A,B)95:0.5,(C,D)0.87:0.4);")
        assert tips == {"A", "B", "C", "D"}

    def test_nhx_annotations_stripped(self):
        tips = etree.tip_set("((A[&x=1],B),(C,D)[q1=0.9]);")
        assert tips == {"A", "B", "C", "D"}


class TestTaxonomyFallbackParity:
    """Taxonomy.is_monophyletic / monophyly_proportion must work without ete3."""

    def test_is_monophyletic_without_ete3(self):
        from markerfinder.taxonomy import is_monophyletic

        tax_map = {
            "B1": {"domain": "Bacteria"}, "B2": {"domain": "Bacteria"},
            "A1": {"domain": "Archaea"}, "A2": {"domain": "Archaea"},
        }
        assert is_monophyletic(
            "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));", "Bacteria", tax_map, "domain"
        ) is True
        assert is_monophyletic(
            "((B1:0.1,A1:0.2),(B2:0.3,A2:0.4));", "Bacteria", tax_map, "domain"
        ) is False

    def test_monophyly_proportion_without_ete3(self):
        from markerfinder.taxonomy import monophyly_proportion

        tax_map = {
            "B1": {"genus": "Bac"}, "B2": {"genus": "Bac"},
            "A1": {"genus": "Arc"}, "A2": {"genus": "Arc"},
            "C1": {"genus": "Cen"}, "C2": {"genus": "Cen"},
        }
        nwk = "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));"
        # Cen tips absent from tree entirely -> only 2 testable taxa.
        prop, n_total, n_mono = monophyly_proportion(nwk, tax_map, "genus")
        assert n_total == 2
        assert n_mono == 2
        assert prop == 1.0

    def test_monophyly_proportion_no_reps_returns_none(self):
        from markerfinder.taxonomy import monophyly_proportion

        tax_map = {"B1": {"genus": "Bac"}, "A1": {"genus": "Arc"}}
        nwk = "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));"
        prop, n_total, n_mono = monophyly_proportion(nwk, tax_map, "genus")
        assert prop is None
        assert n_total == 0 and n_mono == 0

    def test_monophyly_error_when_taxon_absent(self):
        from markerfinder.exceptions import MonophylyError
        from markerfinder.taxonomy import is_monophyletic

        with pytest.raises(MonophylyError):
            is_monophyletic("((A,B),(C,D));", "Zeta", {"A": {"genus": "g"}}, "genus")


class TestEte3Differential:
    """差分测试：纯 Python split-set RF 与 ete3 逐树对比。

    ete3 不可用时 skip——skip 原因含 "NOT EXECUTED"，在 CI 汇总中标记为
    未执行，绝不当作通过。
    """

    _NOT_EXECUTED = (
        "ete3 unavailable in this interpreter — DIFFERENTIAL TEST NOT "
        "EXECUTED (every supported interpreter, Python >=3.10, can import "
        "ete3 once `markerfinder` is imported, so this means a broken "
        "install rather than an unsupported interpreter)"
    )

    @staticmethod
    def _random_resolved_newick(tips, rng):
        # Leaves must be plain strings: etree_fmt unwraps a node as an
        # ``(left, right)`` pair, so a 1-tuple leaf raised
        # "not enough values to unpack" and the differential test died before
        # Comparing anything (it could never pass on any interpreter).
        clades = [t for t in tips]
        while len(clades) > 2:
            i, j = sorted(rng.sample(range(len(clades)), 2), reverse=True)
            a, b = clades.pop(i), clades.pop(j)
            clades.append((a, b))
        (a, b) = clades
        return f"({etree_fmt(a)},{etree_fmt(b)});"

    def test_rf_agrees_with_ete3_when_available(self):
        pytest.importorskip("ete3", reason=self._NOT_EXECUTED)
        from ete3 import Tree as EteTree

        import random

        rng = random.Random(20260920)
        for case in range(30):
            n = rng.randint(4, 12)
            tips = [f"T{k:02d}" for k in range(n)]
            n1 = self._random_resolved_newick(tips, rng)
            n2 = self._random_resolved_newick(tips, rng)
            expected = EteTree(n1).robinson_foulds(EteTree(n2), unrooted_trees=True)
            got_rf, got_norm = etree.rf_distance(n1, n2)
            assert got_rf == expected[0], f"case {case}: {n1} vs {n2}"
            exp_norm = expected[0] / expected[1] if expected[1] > 0 else 0.0
            assert got_norm == pytest.approx(exp_norm), f"case {case}: {n1} vs {n2}"

    def test_quartet_agrees_with_ete3_when_available(self):
        pytest.importorskip("ete3", reason=self._NOT_EXECUTED)
        from ete3 import Tree as EteTree

        import random

        rng = random.Random(20260921)
        tips = [f"T{k:02d}" for k in range(8)]
        for case in range(20):
            nwk = self._random_resolved_newick(tips, rng)
            quartet = tuple(rng.sample(tips, 4))
            t = EteTree(nwk)
            t.prune(list(quartet), preserve_branch_length=True)
            expected = t.write(format=9)
            got = etree.quartet_topology(nwk, quartet)
            # Compare the informative 2|2 split, not the string: root position
            # And side ordering legitimately differ between the two emitters.
            exp_pairing = quartet_pairing(expected)
            got_pairing = quartet_pairing(got)

            assert got_pairing == exp_pairing, (
                f"case {case}: {nwk} {quartet}: {got} != {expected}"
            )


def etree_fmt(node) -> str:
    if isinstance(node, str):
        return node
    a, b = node
    return f"({etree_fmt(a)},{etree_fmt(b)})"


def quartet_pairing(newick: str):
    """Canonicalise a 4-taxon Newick string to its informative 2|2 split.

    A quartet's topology IS its single non-trivial bipartition, so two strings
    agree iff that split agrees — regardless of where the root was placed or
    how the sides are ordered. This must be computed by parsing the topology,
    NOT by string surgery: the previous ``strip(";")``/``split("),(")``
    version could not represent an ete3 rooted emission such as
    ``(T06,((T02,T01),T00))`` (whose informative split is T01,T02 | T00,T06)
    and therefore reported a mismatch between two identical topologies.
    Independent of ``etree.quartet_topology`` on purpose, so the comparison is
    not circular.
    """
    stack: list = [[]]
    token = ""
    for ch in newick:
        if ch in "():,":
            if token.strip():
                stack[-1].append(token.strip())
            token = ""
            if ch == "(":
                stack.append([])
            elif ch == ")":
                if len(stack) < 2:
                    raise AssertionError(f"unbalanced Newick: {newick!r}")
                child = stack.pop()
                stack[-1].append(child)
        else:
            token += ch
    if token.strip():
        stack[-1].append(token.strip())
    root = stack[0][0] if stack[0] else []

    all_tips: set = set()
    subsets: list = []

    def record(node):
        if isinstance(node, str):
            all_tips.add(node)
            return {node}
        s: set = set()
        for c in node:
            s |= record(c)
        subsets.append(s)
        return s

    record(root)
    if len(all_tips) != 4:
        raise AssertionError(f"expected 4 tips, got {sorted(all_tips)} in {newick!r}")
    informative = {
        frozenset([frozenset(s), frozenset(all_tips - s)])
        for s in subsets
        if len(s) == 2 and (all_tips - s)
    }
    if len(informative) != 1:
        raise AssertionError(
            f"no unique 2|2 split in {newick!r}: {informative}"
        )
    return next(iter(informative))
