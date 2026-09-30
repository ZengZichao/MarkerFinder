"""Unit tests for MAD rooting utility.

MAD rooting has no pure-Python equivalent ( covers only the monophyly
*test*), and ``mad_root`` now raises PhyloToolUnavailable instead of echoing its
input when ete3 is unusable. These four cases therefore only exist in an
environment where ete3 imports: without it they are reported NOT EXECUTED rather
than silently passing off the echo path, which is what they used to do here.
"""

import pytest

pytest.importorskip(
    "ete3",
    reason=(
        "ete3 unavailable in this interpreter — MAD ROOTING TESTS NOT "
        "EXECUTED (rooting has no pure-Python fallback; needs 3.10-3.12)"
    ),
)

from markerfinder.utils.mad_root import mad_root


class TestMadRoot:
    def test_small_tree_returns_rooted_newick(self):
        nwk = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
        rooted = mad_root(nwk)
        # Should return a valid Newick string with a root.
        assert rooted.startswith("(")
        assert rooted.strip().endswith(";")
        assert "A" in rooted and "B" in rooted and "C" in rooted and "D" in rooted

    def test_tree_without_branch_lengths(self):
        nwk = "((A,B),(C,D));"
        rooted = mad_root(nwk)
        assert "A" in rooted and "B" in rooted and "C" in rooted and "D" in rooted

    def test_too_few_tips_returns_unchanged(self):
        nwk = "(A:0.1,B:0.2);"
        rooted = mad_root(nwk)
        assert rooted == nwk

    def test_malformed_newick_returns_unchanged(self):
        nwk = "not a tree"
        rooted = mad_root(nwk)
        assert rooted == nwk
