"""Code-review regression tests for ``models/tree.py`` get_tips (#6 B17/B18).

Locks in the fix where ``get_tips`` (ete3-backed) correctly returns tips whose
names are purely numeric (e.g. ``((1,2),(3,4));``) instead of silently
dropping them as if they were bootstrap/support values. Also verifies the
new ``logger.warning`` for duplicate tip names is emitted.
"""

import logging

import pytest

from markerfinder.models.tree import Tree


class TestGetTipsNumeric:
    def test_purely_numeric_tips(self):
        tips = Tree(newick="((1,2),(3,4));").get_tips()
        assert set(tips) == {"1", "2", "3", "4"}
        assert len(tips) == 4

    def test_numeric_tips_not_dropped(self):
        tips = Tree(newick="((1,2),(3,4));").get_tips()
        # None of the four tips should be silently discarded.
        assert all(t in {"1", "2", "3", "4"} for t in tips)
        assert len(tips) == 4

    def test_duplicate_numeric_warns(self, caplog):
        with caplog.at_level(logging.WARNING):
            tips = Tree(newick="((1,1),(2,3));").get_tips()
        assert set(tips) == {"1", "2", "3"}
        assert "duplicate" in caplog.text.lower()
