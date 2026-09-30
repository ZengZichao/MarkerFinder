"""Code-review regression tests for ``utils/tree_utils.py`` (C29).

Locks in the ``calculate_rf_distance`` semantics change:
  - Failure (empty input, missing ete3, incompatible/unparseable trees) now
    returns ``(None, None)`` — a distinct "uncomputable / not applicable"
    sentinel — instead of the old ``(-1, -1.0)``.
  - A genuine identical-tree result is still ``(0, 0.0)``, which must NOT be
    confused with the ``None`` failure sentinel.
"""

import pytest

from markerfinder.utils.tree_utils import calculate_rf_distance


class TestRfDistanceSemantics:
    def test_identical_returns_zero_not_none(self):
        raw, norm = calculate_rf_distance("((A,B),(C,D));", "((A,B),(C,D));")
        assert raw == 0
        assert norm == 0.0
        assert raw is not None  # Distinct from the failure sentinel

    def test_empty_tree_returns_none(self):
        raw, norm = calculate_rf_distance("", "((A,B),(C,D));")
        assert raw is None and norm is None

    def test_both_empty_returns_none(self):
        raw, norm = calculate_rf_distance("", "")
        assert raw is None and norm is None

    def test_incompatible_tree_returns_none(self):
        raw, norm = calculate_rf_distance("not a tree at all", "((A,B),(C,D));")
        assert raw is None and norm is None

    def test_distinct_trees_return_real_distance(self):
        raw, norm = calculate_rf_distance("((A,B),(C,D));", "((A,C),(B,D));")
        assert raw is not None and raw > 0
        assert norm is not None and 0 < norm <= 1.0
