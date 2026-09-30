"""Unit tests for tree utility functions."""

from pathlib import Path

import pytest

from markerfinder.utils.tree_utils import (
    calculate_rf_distance,
    get_quartet_topology,
    merge_newick_files,
    parse_newick_tips,
    strip_nhx_annotations,
)


class TestStripNhxAnnotations:
    def test_strip_simple_nhx(self):
        raw = "((A:0.1,B:0.2)[q1=0.9],C:0.3);"
        assert "[" not in strip_nhx_annotations(raw)

    def test_empty_string(self):
        assert strip_nhx_annotations("") == ""

    def test_strip_quoted_block(self):
        raw = "((A:0.1,B:0.2)'[q1=0.9;pp=1.0]',C:0.3);"
        cleaned = strip_nhx_annotations(raw)
        assert "[" not in cleaned
        assert "'" not in cleaned


class TestParseNewickTips:
    def test_simple_tree(self):
        assert parse_newick_tips("((A,B),(C,D));") == ["A", "B", "C", "D"]

    def test_deduplicates(self):
        assert parse_newick_tips("((A,A),(B,C));") == ["A", "B", "C"]


class TestCalculateRfDistance:
    def test_identical_trees(self):
        raw, norm = calculate_rf_distance("((A,B),(C,D));", "((A,B),(C,D));")
        assert raw == 0
        assert norm == 0.0

    def test_different_trees(self):
        raw, norm = calculate_rf_distance("((A,B),(C,D));", "((A,C),(B,D));")
        assert raw > 0
        assert 0 < norm <= 1.0

    def test_empty_tree(self):
        raw, norm = calculate_rf_distance("", "((A,B),(C,D));")
        assert raw is None
        assert norm is None


class TestGetQuartetTopology:
    def test_quartet_subset(self):
        topo = get_quartet_topology("((A,B),(C,D));", ("A", "B", "C", "D"))
        assert topo is not None

    def test_missing_tip_returns_none(self):
        topo = get_quartet_topology("((A,B),(C,D));", ("A", "B", "C", "X"))
        assert topo is None


class TestMergeNewickFiles:
    def test_merge(self, tmp_path):
        f1 = tmp_path / "a.nwk"
        f2 = tmp_path / "b.nwk"
        f1.write_text("(A,B);\n")
        f2.write_text("(C,D);\n")
        out = tmp_path / "merged.nwk"
        merge_newick_files([str(f1), str(f2)], str(out))
        content = out.read_text()
        assert "(A,B);" in content
        assert "(C,D);" in content
