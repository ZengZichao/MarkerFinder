"""Additional unit tests for taxonomy parsing and monophyly utilities."""

from pathlib import Path

import pytest

from markerfinder.exceptions import MonophylyError, TaxonomyConflictError
from markerfinder.taxonomy import (
    _check_taxonomy_cycle,
    _detect_table_separator,
    _empty_taxonomy,
    detect_scope_rank,
    find_special_identifier,
    get_merged_level_maps,
    is_monophyletic,
    load_taxonomy_table,
    monophyly_proportion,
    parse_custom_levels,
    validate_table_tree_consistency,
)


class TestCustomLevels:
    def test_parse_custom_levels(self):
        cmap, names = parse_custom_levels("kingdom:_k_,subspecies:_ss_")
        assert cmap == {"k": "kingdom", "ss": "subspecies"}
        assert names == ["kingdom", "subspecies"]

    def test_get_merged_level_maps(self):
        a, b, names = get_merged_level_maps("kingdom:_k_")
        assert "k" in a
        assert "kingdom" in names


class TestTableSeparator:
    def test_detect_tab(self, tmp_path):
        path = tmp_path / "taxa.tsv"
        path.write_text("A\td__Bacteria\n")
        assert _detect_table_separator(str(path)) == "\t"

    def test_detect_comma(self, tmp_path):
        path = tmp_path / "taxa.csv"
        path.write_text("A,d__Bacteria\n")
        assert _detect_table_separator(str(path)) == ","


class TestLoadTaxonomyTable:
    def test_load_table(self, tmp_path):
        path = tmp_path / "taxa.tsv"
        path.write_text("A\td__Bacteria;p__Firmicutes\nB\td__Bacteria;p__Proteobacteria\n")
        table = load_taxonomy_table(str(path))
        assert table["A"]["domain"] == "Bacteria"
        assert table["B"]["phylum"] == "Proteobacteria"

    def test_skip_comments_and_empty(self, tmp_path):
        path = tmp_path / "taxa.tsv"
        path.write_text("# comment\n\nA\td__Bacteria\n")
        table = load_taxonomy_table(str(path))
        assert "A" in table

    def test_ignore_malformed(self, tmp_path):
        path = tmp_path / "taxa.tsv"
        path.write_text("A\td__Bacteria\nB\n")
        table = load_taxonomy_table(str(path), ignore_malformed=True)
        assert "A" in table
        assert "B" not in table

    def test_malformed_raises(self, tmp_path):
        path = tmp_path / "taxa.tsv"
        # Use a tip label containing a Unicode bidirectional override character,
        # Which the safety check rejects when ignore_malformed=False.
        path.write_text("A\u202Eevil\td__Bacteria\n")
        with pytest.raises(TaxonomyConflictError):
            load_taxonomy_table(str(path), ignore_malformed=False)

    def test_empty_file_raises(self, tmp_path):
        path = tmp_path / "empty.tsv"
        path.write_text("")
        with pytest.raises(ValueError):
            load_taxonomy_table(str(path))

    def test_header_row_is_not_ingested_as_a_tip(self, tmp_path):
        """``genome_id\\ttaxonomy`` is the shape every spreadsheet export
        produces. It used to be parsed like any other row, which left a phantom
        tip named ``genome_id`` with an all-None taxonomy in the mapping — a tip
        that then joins every tip-set comparison downstream (scope detection for
        the monophyly rank, tree/table cross-validation)."""
        path = tmp_path / "taxa.tsv"
        path.write_text(
            "genome_id\ttaxonomy\n"
            "A\td__Bacteria;p__Firmicutes\n"
            "B\td__Bacteria;p__Proteobacteria\n"
        )
        table = load_taxonomy_table(str(path))
        assert set(table) == {"A", "B"}, table
        assert "genome_id" not in table

    def test_empty_taxonomy_cell_is_still_a_tip(self, tmp_path):
        """A real tip with a MISSING taxonomy must not be eaten by header
        sniffing: it stays in the mapping with an all-None taxonomy."""
        path = tmp_path / "taxa.tsv"
        path.write_text("A\t\nB\td__Bacteria\n")
        table = load_taxonomy_table(str(path))
        assert set(table) == {"A", "B"}, table
        assert all(v is None for v in table["A"].values())

    def test_embedded_format_header_is_recognised(self, tmp_path):
        """--taxonomy-format embedded: the column label has no ``_x_``
        delimiter either, so the same rule must hold."""
        path = tmp_path / "taxa.tsv"
        path.write_text(
            "genome_id\ttaxonomy_segment\n"
            "A\t_d_Bacteria_p_Firmicutes_c_Bacilli_g_Bacillus\n"
        )
        table = load_taxonomy_table(str(path), taxonomy_format="embedded")
        assert set(table) == {"A"}, table


class TestTaxonomyCycle:
    def test_cycle_detected(self):
        existing = {"B": {"domain": "A"}}
        parsed = {"domain": "B"}
        with pytest.raises(TaxonomyConflictError):
            _check_taxonomy_cycle("A", parsed, existing, ignore_malformed=False, line_num=1)

    def test_cycle_ignored(self):
        existing = {"B": {"domain": "A"}}
        parsed = {"domain": "B"}
        # Should not raise when ignore_malformed=True
        _check_taxonomy_cycle("A", parsed, existing, ignore_malformed=True, line_num=1)


class TestValidateTableTreeConsistency:
    def test_consistency(self):
        table = {"A": {}, "B": {}}
        only_table, only_tree = validate_table_tree_consistency(table, {"A", "C"})
        assert only_table == {"B"}
        assert only_tree == {"C"}


class TestMonophyly:
    def test_monophyletic(self):
        tree = "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));"
        tax = {
            "B1": {"domain": "Bacteria", "genus": "B"},
            "B2": {"domain": "Bacteria", "genus": "B"},
            "A1": {"domain": "Bacteria", "genus": "A"},
            "A2": {"domain": "Bacteria", "genus": "A"},
        }
        assert is_monophyletic(tree, "B", tax, "genus") is True

    def test_non_monophyletic(self):
        tree = "((B1:0.1,A1:0.2),(B2:0.3,A2:0.4));"
        tax = {
            "B1": {"domain": "Bacteria", "genus": "B"},
            "B2": {"domain": "Bacteria", "genus": "B"},
            "A1": {"domain": "Bacteria", "genus": "A"},
            "A2": {"domain": "Bacteria", "genus": "A"},
        }
        assert is_monophyletic(tree, "B", tax, "genus") is False

    def test_taxon_not_found_raises(self):
        tree = "((A,B),(C,D));"
        tax = {"A": {"domain": "Bacteria"}}
        with pytest.raises(MonophylyError):
            is_monophyletic(tree, "X", tax, "genus")


class TestMonophylyProportion:
    def test_perfect_proportion(self):
        tree = "((B1,B2),(A1,A2));"
        tax = {
            "B1": {"domain": "B", "genus": "B"},
            "B2": {"domain": "B", "genus": "B"},
            "A1": {"domain": "B", "genus": "A"},
            "A2": {"domain": "B", "genus": "A"},
        }
        prop, total, mono = monophyly_proportion(tree, tax, "genus")
        assert prop == 1.0
        assert total == 2

    def test_no_representatives(self):
        tree = "((A,B),(C,D));"
        tax = {t: {"domain": "B"} for t in ["A", "B", "C", "D"]}
        prop, total, mono = monophyly_proportion(tree, tax, "genus")
        assert prop is None
        assert total == 0


class TestDetectScopeRank:
    def test_scope_genus(self):
        tax = {
            "A": {"domain": "B", "phylum": "P", "genus": "G1"},
            "B": {"domain": "B", "phylum": "P", "genus": "G1"},
        }
        assert detect_scope_rank(tax, ["A", "B"]) == "genus"

    def test_scope_phylum(self):
        tax = {
            "A": {"domain": "B", "phylum": "P", "genus": "G1"},
            "B": {"domain": "B", "phylum": "P", "genus": "G2"},
        }
        assert detect_scope_rank(tax, ["A", "B"]) == "phylum"

    def test_no_scope(self):
        tax = {
            "A": {"domain": "B"},
            "B": {"domain": "A"},
        }
        assert detect_scope_rank(tax, ["A", "B"]) is None


class TestFindSpecialIdentifier:
    def test_luca(self):
        tree = "((A1,A2),(B1,B2));"
        tax = {
            "A1": {"domain": "Archaea"},
            "A2": {"domain": "Archaea"},
            "B1": {"domain": "Bacteria"},
            "B2": {"domain": "Bacteria"},
        }
        node = find_special_identifier(tree, tax, "LUCA")
        assert node is not None

    def test_laca(self):
        tree = "((A1,A2),(B1,B2));"
        tax = {
            "A1": {"domain": "Archaea"},
            "A2": {"domain": "Archaea"},
            "B1": {"domain": "Bacteria"},
            "B2": {"domain": "Bacteria"},
        }
        node = find_special_identifier(tree, tax, "LACA")
        leaves = {n.name for n in node.get_leaves()}
        assert leaves == {"A1", "A2"}

    def test_unknown_identifier_raises(self):
        tree = "(A,B);"
        tax = {"A": {"domain": "Bacteria"}, "B": {"domain": "Bacteria"}}
        with pytest.raises(MonophylyError):
            find_special_identifier(tree, tax, "UNKNOWN")
