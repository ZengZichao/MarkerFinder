"""Unit tests for taxonomy parsing and monophyly checks."""

import pytest

from markerfinder.taxonomy import (
    _empty_taxonomy,
    parse_taxonomy,
    parse_taxonomy_embedded,
    parse_taxonomy_table,
    merge_taxonomy,
    STANDARD_LEVELS,
)
from markerfinder.exceptions import MonophylyError


class TestFormatAParseReverse:
    def test_full_taxonomy(self):
        label = "GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_Cyanobacteriia_o_Cyanobacteriales_f_Prochloraceae_g_Prochloron"
        result = parse_taxonomy_embedded(label, mode="reverse")
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Cyanobacteriota"
        assert result["class"] == "Cyanobacteriia"
        assert result["order"] == "Cyanobacteriales"
        assert result["family"] == "Prochloraceae"
        assert result["genus"] == "Prochloron"
        assert result["species"] is None

    def test_missing_levels(self):
        label = "GB_GCA_000_d_Bacteria_g_Prochloron"
        result = parse_taxonomy_embedded(label, mode="reverse")
        assert result["domain"] == "Bacteria"
        assert result["genus"] == "Prochloron"
        assert result["phylum"] is None
        assert result["species"] is None

    def test_empty_label(self):
        result = parse_taxonomy_embedded("", mode="reverse")
        assert all(v is None for v in result.values())

    def test_no_taxonomy_delimiters(self):
        result = parse_taxonomy_embedded("plain_name_without_taxonomy", mode="reverse")
        assert all(v is None for v in result.values())


class TestFormatAParseGreedy:
    def test_greedy_mode(self):
        label = "GB_GCA_d_Bacteria_p_Cyanobacteriota_g_Prochloron"
        result = parse_taxonomy_embedded(label, mode="greedy")
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Cyanobacteriota"
        assert result["genus"] == "Prochloron"

    def test_greedy_no_d(self):
        label = "GB_GCA_p_Cyanobacteriota_g_Prochloron"
        result = parse_taxonomy_embedded(label, mode="greedy")
        assert result["domain"] is None


class TestFormatAParseSegment:
    def test_segment_mode(self):
        label = "GB_GCA_d_Bacteria_p_Cyanobacteriota_g_Prochloron"
        result = parse_taxonomy_embedded(label, mode="segment")
        assert result["domain"] == "Bacteria"


class TestFormatBParse:
    def test_full_taxonomy(self):
        s = "d__Archaea;p__Thermoproteota;c__Korarchaeia;o__Korarchaeales;f__Korarchaeaceae;g__WALU01;s__"
        result = parse_taxonomy_table(s)
        assert result["domain"] == "Archaea"
        assert result["phylum"] == "Thermoproteota"
        assert result["genus"] == "WALU01"
        assert result["species"] is None

    def test_custom_separator(self):
        s = "d__Bacteria|p__Cyanobacteriota|g__Prochloron"
        result = parse_taxonomy_table(s, sep="|")
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Cyanobacteriota"

    def test_empty_string(self):
        result = parse_taxonomy_table("")
        assert all(v is None for v in result.values())

    def test_invalid_pair(self):
        s = "d__Bacteria;invalid_pair;g__Prochloron"
        result = parse_taxonomy_table(s)
        assert result["domain"] == "Bacteria"
        assert result["genus"] == "Prochloron"

    def test_unknown_level_prefix(self):
        s = "d__Bacteria;z__Unknown"
        result = parse_taxonomy_table(s)
        assert result["domain"] == "Bacteria"


class TestAutoDetectFormat:
    def test_format_b_detection(self):
        result = parse_taxonomy("d__Bacteria;p__Cyanobacteriota")
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Cyanobacteriota"

    def test_format_a_detection(self):
        result = parse_taxonomy("GB_GCA_d_Bacteria_p_Cyanobacteriota")
        assert result["domain"] == "Bacteria"


class TestMergeTaxonomy:
    def test_table_priority(self):
        emb = _empty_taxonomy()
        emb["domain"] = "Bacteria"
        emb["phylum"] = "Firmicutes"

        tab = _empty_taxonomy()
        tab["domain"] = "Bacteria"
        tab["phylum"] = "Cyanobacteriota"

        result = merge_taxonomy(emb, tab, priority="table")
        assert result["phylum"] == "Cyanobacteriota"

    def test_embedded_priority(self):
        emb = _empty_taxonomy()
        emb["domain"] = "Bacteria"
        emb["phylum"] = "Firmicutes"

        tab = _empty_taxonomy()
        tab["domain"] = "Bacteria"
        tab["phylum"] = "Cyanobacteriota"

        result = merge_taxonomy(emb, tab, priority="embedded")
        assert result["phylum"] == "Firmicutes"

    def test_no_conflict(self):
        emb = _empty_taxonomy()
        emb["domain"] = "Bacteria"

        tab = _empty_taxonomy()
        tab["phylum"] = "Cyanobacteriota"

        result = merge_taxonomy(emb, tab)
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Cyanobacteriota"

    def test_both_none(self):
        result = merge_taxonomy(_empty_taxonomy(), _empty_taxonomy())
        assert all(v is None for v in result.values())


class TestEmptyTaxonomy:
    def test_all_none(self):
        tax = _empty_taxonomy()
        for level in STANDARD_LEVELS:
            assert tax[level] is None
