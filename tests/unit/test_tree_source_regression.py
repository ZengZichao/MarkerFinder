"""Behavior-equivalence tests for the consolidated shared utilities.

These pin the exact behaviour of ``markerfinder.utils.tree_source`` and the
shared ``parse_fasta`` so that the duplication-consolidation refactor cannot
silently change priority decisions or FASTA parsing.
"""

import pytest

from markerfinder.utils.tree_source import (
    prioritize_tree_source,
    coalescent_is_degenerate,
    DEFAULT_CONCAT_LABEL,
)
from markerfinder.utils.io import parse_fasta


class TestPrioritizeTreeSource:
    def test_astral_wins_regardless_of_concat(self):
        assert prioritize_tree_source("astral", True) == "astral"
        assert prioritize_tree_source("astral", False) == "astral"

    def test_consensus_wins_regardless_of_concat(self):
        assert prioritize_tree_source("consensus", True) == "consensus"
        assert prioritize_tree_source("consensus", False) == "consensus"

    def test_concat_source_key_passes_through(self):
        # "concat" is itself a proper merged-tree tier (rank >= concat).
        assert prioritize_tree_source("concat", True) == "concat"

    def test_no_coalescent_falls_back_to_concat_label(self):
        assert prioritize_tree_source(None, True) == DEFAULT_CONCAT_LABEL
        assert prioritize_tree_source("none", True) == DEFAULT_CONCAT_LABEL

    def test_no_coalescent_no_concat_is_none(self):
        assert prioritize_tree_source(None, False) == "none"
        assert prioritize_tree_source("none", False) == "none"

    def test_degenerate_first_gene_tree_without_concat(self):
        assert prioritize_tree_source("first_gene_tree", False) == "first_gene_tree"

    def test_degenerate_first_gene_tree_prefers_concat(self):
        # A concatenation tree is more reliable than the single-gene-tree fallback.
        assert prioritize_tree_source("first_gene_tree", True) == DEFAULT_CONCAT_LABEL

    def test_concat_label_override(self):
        assert prioritize_tree_source(None, True, concat_label="concat") == "concat"
        assert prioritize_tree_source("first_gene_tree", True, concat_label="concat") == "concat"

    def test_coalescent_is_degenerate(self):
        assert coalescent_is_degenerate("first_gene_tree") is True
        assert coalescent_is_degenerate("astral") is False
        assert coalescent_is_degenerate(None) is False


class TestParseFasta:
    def test_parses_first_token_as_id(self, tmp_path):
        p = tmp_path / "x.faa"
        p.write_text(">seq1 description here\nACGT\nACGT\n>seq2\nTGCA\n")
        seqs = parse_fasta(str(p))
        assert seqs == {"seq1": "ACGTACGT", "seq2": "TGCA"}

    def test_missing_file_returns_empty(self):
        assert parse_fasta("/no/such/file.faa") == {}

    def test_empty_sequence_allowed(self, tmp_path):
        p = tmp_path / "y.faa"
        p.write_text(">seq1\n\n>seq2\nACGT\n")
        seqs = parse_fasta(str(p))
        assert seqs == {"seq1": "", "seq2": "ACGT"}
