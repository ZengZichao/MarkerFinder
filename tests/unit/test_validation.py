"""Unit tests for input validation module."""

import pytest
import tempfile
from pathlib import Path

from markerfinder.validation import (
    validate_newick_string,
    detect_tree_format,
    count_trees_in_file,
    parse_nhx_annotations,
    strip_nhx_annotations,
    detect_sequence_format,
    validate_sequence_file,
    validate_tree_file,
    load_tree,
    cross_validate,
    ValidationReport,
)
from markerfinder.exceptions import TreeValidationError, SequenceValidationError


class TestNewickValidation:
    def test_valid_simple(self):
        errors, _ = validate_newick_string("((A:0.1,B:0.2),(C:0.3,D:0.4));")
        critical = [e for e in errors if "Negative" in e or "Duplicate tip" in e]
        assert not critical

    def test_unbalanced_parens(self):
        errors, _ = validate_newick_string("((A,B);")
        assert any("Unbalanced" in e for e in errors)

    def test_negative_branch(self):
        errors, _ = validate_newick_string("((A:-1,B),(C,D));")
        assert any("Negative" in e for e in errors)

    def test_duplicate_tips(self):
        errors, _ = validate_newick_string("((A,B),(A,D));")
        assert any("Duplicate tip" in e for e in errors)

    def test_empty_string(self):
        errors, _ = validate_newick_string("")
        assert len(errors) > 0


class TestTreeFormatDetection:
    def test_newick_extension(self, tmp_path):
        f = tmp_path / "tree.nwk"
        f.write_text("(A,B);", encoding="utf-8")
        assert detect_tree_format(str(f)) == "newick"

    def test_nexus_extension(self, tmp_path):
        f = tmp_path / "tree.nex"
        f.write_text("#NEXUS\n", encoding="utf-8")
        assert detect_tree_format(str(f)) == "nexus"

    def test_nexus_content(self, tmp_path):
        f = tmp_path / "tree.txt"
        f.write_text("#NEXUS\ntree t1 = (A,B);", encoding="utf-8")
        assert detect_tree_format(str(f)) == "nexus"


class TestCountTrees:
    def test_single_newick(self, tmp_path):
        f = tmp_path / "tree.nwk"
        f.write_text("(A,B);\n", encoding="utf-8")
        assert count_trees_in_file(str(f), "newick") == 1

    def test_multiple_newick(self, tmp_path):
        f = tmp_path / "trees.nwk"
        f.write_text("(A,B);\n(C,D);\n", encoding="utf-8")
        assert count_trees_in_file(str(f), "newick") == 2


class TestNHXAnnotations:
    def test_parse_nhx(self):
        newick = "((A[&&NHX:B=95:taxid=123],B),C);"
        ann = parse_nhx_annotations(newick)
        assert len(ann) > 0

    def test_strip_nhx(self):
        newick = "((A[&&NHX:B=95],B),C);"
        stripped = strip_nhx_annotations(newick)
        assert "&&NHX" not in stripped


class TestSequenceFormat:
    def test_fasta(self, tmp_path):
        f = tmp_path / "seq.fasta"
        f.write_text(">seq1\nACGT\n", encoding="utf-8")
        assert detect_sequence_format(str(f)) == "fasta"

    def test_fastq_extension(self, tmp_path):
        f = tmp_path / "seq.fq"
        f.write_text("@seq1\nACGT\n+\n!!!!\n", encoding="utf-8")
        assert detect_sequence_format(str(f)) == "fastq"


class TestSequenceValidation:
    def test_valid_fasta(self, tmp_path):
        f = tmp_path / "seq.fasta"
        f.write_text(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n", encoding="utf-8")
        seqs, errors = validate_sequence_file(str(f), mol_type="DNA")
        assert "seq1" in seqs
        assert "seq2" in seqs

    def test_duplicate_ids(self, tmp_path):
        f = tmp_path / "dup.fasta"
        f.write_text(">seq1\nACGT\n>seq1\nTGCA\n", encoding="utf-8")
        with pytest.raises(SequenceValidationError, match="Duplicate"):
            validate_sequence_file(str(f))

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.fasta"
        f.write_text("", encoding="utf-8")
        with pytest.raises(SequenceValidationError, match="empty"):
            validate_sequence_file(str(f))

    def test_invalid_alphabet(self, tmp_path):
        f = tmp_path / "bad.fasta"
        f.write_text(">seq1\nACGTXYZ\n", encoding="utf-8")
        seqs, errors = validate_sequence_file(str(f), mol_type="DNA")
        assert any("Invalid" in e for e in errors)

    def test_length_inconsistency(self, tmp_path):
        f = tmp_path / "uneven.fasta"
        f.write_text(">seq1\nACGT\n>seq2\nACGTACGT\n", encoding="utf-8")
        seqs, errors = validate_sequence_file(str(f), mol_type="DNA")
        assert any("length" in e.lower() or "consistent" in e.lower() for e in errors)


class TestTreeFileValidation:
    def test_valid_file(self, tmp_path):
        f = tmp_path / "tree.nwk"
        f.write_text("((A:0.1,B:0.2),(C:0.3,D:0.4));", encoding="utf-8")
        newick, errors = validate_tree_file(str(f))
        assert "A" in newick

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.nwk"
        f.write_text("", encoding="utf-8")
        with pytest.raises(TreeValidationError, match="empty"):
            validate_tree_file(str(f))

    def test_negative_branch_file(self, tmp_path):
        f = tmp_path / "neg.nwk"
        f.write_text("((A:-1,B),(C,D));", encoding="utf-8")
        with pytest.raises(TreeValidationError, match="Negative"):
            validate_tree_file(str(f))

    def test_not_found(self):
        with pytest.raises(FileNotFoundError):
            validate_tree_file("/nonexistent/tree.nwk")

    def test_directory_instead_of_file(self, tmp_path):
        with pytest.raises(TreeValidationError, match="directory"):
            validate_tree_file(str(tmp_path))


class TestLoadTree:
    def test_load_valid(self, tmp_path):
        f = tmp_path / "tree.nwk"
        f.write_text("((A,B),(C,D));", encoding="utf-8")
        tree = load_tree(str(f), validate=True)
        assert tree.newick != ""

    def test_load_no_validate(self, tmp_path):
        f = tmp_path / "tree.nwk"
        f.write_text("((A,B),(C,D));", encoding="utf-8")
        tree = load_tree(str(f), validate=False)
        assert "A" in tree.newick


class TestValidationReport:
    def test_initial_state(self):
        r = ValidationReport()
        assert r.is_valid is True
        assert len(r.tree_tips) == 0
