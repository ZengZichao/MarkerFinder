"""Unit tests for GTDB-TK marker loading utilities."""

import pytest
from pathlib import Path

from markerfinder.models.genome import GeneState
from markerfinder.utils.gtdb_tk_markers import (
    _parse_fasta,
    load_gtdb_markers,
    split_user_msa,
)


def write_fasta(path: Path, records: dict):
    lines = []
    for header, seq in records.items():
        lines.append(f">{header}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class TestParseFasta:
    def test_simple(self, tmp_path):
        f = tmp_path / "test.faa"
        write_fasta(f, {"G1": "ACDEF", "G2": "GHIKL"})
        seqs = _parse_fasta(str(f))
        assert seqs == {"G1": "ACDEF", "G2": "GHIKL"}

    def test_missing_file(self):
        seqs = _parse_fasta("/nonexistent/path.faa")
        assert seqs == {}


class TestLoadGtdbMarkers:
    def test_load_single_marker(self, tmp_path):
        d = tmp_path / "gtdb_markers"
        d.mkdir()
        write_fasta(d / "COG001.faa", {"G1": "ACDEF", "G2": "GHIKL"})

        matrix, mseq, occ = load_gtdb_markers(str(d), genome_ids=["G1", "G2"])

        assert matrix.cogs == ["COG001"]
        assert matrix.get("G1", "COG001") == GeneState.SINGLE_COPY
        assert matrix.get("G2", "COG001") == GeneState.SINGLE_COPY
        assert occ["COG001"] == 1.0
        assert len(mseq["COG001"]) == 2

    def test_missing_genome_reported_absent(self, tmp_path):
        d = tmp_path / "gtdb_markers"
        d.mkdir()
        write_fasta(d / "COG001.faa", {"G1": "ACDEF"})

        matrix, _, occ = load_gtdb_markers(str(d), genome_ids=["G1", "G2"])

        assert matrix.get("G1", "COG001") == GeneState.SINGLE_COPY
        assert matrix.get("G2", "COG001") == GeneState.ABSENT
        assert occ["COG001"] == 0.5

    def test_empty_directory_raises(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        with pytest.raises(FileNotFoundError):
            load_gtdb_markers(str(d))


class TestSplitUserMsa:
    def test_split(self, tmp_path):
        msa = tmp_path / "msa.faa"
        # Two genomes, two markers back-to-back (lengths 3 and 5).
        msa.write_text(
            ">G1\nABCDEFGH\n>G2\n12345678\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "split"
        split_user_msa(str(msa), [("M1", 3), ("M2", 5)], str(out_dir))

        m1 = (out_dir / "M1.faa").read_text(encoding="utf-8")
        m2 = (out_dir / "M2.faa").read_text(encoding="utf-8")
        assert ">G1\nABC\n>G2\n123" in m1
        assert ">G1\nDEFGH\n>G2\n45678" in m2
