import pytest
from pathlib import Path

from markerfinder.utils.io import load_genomes_from_directory, write_fasta
from markerfinder.models.genome import GenomeSource


class TestLoadGenomes:
    def test_load_faa(self, tmp_path):
        (tmp_path / "g1.faa").write_text(">s1\nACGT\n")
        (tmp_path / "g2.faa").write_text(">s2\nTGCA\n")
        gs = load_genomes_from_directory(str(tmp_path))
        assert len(gs) == 2

    def test_detect_mag(self, tmp_path):
        (tmp_path / "mag_bin1.faa").write_text(">s1\nACGT\n")
        gs = load_genomes_from_directory(str(tmp_path))
        assert gs[0].source_type == GenomeSource.MAG

    def test_empty_dir(self, tmp_path):
        empty = tmp_path / "empty_subdir"
        empty.mkdir()
        gs = load_genomes_from_directory(str(empty))
        assert gs == []

    def test_nonexistent(self):
        with pytest.raises(FileNotFoundError):
            load_genomes_from_directory("/no/such/dir")


class TestWriteFasta:
    def test_basic(self, tmp_path):
        out = tmp_path / "o.fasta"
        write_fasta({"s1": "ACGT"}, str(out))
        assert ">s1" in out.read_text()

    def test_wrap(self, tmp_path):
        out = tmp_path / "o.fasta"
        write_fasta({"s1": "A" * 200}, str(out))
        lines = [l for l in out.read_text().split("\n") if l and not l.startswith(">")]
        for l in lines:
            assert len(l) <= 80
