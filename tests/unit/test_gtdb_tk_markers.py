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


class TestGeneTreeFastTreeUnavailable:
    """The "FastTree unavailable" branch must degrade, not crash.

    The summary warning interpolated ``ft_last_err`` -- a name that was never
    assigned on any code path. So the branch taken whenever fasttree is simply
    not installed, which is the most ordinary reason to lack a gene tree, raised
    ``NameError`` from inside the logger call whose whole job was to describe
    that degradation. A recoverable loss became a crash.

    Must-fail control: with ``ft_last_err`` unbound, the assertion on "not on
    PATH" below is unreachable because the call raises first.
    """

    @staticmethod
    def _seqs(n=4):
        return [
            {"genome_id": f"G{i}", "seq": "ACDEFGHIKLMNPQRSTVWY" * 2}
            for i in range(n)
        ]

    def _run_without_fasttree(self, tmp_path, monkeypatch, caplog):
        import subprocess as real_subprocess

        from markerfinder.utils import gtdb_tk_markers as g

        marker_id = "M1"
        trimmed = tmp_path / f"{marker_id}.aln.trim.faa"

        def fake_run(cmd, **kwargs):
            exe = cmd[0]
            if exe == "mafft":
                # MAFFT writes the alignment into the file object it was given.
                kwargs["stdout"].write(">G0\nACDE\n")
                return real_subprocess.CompletedProcess(cmd, 0)
            if exe == "trimal":
                trimmed.write_text(">G0\nACDE\n", encoding="utf-8")
                return real_subprocess.CompletedProcess(cmd, 0)
            # Both "fasttree" and "FastTree" are absent.
            raise FileNotFoundError(2, "No such file or directory", exe)

        monkeypatch.setattr(g.subprocess, "run", fake_run)
        with caplog.at_level("WARNING", logger=g.logger.name):
            return g._build_one_gene_tree(marker_id, self._seqs(), str(tmp_path))

    def test_missing_fasttree_returns_no_tree_instead_of_raising(
        self, tmp_path, monkeypatch, caplog
    ):
        marker_id, path = self._run_without_fasttree(tmp_path, monkeypatch, caplog)
        assert marker_id == "M1"
        assert path is None, "no tree can be produced without FastTree"

    def test_the_degradation_names_the_real_cause(self, tmp_path, monkeypatch, caplog):
        self._run_without_fasttree(tmp_path, monkeypatch, caplog)
        warned = [
            r.getMessage() for r in caplog.records if "FastTree unavailable" in r.getMessage()
        ]
        assert warned, (
            "a skipped gene tree must be reported, got: "
            f"{[r.getMessage() for r in caplog.records]}"
        )
        assert "not on PATH" in warned[0], warned[0]
