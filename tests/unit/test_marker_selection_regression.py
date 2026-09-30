"""Code-review regression test for ``modules/marker_selection.py`` extract_marker_sequences (#6 D39).

The fix for #6 D39 was documentation-only (the homologous multi-copy singleton
selection strategy); behaviour is unchanged. This regression test proves the
function still returns the expected per-genome representative sequences:
  - one entry per (marker, genome), keyed by the bare genome id;
  - when a genome hits several homologous copies (multi-copy), only the
    highest-bitscore representative is kept (the singleton strategy).

``hmmsearch`` is mocked so the test needs no HMMER binary.
"""

from unittest.mock import patch

import pytest

from markerfinder.config import SelectionConfig
from markerfinder.models.genome import GeneState, Genome, OccupancyMatrix
from markerfinder.models.marker import SelectedMarkerSet
from markerfinder.modules.marker_selection import AdaptiveMarkerSelectionModule

FAKE_HITS = {
    "prot1": [
        {
            "target": "prot1",
            "score": 100.0,
            "evalue": 1e-10,
            "ali_start": 1,
            "ali_end": 10,
        }
    ],
}

MULTI_HITS = {
    "prot1": [
        {
            "target": "prot1",
            "score": 100.0,
            "evalue": 1e-10,
            "ali_start": 1,
            "ali_end": 10,
        }
    ],
    "prot2": [
        {
            "target": "prot2",
            "score": 50.0,
            "evalue": 1e-5,
            "ali_start": 1,
            "ali_end": 10,
        }
    ],
}


@pytest.fixture
def env(tmp_path):
    # A (dummy) HMM profile so hmm_path_for resolves a path for COG1.
    (tmp_path / "COG1.hmm").write_text("dummy")
    faa = tmp_path / "g1.faa"
    faa.write_text(">prot1\nMKALLIILFQAVLT\n>prot2\nAAAAAAAAAAAAAAAAAAAA\n")
    genome = Genome(
        id="g1", fasta_path=str(faa), protein_fasta_path=str(faa)
    )
    matrix = OccupancyMatrix(genomes=["g1"], cogs=["COG1"])
    matrix.set("g1", "COG1", GeneState.SINGLE_COPY)
    marker_set = SelectedMarkerSet(markers=["COG1"])
    mod = AdaptiveMarkerSelectionModule(
        SelectionConfig(marker_hmm_dir=str(tmp_path))
    )
    return mod, genome, matrix, marker_set, str(tmp_path)


class TestExtractMarkerSequences:
    def test_returns_expected_sequences(self, env):
        mod, genome, matrix, marker_set, td = env
        with patch(
            "markerfinder.modules.marker_selection._run_hmmsearch",
            return_value=FAKE_HITS,
        ):
            mg = mod.extract_marker_sequences(
                [genome], marker_set, matrix, tmp_dir=td
            )
        assert "COG1" in mg
        seqs = mg["COG1"]
        assert len(seqs) == 1
        entry = seqs[0]
        assert entry["id"] == "g1"
        assert entry["target"] == "prot1"
        # Ali_start=1, ali_end=10 -> first 10 residues of prot1.
        assert entry["seq"] == "MKALLIILFQ"

    def test_single_representative_for_multi_copy(self, env):
        mod, genome, matrix, marker_set, td = env
        with patch(
            "markerfinder.modules.marker_selection._run_hmmsearch",
            return_value=MULTI_HITS,
        ):
            mg = mod.extract_marker_sequences(
                [genome], marker_set, matrix, tmp_dir=td
            )
        seqs = mg["COG1"]
        # Only the highest-scoring representative (prot1, score 100) is kept.
        assert len(seqs) == 1
        assert seqs[0]["target"] == "prot1"
