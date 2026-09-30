import pytest
from markerfinder.models.genome import (
    Genome, GenomeSource, GenomeQuality, GeneState, OccupancyMatrix,
)
from markerfinder.models.marker import (
    MarkerLevel, MarkerQualityScore, SelectedMarkerSet, SelectionStrategy,
    ResolutionPreset, RESOLUTION_PRESETS,
)
from markerfinder.models.alignment import PartitionEntry, PartitionFile, AlignmentResult
from markerfinder.models.tree import Tree, TreeRecommendation
from markerfinder.models.report import PipelineResult, PhaseContext, RuntimeInfo
from markerfinder.models.pipeline_types import (
    NoMarkerAvailableError, BlastHit,
    PhyloStepResult, ConflictType, HGTEvaluation, HGTReport,
)


class TestGenomeSource:
    def test_values(self):
        assert GenomeSource.ISOLATE.value == "isolate"
        assert GenomeSource.MAG.value == "mag"


class TestGenome:
    def test_creation(self):
        g = Genome(id="g1", fasta_path="g1.faa")
        assert g.id == "g1"
        assert g.source_type == GenomeSource.ISOLATE

    def test_mag(self):
        g = Genome(id="m1", fasta_path="m1.faa", source_type=GenomeSource.MAG)
        assert g.source_type == GenomeSource.MAG

    def test_codon_usage_from_fasta(self, tmp_path):
        fa = tmp_path / "test.faa"
        fa.write_text(">seq1\nATGAAATTT\n>seq2\nATGCCCGGG\n")
        g = Genome(id="g1", fasta_path=str(fa))
        codon = g.get_background_codon_usage()
        assert len(codon) > 0
        assert all(v > 0 for v in codon.values())

    def test_tetranucleotide_from_fasta(self, tmp_path):
        fa = tmp_path / "test.faa"
        fa.write_text(">seq1\nATGAAATTTCCC\n>seq2\nATGCCCGGGAAA\n")
        g = Genome(id="g1", fasta_path=str(fa))
        tetra = g.get_tetranucleotide_signature()
        assert len(tetra) > 0

    def test_codon_usage_missing_file(self):
        g = Genome(id="g1", fasta_path="/nonexistent/path.faa")
        assert g.get_background_codon_usage() == {}

    def test_tetranucleotide_missing_file(self):
        g = Genome(id="g1", fasta_path="/nonexistent/path.faa")
        assert g.get_tetranucleotide_signature() == []


class TestGenomeQuality:
    def test_default(self):
        gq = GenomeQuality()
        assert gq.completeness == 0.0


class TestOccupancyMatrix:
    def test_creation(self):
        m = OccupancyMatrix(genomes=["g1"], cogs=["C1"])
        assert m.n_genomes == 1
        assert m.n_cogs == 1

    def test_get_set(self):
        m = OccupancyMatrix(genomes=["g1"], cogs=["C1"])
        assert m.get("g1", "C1") == GeneState.ABSENT
        m.set("g1", "C1", GeneState.SINGLE_COPY)
        assert m.get("g1", "C1") == GeneState.SINGLE_COPY


class TestTree:
    def test_creation(self):
        t = Tree(newick="(A,B);")
        assert t.to_newick() == "(A,B);"

    def test_read(self, tmp_path):
        p = tmp_path / "t.nwk"
        p.write_text("(A,B,C);")
        t = Tree.read(str(p))
        assert t.newick == "(A,B,C);"


class TestPartitionFile:
    def test_nexus(self, tmp_path):
        pf = PartitionFile(entries=[PartitionEntry(name="C1", start=1, end=300)])
        out = tmp_path / "p.nex"
        pf.write_nexus(str(out))
        assert "charset C1 = 1-300;" in out.read_text()


class TestHGTEvaluation:
    def test_default(self):
        ev = HGTEvaluation()
        assert ev.overall_risk == 0.0
        assert ev.level == MarkerLevel.LEVEL_1


class TestPipelineResult:
    def test_default(self):
        pr = PipelineResult()
        assert pr.runtime.duration == 0.0
