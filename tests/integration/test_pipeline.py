"""End-to-end integration tests for MarkerFinderPipeline.

These tests run the full pipeline with all external tools (MAFFT, trimAl,
FastTree, IQ-TREE3, ASTRAL-III, HMMER) mocked via ``subprocess.run``.
They use tiny synthetic datasets and are marked with ``@pytest.mark.integration``
so they can be skipped with ``pytest -m "not integration"`` if desired.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from markerfinder.config import (
    HGTConfig,
    PhylogeneticConfig,
    PipelineConfig,
    SelectionConfig,
    TaxonomyConfig,
)
from markerfinder.models.genome import Genome
from markerfinder.pipeline import MarkerFinderPipeline
from markerfinder.utils.io import load_genomes_from_directory


pytestmark = pytest.mark.integration


def _aln_fasta(ids):
    seq = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(f">{sid}\n{seq}\n" for sid in ids)


def _tree_newick(ids):
    if len(ids) == 2:
        return f"({ids[0]},{ids[1]});"
    return f"(({ids[0]},{ids[1]}),({ids[2]},{ids[3]}));"


def _domtblout(genome_fasta: Path, marker_id: str, score: float = 100.0):
    """Create a minimal hmmsearch domtblout with one hit (single-copy)."""
    lines = ["# hmmsearch domtblout"]
    for line in genome_fasta.read_text().splitlines():
        if line.startswith(">"):
            seq_id = line[1:].split()[0]
            lines.append(
                f"{seq_id}\t-\t{seq_id}\t-\t0\t0\t1e-10\t{score}\t"
                f"0\t0\t1e-10\t{score}\t1\t100\t1\t20\t1\t20\t0\t0\t0\t-"
            )
            break  # Only one hit per genome so it is counted as single-copy
    return "\n".join(lines) + "\n"


def _mock_subprocess_run(tmp_path: Path, genome_ids):
    """Return a side-effect for ``subprocess.run`` that mimics external tools."""
    ids = genome_ids

    def runner(cmd, **kwargs):
        exe = cmd[0].lower()

        if exe == "hmmsearch":
            domtblout = cmd[cmd.index("--domtblout") + 1]
            fasta_path = Path(cmd[-1])
            Path(domtblout).write_text(_domtblout(fasta_path, cmd[-2]), encoding="utf-8")
            return MagicMock(returncode=0)

        if "mafft" in exe:
            out = kwargs.get("stdout")
            aln = _aln_fasta(ids)
            if out is not None:
                out.write(aln)
            return MagicMock(returncode=0)

        if exe == "trimal":
            in_path = Path(cmd[cmd.index("-in") + 1])
            out_path = Path(cmd[cmd.index("-out") + 1])
            if in_path.exists():
                out_path.write_text(in_path.read_text(), encoding="utf-8")
            return MagicMock(returncode=0)

        if exe == "fasttree" or cmd[0] == "FastTree":
            out = cmd[cmd.index("-out") + 1]
            Path(out).write_text(_tree_newick(ids), encoding="utf-8")
            return MagicMock(returncode=0)

        if "iqtree" in exe or exe == "iqtree3":
            prefix = cmd[cmd.index("-pre") + 1]
            Path(f"{prefix}.treefile").write_text(_tree_newick(ids), encoding="utf-8")
            return MagicMock(returncode=0)

        if exe == "astral":
            out = cmd[cmd.index("-o") + 1]
            Path(out).write_text(_tree_newick(ids), encoding="utf-8")
            return MagicMock(returncode=0)

        # Default: success without side effects.
        return MagicMock(returncode=0)

    return runner


@pytest.fixture
def genome_ids():
    return ["g1", "g2", "g3", "g4"]


@pytest.fixture
def make_input_dir(tmp_path, genome_ids):
    def _make():
        indir = tmp_path / "input"
        indir.mkdir()
        for gid in genome_ids:
            (indir / f"{gid}.faa").write_text(
                f">p1\nACDEFGHIKLMNPQRSTVWY\n>p2\nACDEFGHIKLMNPQRSTVWY\n",
                encoding="utf-8",
            )
        return indir
    return _make


class TestFullPipelineGtdbTk:
    def test_gtdb_tk_pipeline_runs(self, tmp_path, make_input_dir, genome_ids):
        indir = make_input_dir()
        gtdb_dir = tmp_path / "gtdb_markers"
        gtdb_dir.mkdir()
        for mid in ("M1", "M2"):
            lines = [f">{gid}\nACDEFGHIKLMNPQRSTVWY" for gid in genome_ids]
            (gtdb_dir / f"{mid}.faa").write_text("\n".join(lines) + "\n", encoding="utf-8")

        species_tree = tmp_path / "tree.nwk"
        species_tree.write_text(_tree_newick(genome_ids), encoding="utf-8")

        cfg = PipelineConfig(
            input_dir=str(indir),
            output_dir=str(tmp_path / "output"),
            tmp_dir=str(tmp_path / "tmp"),
            output_prefix="markerfinder",
            selection_config=SelectionConfig(
                marker_mode="gtdb_tk",
                gtdb_markers_dir=str(gtdb_dir),
                species_tree=str(species_tree),
                min_occupancy=0.0,
                max_markers=10,
            ),
        )
        cfg.phylo_config.fast_mode = True
        cfg.phylo_config.use_fasttree = True
        cfg.phylo_config.coalescent_mode = "post-filter"
        cfg.report_config.output_dir = cfg.output_dir
        cfg.report_config.output_prefix = cfg.output_prefix

        pipeline = MarkerFinderPipeline(cfg)
        genomes = load_genomes_from_directory(str(indir))

        with patch("subprocess.run", side_effect=_mock_subprocess_run(tmp_path, genome_ids)):
            with patch("markerfinder.utils.gtdb_tk_markers.shutil.which", return_value="/usr/bin/tool"):
                result = pipeline.run(genomes)

        assert result.usable_marker_count > 0
        assert result.species_tree_source == "astral"
        out = Path(cfg.output_dir)
        assert (out / "Phase5_reports" / "markerfinder.report.html").exists()
        assert (out / "Phase5_reports" / "markerfinder.pipeline_summary.txt").exists()
        assert (out / "Phase4_trees" / "markerfinder.species_tree_concat.newick").exists()
        assert (out / "Phase4_trees" / "markerfinder.species_tree_astral.newick").exists()
        assert (out / "Phase4_alignments" / "markerfinder.partition.nex").exists()
        assert (out / "Phase5_metadata" / "run_config.json").exists()


class TestFullPipelineHmm:
    def test_hmm_pipeline_runs(self, tmp_path, make_input_dir, genome_ids):
        indir = make_input_dir()
        hmm_dir = tmp_path / "hmms"
        hmm_dir.mkdir()
        for mid in ("M1", "M2"):
            (hmm_dir / f"{mid}.hmm").write_text("HMMER3/b [3.1b2]\n", encoding="utf-8")

        taxonomy = tmp_path / "taxa.tsv"
        taxonomy.write_text(
            "g1\td__Bacteria;p__P;c__C1;o__O;f__F;g__G1\n"
            "g2\td__Bacteria;p__P;c__C1;o__O;f__F;g__G1\n"
            "g3\td__Bacteria;p__P;c__C2;o__O;f__F;g__G2\n"
            "g4\td__Bacteria;p__P;c__C2;o__O;f__F;g__G2\n",
            encoding="utf-8",
        )

        cfg = PipelineConfig(
            input_dir=str(indir),
            output_dir=str(tmp_path / "output"),
            tmp_dir=str(tmp_path / "tmp"),
            output_prefix="markerfinder",
            selection_config=SelectionConfig(
                marker_mode="hmm",
                marker_hmm_dir=str(hmm_dir),
                min_occupancy=0.0,
                max_markers=10,
                min_hmm_score=20.0,
            ),
            taxonomy_config=TaxonomyConfig(
                taxonomy_table=str(taxonomy),
                taxonomy_format="table",
            ),
        )
        cfg.phylo_config.fast_mode = True
        cfg.phylo_config.use_fasttree = True
        cfg.phylo_config.coalescent_mode = "post-filter"
        cfg.report_config.output_dir = cfg.output_dir
        cfg.report_config.output_prefix = cfg.output_prefix

        pipeline = MarkerFinderPipeline(cfg)
        genomes = load_genomes_from_directory(str(indir))
        table_taxa = {
            gid: {
                "domain": "Bacteria", "phylum": "P", "class": "C1" if gid in ("g1", "g2") else "C2",
                "order": "O", "family": "F", "genus": "G1" if gid in ("g1", "g2") else "G2",
            }
            for gid in genome_ids
        }

        with patch("subprocess.run", side_effect=_mock_subprocess_run(tmp_path, genome_ids)):
            result = pipeline.run(genomes, table_taxa=table_taxa, taxonomy_origin="external table")

        assert result.usable_marker_count > 0
        assert result.species_tree_source == "astral"
        out = Path(cfg.output_dir)
        assert (out / "Phase5_reports" / "markerfinder.report.html").exists()
        assert (out / "Phase5_reports" / "markerfinder.hgt_evaluation.tsv").exists()
