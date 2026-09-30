import pytest
from pathlib import Path

from markerfinder.config import PipelineConfig, SelectionConfig
from markerfinder.models.genome import Genome, GenomeSource
from markerfinder.models.report import PhaseContext
from markerfinder.pipeline import MarkerFinderPipeline
from markerfinder.utils.io import load_genomes_from_directory


class TestMarkerFinderPipeline:
    def test_init(self):
        cfg = PipelineConfig()
        p = MarkerFinderPipeline(cfg)
        assert p.config is not None

    def test_run_without_gtdb_dir_raises(self, tmp_path):
        """Without --gtdb-markers-dir, gtdb_tk path raises ValueError."""
        for i in range(3):
            (tmp_path / f"g{i}.faa").write_text(f">seq{i}\nACGTACGT\n")
        cfg = PipelineConfig(
            input_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            tmp_dir=str(tmp_path / "tmp"),
            selection_config=SelectionConfig(min_occupancy=0.0, max_markers=5),
        )
        p = MarkerFinderPipeline(cfg)
        genomes = load_genomes_from_directory(str(tmp_path))
        with pytest.raises(ValueError):
            p.run_scan(genomes, PhaseContext())
