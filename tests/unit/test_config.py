"""Unit tests for configuration dataclasses."""

import pytest

from markerfinder.config import (
    AlignerConfig,
    HGTConfig,
    HeterogeneityConfig,
    MAGConfig,
    OrthologConfig,
    PhylogeneticConfig,
    PipelineConfig,
    ReportConfig,
    SelectionConfig,
)
from markerfinder.models.marker import SelectionStrategy


class TestSelectionConfig:
    def test_defaults(self):
        c = SelectionConfig()
        assert c.min_occupancy == 0.75
        assert c.max_markers == 60
        assert c.strategy == SelectionStrategy.GREEDY

    def test_custom(self):
        c = SelectionConfig(min_occupancy=0.5, max_markers=100)
        assert c.min_occupancy == 0.5
        assert c.max_markers == 100


class TestHGTConfig:
    def test_defaults(self):
        c = HGTConfig()
        assert c.enable_phylogenetic is True

    def test_level_thresholds(self):
        c = HGTConfig()
        assert c.level_thresholds['level1_max'] == 0.25
        assert c.level_thresholds['level2_max'] == 0.60


class TestMAGConfig:
    def test_defaults(self):
        c = MAGConfig()
        assert c.min_completeness == 50.0
        assert c.quality_weighted is True


class TestAlignerConfig:
    def test_defaults(self):
        c = AlignerConfig()
        assert c.cpus == 1
        assert c.output_prefix == 'markerfinder'


class TestPhylogeneticConfig:
    def test_defaults(self):
        c = PhylogeneticConfig()
        assert c.ufboot_replicates == 1000
        assert c.coalescent_mode == 'post-filter'
        assert c.use_fasttree is True


class TestReportConfig:
    def test_defaults(self):
        c = ReportConfig()
        assert c.report_format == 'html'


class TestOrthologConfig:
    def test_defaults(self):
        c = OrthologConfig()
        assert c.bbh_evalue == 1e-10
        assert c.bbh_identity == 30.0
        assert c.use_diamond is True


class TestHeterogeneityConfig:
    def test_defaults(self):
        c = HeterogeneityConfig()
        assert c.min_allele_freq == 0.2
        assert c.identity_threshold == 0.99


class TestPipelineConfig:
    def test_defaults(self):
        c = PipelineConfig()
        assert c.mode == 'standard'
        assert c.verify_db is True

    def test_to_dict(self):
        c = PipelineConfig(input_dir="/tmp/in", output_dir="/tmp/out")
        d = c.to_dict()
        assert d['input_dir'] == '/tmp/in'
        assert d['output_dir'] == '/tmp/out'
        assert 'mag_config' in d

    def test_sub_configs(self):
        c = PipelineConfig()
        assert isinstance(c.mag_config, MAGConfig)
        assert isinstance(c.selection_config, SelectionConfig)
        assert isinstance(c.hgt_config, HGTConfig)
        assert isinstance(c.aligner_config, AlignerConfig)
        assert isinstance(c.phylo_config, PhylogeneticConfig)
        assert isinstance(c.report_config, ReportConfig)


class TestOneAnswerPerRunLocation:
    """``PipelineConfig`` must be the single source of "where does this run write".

    Each section dataclass used to keep its own ``output_dir`` / ``tmp_dir`` /
    ``output_prefix`` default, so building a config as a library user does —
    ``PipelineConfig(output_dir="run1")`` — sent the reports, trees and partition
    file to ``./output`` (the ReportConfig default) while the run's own snapshot
    went to ``run1``. Tests that did this quietly littered the repository root.
    """

    def test_sections_inherit_the_run_locations(self, tmp_path):
        out = str(tmp_path / "run1")
        cfg = PipelineConfig(output_dir=out, tmp_dir=str(tmp_path / "scratch"),
                             input_dir=str(tmp_path))
        for name, section in (("report_config", cfg.report_config),
                              ("phylo_config", cfg.phylo_config),
                              ("selection_config", cfg.selection_config),
                              ("aligner_config", cfg.aligner_config),
                              ("mag_config", cfg.mag_config)):
            for field in ("output_dir", "tmp_dir", "input_dir"):
                if hasattr(section, field):
                    assert getattr(section, field) == {
                        "output_dir": out,
                        "tmp_dir": str(tmp_path / "scratch"),
                        "input_dir": str(tmp_path),
                    }[field], f"{name}.{field} disagrees with the pipeline"

    def test_values_set_after_construction_are_respected(self, tmp_path):
        """Propagation happens at construction, not on every read."""
        cfg = PipelineConfig(output_dir=str(tmp_path / "run1"))
        cfg.report_config.output_dir = str(tmp_path / "elsewhere")
        assert cfg.to_dict()["report_config"]["output_dir"] == (
            str(tmp_path / "elsewhere"))

    def test_control_default_config_keeps_one_output_dir(self):
        """The shipped default is still ``./output``; what must not happen is a
        section disagreeing with it, because that is what scattered files into a
        directory the caller never named."""
        cfg = PipelineConfig()
        sections = [name for name in dir(cfg)
                    if name.endswith("_config")
                    and hasattr(getattr(cfg, name), "output_dir")]
        assert sections, "no section carries output_dir: this control is vacuous"
        for name in sections:
            assert getattr(cfg, name).output_dir == cfg.output_dir, name
