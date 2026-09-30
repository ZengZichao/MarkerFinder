import pytest
from unittest.mock import patch, MagicMock

from markerfinder.config import MAGConfig
from markerfinder.models.genome import Genome, GenomeQuality, GenomeSource, OccupancyMatrix
from markerfinder.models.pipeline_types import PreprocessingResult
from markerfinder.modules.mag_optimization import (
    MAGOptimizationModule, MAGQualityPreprocessor, SparseSupermatrixOptimizer,
    MAGHeterogeneityHandler, _parse_checkm_results, _run_checkm,
)
from markerfinder.config import HeterogeneityConfig


class TestMAGQualityPreprocessor:
    def test_stratify(self):
        p = MAGQualityPreprocessor(MAGConfig())
        qr = {
            "iso": p._assess_quality([Genome(id="iso", fasta_path="x", source_type=GenomeSource.ISOLATE)])["iso"],
            "mag": p._assess_quality([Genome(id="mag", fasta_path="x", source_type=GenomeSource.MAG)])["mag"],
        }
        layers = p._stratify_genomes(qr)
        assert "iso" in layers["high"]

    def test_adaptive_params(self):
        p = MAGQualityPreprocessor(MAGConfig())
        params = p._compute_adaptive_parameters({"high": ["g1", "g2"], "medium": [], "low": [], "contaminated": []})
        assert params.min_occupancy == 0.75

    def test_quality_source_default_for_protein_input(self):
        p = MAGQualityPreprocessor(MAGConfig())
        result = p.preprocess([Genome(id="iso", fasta_path="x", source_type=GenomeSource.ISOLATE)])
        assert result.quality_results.quality_source == "default_no_nucleotide"

    def test_assess_quality_routes_to_checkm(self):
        p = MAGQualityPreprocessor(MAGConfig(input_dir="/tmp"))
        with patch("markerfinder.modules.mag_optimization._run_checkm") as mock_checkm:
            mock_checkm.return_value = {
                "g1": GenomeQuality(completeness=95.0, contamination=2.0, quality_score=85.0)
            }
            result = p._assess_quality([
                Genome(id="g1", fasta_path="g1.faa", nucleotide_fasta_path="g1.fna")
            ])
            mock_checkm.assert_called_once()
            assert result["g1"].completeness == 95.0

    def test_skip_checkm_uses_fallback(self):
        p = MAGQualityPreprocessor(MAGConfig(input_dir="/tmp", skip_checkm=True))
        with patch("markerfinder.modules.mag_optimization._run_checkm") as mock_checkm:
            result = p._assess_quality([
                Genome(id="g1", fasta_path="g1.faa", source_type=GenomeSource.ISOLATE)
            ])
            mock_checkm.assert_not_called()
            assert "__skipped__" in result
            assert result["g1"].completeness == 98.0

    def test_skip_checkm_quality_source(self):
        p = MAGQualityPreprocessor(MAGConfig(skip_checkm=True))
        result = p.preprocess([
            Genome(id="g1", fasta_path="g1.faa", source_type=GenomeSource.ISOLATE)
        ])
        assert result.quality_results.quality_source == "skipped"


class TestRunCheckMEmptyParse:
    def test_empty_parse_returns_sentinel_not_empty(self, tmp_path):
        """CheckM 命令成功但输出解析为空时, 不应返回 {} (质量信息全丢、
        自适应参数误判低质量), 而应带 __no_checkm__ sentinel 使用默认估计."""
        config = MAGConfig(input_dir="/tmp", tmp_dir=str(tmp_path))
        # 预创建空的 checkm 输出文件(解析结果为空 dict)
        (tmp_path / "checkm.tsv").write_text("")
        genome = Genome(
            id="g1", fasta_path="g1.fna",
            nucleotide_fasta_path="g1.fna", source_type=GenomeSource.MAG,
        )
        with patch("subprocess.run", return_value=MagicMock(returncode=0)):
            result = _run_checkm([genome], config)
        # 关键: 不再是空 dict —— 旧行为会静默丢失质量信息并误判低质量
        assert result != {}
        assert "__no_checkm__" in result
        assert "g1" in result
        # 使用默认质量估计(MAG -> 75/5/50), 而非把质量当作 0 / 缺失
        assert result["g1"].completeness > 0


class TestParseCheckMResults:
    def test_parse_valid_tsv(self, tmp_path):
        tsv = tmp_path / "checkm.tsv"
        tsv.write_text(
            "Bin Id\tMarker Lineage\t# Genomes\t# Markers\t# Marker Sets\t"
            "0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain Heterogeneity\n"
            "g1\td__Archaea\t1\t100\t50\t0\t0\t0\t0\t0\t0\t98.5\t1.2\t0.0\n"
        )
        genomes = [Genome(id="g1", fasta_path="g1.faa")]
        results = _parse_checkm_results(str(tsv), genomes)
        assert "g1" in results
        assert results["g1"].completeness == 98.5
        assert results["g1"].contamination == 1.2
        assert results["g1"].quality_score == 98.5 - 5 * 1.2


class TestMAGOptimizationModule:
    def test_run(self):
        mod = MAGOptimizationModule(MAGConfig())
        r = mod.run([Genome(id="g1", fasta_path="x")])
        assert isinstance(r, PreprocessingResult)


class TestMAGHeterogeneityHandler:
    def test_single_hit(self):
        h = MAGHeterogeneityHandler(HeterogeneityConfig())
        r = h.handle_marker_heterogeneity("m", "c", ["hit1"])
        assert r.n_alleles == 1
