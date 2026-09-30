from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from markerfinder.phases import canonical_log_tag
from markerfinder.config import HeterogeneityConfig, MAGConfig
from markerfinder.models.genome import GeneState, Genome, GenomeQuality, GenomeSource, OccupancyMatrix
from markerfinder.models.pipeline_types import (
    AdaptiveParams,
    HeterogeneityResult,
    OptimizedMatrix,
    PreprocessingResult,
    QualityData,
)
import logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# 默认质量估计 (无 CheckM 时使用) —— 单一事实源.
# 对应 ``_default_quality_estimates`` 中的硬编码魔法数字.
# ---------------------------------------------------------------------------
DEFAULT_MAG_COMPLETENESS = 75.0
DEFAULT_MAG_CONTAMINATION = 5.0
DEFAULT_MAG_QUALITY_SCORE = 50.0
DEFAULT_NON_MAG_COMPLETENESS = 98.0
DEFAULT_NON_MAG_CONTAMINATION = 1.0
DEFAULT_NON_MAG_QUALITY_SCORE = 93.0

# ---------------------------------------------------------------------------
# 自适应参数阈值 (MAGQualityPreprocessor._compute_adaptive_parameters) ——
# 单一事实源. 按 high / medium / low 三层质量分层计算 MAGOptimization 参数时
# 使用的硬编码魔法数字.
# ---------------------------------------------------------------------------
ADAPTIVE_HIGH_RATIO_THRESHOLD = 0.7    # High_ratio 超过该值 -> high 分支
ADAPTIVE_MEDIUM_RATIO_THRESHOLD = 0.3  # High_ratio 超过该值 -> medium 分支 (否则 low)
ADAPTIVE_HIGH_MIN_HMM_SCORE = 30.0
ADAPTIVE_HIGH_MIN_OCCUPANCY = 0.75
ADAPTIVE_HIGH_MAX_MARKERS = 60
ADAPTIVE_MEDIUM_MIN_HMM_SCORE = 20.0
ADAPTIVE_MEDIUM_MIN_OCCUPANCY = 0.55
ADAPTIVE_MEDIUM_MAX_MARKERS = 80
ADAPTIVE_LOW_MIN_HMM_SCORE = 15.0
ADAPTIVE_LOW_MIN_OCCUPANCY = 0.35
ADAPTIVE_LOW_MAX_MARKERS = 150


def _parse_checkm_results(
    quality_file: str, genomes: List[Genome]
) -> Dict[str, GenomeQuality]:
    """解析 CheckM (lineage_wf) --tablename 输出。"""
    results: Dict[str, GenomeQuality] = {}
    genome_ids = {g.id for g in genomes}

    try:
        with open(quality_file, encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # CheckM 以文件 stem 作为 Bin Id
                name = row.get("Bin Id", "")
                if name not in genome_ids:
                    continue
                completeness = float(row.get("Completeness", 0))
                contamination = float(row.get("Contamination", 0))
                quality_score = completeness - 5 * contamination
                results[name] = GenomeQuality(
                    completeness=completeness,
                    contamination=contamination,
                    quality_score=quality_score,
                )
    except Exception as e:
        logger.warning(f"  Failed to parse CheckM results: {e}")

    return results


def _no_checkm_fallback(genomes: List[Genome], sentinel: str = "__no_checkm__") -> Dict[str, GenomeQuality]:
    """无 CheckM 可用时的统一回退。"""
    out = _default_quality_estimates(genomes)
    out[sentinel] = GenomeQuality(completeness=0.0, contamination=0.0, quality_score=0.0)
    return out


def _run_checkm(genomes: List[Genome], config: MAGConfig) -> Dict[str, GenomeQuality]:
    """运行 CheckM 评估基因组质量。

    CheckM (lineage_wf) 依赖其参考数据集（HMM + 基因组树，约 1.4 GB），
    需通过 ``checkm data`` 命令预先下载并配置；它不是零配置工具。
    该步骤基于核酸 contig，因此需要 effective_nucleotide_path。
    纯蛋白输入时直接回退默认估计。
    """
    # 仅当所有 genome 都没有核酸路径时,才回退默认估计(蛋白-only 输入).
    has_nucleotide = any(g.effective_nucleotide_path is not None for g in genomes)
    if not has_nucleotide:
        logger.info(
            "  Protein-only input (.faa): skipping CheckM "
            "(requires nucleotide .fna). Using default quality estimates."
        )
        results: Dict[str, GenomeQuality] = {}
        results["__protein_input__"] = GenomeQuality(
            completeness=0.0, contamination=0.0, quality_score=0.0)
        default = _default_quality_estimates(genomes)
        results.update(default)
        return results

    if not config.input_dir:
        logger.warning(
            "  Nucleotide sequences available but input_dir not set; "
            "cannot run CheckM. Using default quality estimates."
        )
        return _no_checkm_fallback(genomes)

    output_dir = os.path.join(config.tmp_dir, "checkm_output")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    quality_file = os.path.join(config.tmp_dir, "checkm.tsv")

    cmd = [
        "checkm", "lineage_wf",
        "-t", str(config.cpus),
        "-x", "fna",
        "--tab_table",
        "-f", quality_file,
        config.input_dir,
        output_dir,
    ]

    try:
        logger.info("  Running CheckM...")
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=3600)
    except FileNotFoundError:
        logger.warning("  checkm not found. Using default quality estimates.")
        return _no_checkm_fallback(genomes)
    except subprocess.CalledProcessError as e:
        logger.warning(f"  CheckM failed: {e.stderr[:300]}")
        return _no_checkm_fallback(genomes)
    except subprocess.TimeoutExpired:
        logger.warning("  CheckM timed out")
        return _no_checkm_fallback(genomes)

    if Path(quality_file).exists():
        parsed = _parse_checkm_results(quality_file, genomes)
        # CheckM 命令成功(退出码 0)但解析为空(输出文件为空 / 格式不符 / 无匹配
        # Genome)时, 绝不能回退为 {} —— 空 dict 会让 preprocess 把质量信息"全丢",
        # 且自适应参数因总质量数为 0 而误判为低质量分支. 显式标记 unknown:
        # 带上 __no_checkm__ sentinel 并使用默认质量估计, 使下游走
        # Quality_source="default_no_checkm" 而非"低质量".
        if not parsed:
            logger.warning(
                "  CheckM ran but produced no parseable quality results; "
                "using default quality estimates (marked as unknown)."
            )
            return _no_checkm_fallback(genomes)
        return parsed

    logger.warning("  CheckM output not found. Using default quality estimates.")
    return _no_checkm_fallback(genomes)


def _default_quality_estimates(genomes: List[Genome]) -> Dict[str, GenomeQuality]:
    """默认质量估计（无 CheckM 时使用）。"""
    results: Dict[str, GenomeQuality] = {}
    for g in genomes:
        if g.source_type == GenomeSource.MAG:
            results[g.id] = GenomeQuality(
                completeness=DEFAULT_MAG_COMPLETENESS,
                contamination=DEFAULT_MAG_CONTAMINATION,
                quality_score=DEFAULT_MAG_QUALITY_SCORE,
            )
        else:
            results[g.id] = GenomeQuality(
                completeness=DEFAULT_NON_MAG_COMPLETENESS,
                contamination=DEFAULT_NON_MAG_CONTAMINATION,
                quality_score=DEFAULT_NON_MAG_QUALITY_SCORE,
            )
    return results


class MAGHeterogeneityHandler:
    """Reserved interface for MAG marker heterogeneity handling.

    Currently not used by the main pipeline (single best hit is kept per
    marker), but retained so future versions can resolve within-genome
    marker duplicates without breaking the public API/tests.
    """
    def __init__(self, config: HeterogeneityConfig):
        self.config = config

    def handle_marker_heterogeneity(self, mag_id, marker_id, raw_hits, assembly_graph=None):
        if len(raw_hits) <= 1:
            return HeterogeneityResult(
                selected_sequence=raw_hits[0] if raw_hits else None,
                heterogeneity_score=0.0, n_alleles=len(raw_hits), confidence="high",
            )
        return HeterogeneityResult(
            selected_sequence=raw_hits[0], heterogeneity_score=0.0,
            n_alleles=len(raw_hits), confidence="low",
        )


class MAGQualityPreprocessor:
    def __init__(self, config: MAGConfig):
        self.config = config

    def preprocess(self, genomes: List[Genome]) -> PreprocessingResult:
        logger.info(f"{canonical_log_tag('0')} Quality-Aware Preprocessing")
        quality_dict = self._assess_quality(genomes)
        # Quality_source 由 _assess_quality 返回的 sentinel 决定:
        # - '__precomputed__': 用户提供了 --checkm-results;
        # - '__protein_input__': 蛋白-only.faa 输入主动跳过 CheckM(默认估计);
        # - '__skipped__': 用户显式通过 --skip-checkm 跳过 CheckM;
        # - '__no_checkm__': 无 CheckM 结果可用(未安装/失败/无质量文件);
        # - 其他: 真实 CheckM 质量值.
        if "__precomputed__" in quality_dict:
            quality_dict.pop("__precomputed__")
            quality_source = "precomputed"
        elif "__protein_input__" in quality_dict:
            quality_dict.pop("__protein_input__")
            quality_source = "default_no_nucleotide"
        elif "__skipped__" in quality_dict:
            quality_dict.pop("__skipped__")
            quality_source = "skipped"
        elif "__no_checkm__" in quality_dict:
            quality_dict.pop("__no_checkm__")
            quality_source = "default_no_checkm"
        else:
            quality_source = "checkm"
        # Ensure each GenomeQuality carries its genome id for reporting.
        for gid, gq in quality_dict.items():
            if gq.genome_id is None:
                gq.genome_id = gid
        quality_data = QualityData(genome_qualities=list(quality_dict.values()),
                                  quality_source=quality_source)
        layers = self._stratify_genomes(quality_dict)
        adaptive_params = self._compute_adaptive_parameters(layers)

        logger.info(f"  Quality layers:")
        logger.info(f"    High quality  (>90%, <5%):  {len(layers['high'])}")
        logger.info(f"    Medium quality (70-90%, <10%): {len(layers['medium'])}")
        logger.info(f"    Low quality  (<70%): {len(layers['low'])}")
        logger.info(f"    Contaminated (>10%): {len(layers['contaminated'])}")

        return PreprocessingResult(
            quality_results=quality_data, adaptive_params=adaptive_params, layers=layers,
        )

    def _assess_quality(self, genomes: List[Genome]) -> Dict[str, GenomeQuality]:
        """评估基因组质量。

        优先使用用户提供的 CheckM 结果文件;
        若用户显式 --skip-checkm, 直接回退默认估计;
        否则调用 CheckM;
        纯蛋白输入或 CheckM 未安装/失败时回退默认估计.
        """
        has_checkm_file = self.config.checkm_results and Path(self.config.checkm_results).exists()
        if has_checkm_file:
            logger.info(f"  Using pre-computed CheckM results: {self.config.checkm_results}")
            results = _parse_checkm_results(self.config.checkm_results, genomes)
            # Sentinel so preprocess can report quality_source = precomputed
            results["__precomputed__"] = GenomeQuality(completeness=0.0, contamination=0.0, quality_score=0.0)
            return results

        if self.config.skip_checkm:
            logger.info("  --skip-checkm enabled; skipping CheckM and using default quality estimates.")
            return _no_checkm_fallback(genomes, sentinel="__skipped__")

        return _run_checkm(genomes, self.config)

    def _stratify_genomes(self, quality_results: Dict[str, GenomeQuality]) -> Dict[str, List[str]]:
        layers: Dict[str, List[str]] = {"high": [], "medium": [], "low": [], "contaminated": []}
        for gid, q in quality_results.items():
            if q.contamination >= 10:
                layers["contaminated"].append(gid)
            elif q.completeness > 90 and q.contamination < 5:
                layers["high"].append(gid)
            elif q.completeness >= 70 and q.contamination < 10:
                layers["medium"].append(gid)
            else:
                layers["low"].append(gid)
        return layers

    def _compute_adaptive_parameters(self, layers: Dict[str, List[str]]) -> AdaptiveParams:
        total = sum(len(v) for v in layers.values())
        high_ratio = len(layers["high"]) / total if total > 0 else 0
        if high_ratio > ADAPTIVE_HIGH_RATIO_THRESHOLD:
            return AdaptiveParams(
                min_hmm_score=ADAPTIVE_HIGH_MIN_HMM_SCORE,
                min_occupancy=ADAPTIVE_HIGH_MIN_OCCUPANCY,
                max_markers=ADAPTIVE_HIGH_MAX_MARKERS,
                missing_data_strategy="strict_gap", quality_weighted=True,
            )
        elif high_ratio > ADAPTIVE_MEDIUM_RATIO_THRESHOLD:
            return AdaptiveParams(
                min_hmm_score=ADAPTIVE_MEDIUM_MIN_HMM_SCORE,
                min_occupancy=ADAPTIVE_MEDIUM_MIN_OCCUPANCY,
                max_markers=ADAPTIVE_MEDIUM_MAX_MARKERS,
                missing_data_strategy="gap", quality_weighted=True,
            )
        else:
            return AdaptiveParams(
                min_hmm_score=ADAPTIVE_LOW_MIN_HMM_SCORE,
                min_occupancy=ADAPTIVE_LOW_MIN_OCCUPANCY,
                max_markers=ADAPTIVE_LOW_MAX_MARKERS,
                missing_data_strategy="gap", quality_weighted=True, allow_partial_hits=True,
            )


class SparseSupermatrixOptimizer:
    """Reserved interface for sparse supermatrix optimization on MAG-heavy data.

    Currently returns the input matrix unchanged; the main pipeline relies on
    the marker-selection strategies (especially SPARSE_OPTIMIZED) to handle
    incomplete genomes. Retained for future extensions.
    """
    def __init__(self, config: MAGConfig):
        self.config = config

    def optimize(self, occupancy_matrix, quality_data):
        return OptimizedMatrix(
            matrix=occupancy_matrix,
            selected_markers=list(occupancy_matrix.cogs),
            selected_species=list(occupancy_matrix.genomes),
        )


class MAGOptimizationModule:
    def __init__(self, config: MAGConfig):
        self.config = config
        self.preprocessor = MAGQualityPreprocessor(config)
        self.heterogeneity_handler = MAGHeterogeneityHandler(HeterogeneityConfig())

    def run(self, genomes: List[Genome]) -> PreprocessingResult:
        return self.preprocessor.preprocess(genomes)
