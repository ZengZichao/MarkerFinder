"""Global configuration dataclasses for MarkerFinder pipeline."""

import os
import shutil
import tempfile
import uuid
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

from markerfinder.models.marker import SelectionStrategy


# ---------------------------------------------------------------------------
# HGT 风险分桶 / 分级阈值 —— 单一事实源.
#
# 这些常量同时被以下位置引用, 避免 0.25 / 0.60 等 HGT 业务阈值散落为魔法数字:
# * ``HGTConfig.level_thresholds`` 的默认值 (config.py 内部);
# * 报告生成模块 report_generator 的 HGT 风险分桶边界 (risk_buckets);
# * CLI 入口 __main__.py 构造 HGTConfig 时的 level2_max 兜底值.
# 修改 HGT 分级边界时, 只需改这里一处, 全链路取值自动保持一致.
# ---------------------------------------------------------------------------
HGT_LEVEL1_MAX = 0.25  # 风险 < 该值 -> Level 1 (clean)
HGT_LEVEL2_MAX = 0.60  # 风险 < 该值 -> Level 2 (suspicious); 否则 Level 3 (excluded)


# Per-run scratch root: {system_tmp}/markerfinder-{8-hex}. Keeps the domtblout /
# Concat files of concurrent processes apart, and is cleaned up when the run
# Ends.
#
# Deliberately a function, not a module-level constant: the constant was
# Computed *and created* at import time, so merely importing markerfinder left a
# Directory behind (in the working tree, whenever the system temp dir is not
# Writable), and every configuration built in one process shared that single
# Directory. A default_factory gives each configuration its own, and creating it
# Is left to the code that is about to write there (Pipeline._prepare_context).
def _default_tmp_dir() -> str:
    return os.path.join(tempfile.gettempdir(),
                        f"markerfinder-{uuid.uuid4().hex[:8]}")


@dataclass
class BaseConfig:
    """Shared base for configuration dataclasses that expose a CPU count.

    Centralises the previously duplicated ``cpus: int = 1`` field so the thread
    count is configured in exactly one place.
    """
    cpus: int = 1


@dataclass
class SelectionConfig(BaseConfig):
    """标记基因选择配置"""
    marker_mode: str = "gtdb_tk"
    marker_hmm_dir: str = ""
    # 当未显式指定 marker_hmm_dir 时如何定位标记 HMM: auto=按约定路径在 db/ 下
    # 自动发现(ar53/bac120 按输入基因组域选取); cog=仅看 cog_hmm; gtdb=仅看 gtdb_markers;
    # None=不自动发现(必须显式 --marker-hmm-dir).
    marker_db_source: str = "auto"
    min_hmm_score: float = 20.0
    min_occupancy: float = 0.75
    # 用户通过 --min-occupancy / --min-marker-coverage 显式指定的占用率下限.
    # None 表示未指定(沿用 Phase 0 自适应分层阈值); 非 None 时覆盖自适应值,
    # 保证显式用户输入优先于自动推断.
    user_min_occupancy: Optional[float] = None
    max_markers: int = 60
    # 与 user_min_occupancy 同理: 用户通过 --max-markers(或显式 --marker-preset)
    # 给出的标记预算. None 表示未指定(沿用 Phase 0 自适应预算); 非 None 时覆盖
    # 自适应值 —— 否则 --max-markers 永远被 adaptive_params.max_markers 淹没,
    # 成为一个只被记录、不被执行的旋钮.
    user_max_markers: Optional[int] = None
    strategy: SelectionStrategy = SelectionStrategy.GREEDY
    gtdb_markers_dir: str = ""
    species_tree: str = ""
    tmp_dir: str = field(default_factory=_default_tmp_dir)
    quality_weighted: bool = True


@dataclass
class OrthologConfig(BaseConfig):
    """直系同源物选择配置"""
    bbh_evalue: float = 1e-10
    bbh_identity: float = 30.0
    bbh_coverage: float = 0.6
    graph_clustering: bool = True
    use_diamond: bool = True
    tmp_dir: str = field(default_factory=_default_tmp_dir)


@dataclass
class HGTConfig(BaseConfig):
    """HGT过滤配置

    自包含默认值: 启用系统发育(Phylogenetic) 步骤.
    系统发育步骤基于基因树-物种树不一致性(MAD 定根 + 单系比例).
    不依赖任何外部大型比对库(如 NCBI nr), 整套流程自包含运行.
    """
    enable_phylogenetic: bool = True

    tmp_dir: str = field(default_factory=_default_tmp_dir)

    phylogenetic_threshold: float = 0.5

    # Phylogenetic step monophyly screen (MAD rooting + scope-based monophyly
    # Proportion). Only used when a taxonomy table is supplied. The tree's
    # Taxonomic *scope* is inferred from its tips (deepest rank all tips share)
    # And, under the default ``auto``, the monophyly proportion is computed at
    # The rank ONE LEVEL BELOW the scope (domain->phylum,..., genus->species).
    # A marker whose gene tree disagrees with the taxonomy (low monophyly
    # Proportion) is flagged as HGT-prone. Naming a rank instead of ``auto``
    # Measures at that rank; the scope is then reported but never applied.
    monophyly_rank: str = "auto"
    monophyly_threshold: float = 0.5

    level_thresholds: Dict[str, float] = field(default_factory=lambda: {
        "level1_max": HGT_LEVEL1_MAX,
        "level2_max": HGT_LEVEL2_MAX,
    })

    # 远缘数据集上的宽容模式. 当输入基因组跨多个高阶分类单元(属/科/目),
    # 系统发育 distant 基因组间风险评分饱和,会把所有标记归为 等级 3 剔除.
    # 开启此模式后, level2_max 放宽到 ``level2_max_far`` 以保留足够标记建树.
    # 默认关闭: 此模式会悄悄放宽 HGT 放行尺度, 必须由用户在配置/CLI 显式开启
    # (``--hgt-adaptive-thresholds``), 避免在无感知情况下改变最终标记集.
    adaptive_far_thresholds: bool = False
    far_distance_genera_ratio: float = 0.70  # N_genera/n_genomes 触发远缘模式
    far_distance_min_orders: int = 2         # N_orders 触发远缘模式
    level2_max_far: float = 0.95              # 远缘模式下的 level2_max

    # 聚合权重: 多 HGT 信号(远缘 RF 距离、quartet 一致性)合成为单一风险分时
    # 的权重. 默认等权 {"rf": 0.5, "quartet": 0.5} —— 在 HGT 风险评分尚无
    # 大规模敏感性分析结论前, 等权是最保守且不偏向任一信号的中性选择.
    # 该字典被 ``combine_hgt_scores`` 读取, 用户可在配置/代码显式覆盖.
    hgt_score_weights: Dict[str, float] = field(
        default_factory=lambda: {"rf": 0.5, "quartet": 0.5}
    )

    # Screening criterion mode. Default 'risk'
    # NEVER auto-switches; consistency/hybrid additionally require the
    # Prerequisite gate. --min-informative-sites is consumed by 's
    # Informativeness screen (provisional, 0 = disabled).
    hgt_mode: str = "risk"
    consistency_stringency: int = 1
    min_informative_sites: int = 0


@dataclass
class HeterogeneityConfig:
    """Composition / GC heterogeneity diagnostics.

    Promoted from a "reserved interface" to a real configuration. The screen
    is DEFAULT-OFF: it produces parallel evidence columns only — the
    composition metrics never merge into overall_risk/overall_score.
    """
    min_allele_freq: float = 0.2
    consensus_threshold: float = 0.7
    minor_allele_threshold: float = 0.3
    identity_threshold: float = 0.99
    # Composition screen switch + thresholds.
    enable_composition_screen: bool = False
    composition_outlier_metric: str = "rcv"      # "rcv" (default)
    composition_warn_threshold: float = 0.15     # Provisional, pending


@dataclass
class MAGConfig(BaseConfig):
    """MAG优化配置"""
    checkm_results: Optional[str] = None
    skip_checkm: bool = False
    min_completeness: float = 50.0
    min_marker_coverage: float = 0.3
    quality_weighted: bool = True
    tmp_dir: str = field(default_factory=_default_tmp_dir)
    input_dir: str = ""
    max_markers: int = 150


@dataclass
class AlignerConfig(BaseConfig):
    """比对配置"""
    tmp_dir: str = field(default_factory=_default_tmp_dir)
    output_prefix: str = "markerfinder"


@dataclass
class PhylogeneticConfig(BaseConfig):
    """系统发育推断配置"""
    tmp_dir: str = field(default_factory=_default_tmp_dir)
    output_prefix: str = "markerfinder"
    ufboot_replicates: int = 1000
    fast_mode: bool = False
    use_fasttree: bool = True
    gene_tree_builder: Optional[str] = None  # "fasttree" or "iqtree"
    iqtree_timeout: int = 3600
    astral_timeout: int = 1800
    # Optional minimum average local-support (UFBOOT / local support) for a
    # Gene tree to be admitted into the coalescent species-tree inference.
    # ``None`` (default) disables the support gate so behaviour is unchanged:
    # Only the ``n_tips >= 4`` filter is applied. Set to a float (e.g. 0.0 or
    # 50.0) to also drop low-support gene trees that would degrade ASTRAL-III.
    min_gene_tree_support: Optional[float] = None
    coalescent_mode: str = "post-filter"
    # Recommendation band thresholds, promoted from the
    # In-method literals 0.1/0.3 so they are configurable and reportable.
    recommend_rf_low: float = 0.1
    recommend_rf_moderate: float = 0.3


@dataclass
class ReportConfig:
    """报告配置"""
    output_dir: str = "./output"
    output_prefix: str = "markerfinder"
    template_dir: str = "markerfinder/templates"
    report_format: str = "html"
    save_intermediates: bool = False
    # HGT 风险分桶边界, 与 HGTConfig.level_thresholds 同源(由 CLI --hgt-threshold
    # 构造), 供报告的风险分布图使用; 默认值取全局常量, 避免报告与实际分级漂移.
    level1_max: float = HGT_LEVEL1_MAX
    level2_max: float = HGT_LEVEL2_MAX
    # Pipeline_summary 在覆盖率低于此值时于顶部显著警示.
    # CLI 对应 --require-evidence-coverage ( 接线); 默认仅告警不阻断.
    require_evidence_coverage: float = 0.5
    # Assertion layer switches. strict=on by default; the
    # Bypass flag is a deliberate, reported escape hatch.
    strict_assertions: bool = True
    allow_assertion_failure: bool = False
    # Threshold scan switch (writes
    # Phase5_reports/<prefix>.threshold_scan.tsv).
    hgt_scan: bool = False
    # Minimum acceptable cross-band Jaccard of the selected
    # Marker sets. Previously only a CLI flag with no consumer, so
    # ``--scan-stability-min`` silently did nothing.
    scan_stability_min: float = 0.6
    # Composition / GC diagnostics as PARALLEL evidence
    # Columns. Default off. These metrics are never merged into
    # Overall_risk / overall_score (locked by test).
    composition_screen: bool = False
    # The report has to say which criterion adjudicated. The generator
    # Cannot infer it from the numbers (a risk score looks the same in both
    # Modes), so the mode is passed through explicitly.
    hgt_mode: str = "risk"
    # Optional "marker_id<TAB>functional_category" TSV. None/missing file
    # => the column renders NA with a report note (never a blank lie).
    cog_category_map: Optional[str] = None
    # Path to the taxonomy must-pass baseline. None = the
    # Gate is off, and the run says so explicitly instead of implying a check
    # Happened.
    taxonomy_mustpass: Optional[str] = None


@dataclass
class TaxonomyConfig:
    """分类学相关配置"""
    taxonomy_table: Optional[str] = None
    taxonomy_format: str = "table"
    taxonomy_source_priority: str = "table"
    taxonomy_delimiter_mode: str = "reverse"
    table_sep: Optional[str] = None
    ignore_malformed: bool = False
    taxonomy_levels: Optional[str] = None


@dataclass
class PipelineConfig(BaseConfig):
    """流水线主配置"""
    input_dir: str = ""
    output_dir: str = "./output"
    output_prefix: str = "markerfinder"

    tmp_dir: str = field(default_factory=_default_tmp_dir)
    # Tmp_dir 是否为本进程自动生成的每趟独立目录. True(默认, 未传 --tmp-dir)时
    # 流程结束自动清理; 用户显式指定 --tmp-dir 时不自动清理(与 --tmp-dir 帮助
    # 文本的承诺一致, 避免误删用户目录).
    tmp_dir_auto: bool = True

    keep_tmp: bool = False  # 流程结束时保留 tmp 目录(用于调试)

    save_intermediates: bool = False  # 将有用中间产物复制到 output/Phase4_intermediate/

    mode: str = "standard"

    mag_config: MAGConfig = field(default_factory=MAGConfig)
    selection_config: SelectionConfig = field(default_factory=SelectionConfig)
    hgt_config: HGTConfig = field(default_factory=HGTConfig)
    ortholog_config: OrthologConfig = field(default_factory=OrthologConfig)
    aligner_config: AlignerConfig = field(default_factory=AlignerConfig)
    phylo_config: PhylogeneticConfig = field(default_factory=PhylogeneticConfig)
    # Composition / GC diagnostics thresholds.
    # ``pipeline.py`` has read ``getattr(self.config, "heterogeneity_config",...)``
    # Since that feature landed, but no such field existed, so the chain always
    # Resolved to None and ``composition_warn_threshold`` could never be anything
    # But the inline fallback -- a consumer wired to an object nobody provides.
    # Declaring it here makes the read real; ``to_dict`` then carries the value
    # Into run_config.json, so a reproduced run keeps the same threshold.
    heterogeneity_config: HeterogeneityConfig = field(
        default_factory=HeterogeneityConfig
    )
    report_config: ReportConfig = field(default_factory=ReportConfig)
    taxonomy_config: TaxonomyConfig = field(default_factory=TaxonomyConfig)

    db_dir: str = "./db"
    databases: Dict[str, str] = field(default_factory=dict)
    verify_db: bool = True

    force: bool = False
    no_clobber: bool = False

    sequences_path: Optional[str] = None
    mol_type: Optional[str] = None
    skip_length_check: bool = False
    strip_annotations: bool = False

    # Per-run locations and resource limits that must mean the same thing in
    # Every section. ``PipelineConfig(output_dir="run1")`` used to leave
    # ``report_config.output_dir`` at its own default ``"./output"``, so the
    # Reports, trees and partition file landed in the current directory while
    # The run's own snapshot landed in ``run1`` — two answers to "where does
    # This run write". The parent now propagates them at construction.
    _PROPAGATED_FIELDS = ("input_dir", "output_dir", "output_prefix", "tmp_dir",
                          "cpus")

    def __post_init__(self) -> None:
        import dataclasses

        for field in dataclasses.fields(self):
            section = getattr(self, field.name, None)
            if not dataclasses.is_dataclass(section) or isinstance(section, type):
                continue
            for name in self._PROPAGATED_FIELDS:
                if not dataclasses.fields(section) or not any(
                        f.name == name for f in dataclasses.fields(section)):
                    continue
                value = getattr(self, name, None)
                if value not in (None, ""):
                    setattr(section, name, value)

    def to_dict(self) -> Dict[str, Any]:
        """序列化为字典，用于结果快照"""
        from dataclasses import asdict
        return asdict(self)
