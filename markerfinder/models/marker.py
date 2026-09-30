from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


class MarkerLevel(Enum):
    LEVEL_1 = "level_1"
    LEVEL_2 = "level_2"
    LEVEL_3 = "level_3"
    # 外部工具缺失/失败、或系统发育步骤无法运行(无参考物种树/无 taxonomy 表)
    # 时, 标记的 HGT 状态**未知**而非被误判为 Level 3(高 HGT 风险并排除).
    # 用 UNKNOWN 显式标记, 避免"假装已排除".
    UNKNOWN = "unknown"


# Grading / normalization magic numbers promoted to module-level constants so
# They are defined in exactly one place (see MarkerQualityScore.overall_score
# And.level).
LEVEL_1_MIN_SCORE = 0.8
LEVEL_2_MIN_SCORE = 0.5
# Length-ratio coherence above this value saturates; below it floors out.
LENGTH_RATIO_SATURATION = 1.5
LENGTH_RATIO_LOW_FLOOR = 0.5
# Bitscore scale on which hmm_score is normalized (scores near/above this
# Saturate the hmm component of the quality score).
HMM_SCORE_NORMALIZER = 100.0


class SelectionStrategy(Enum):
    GREEDY = "greedy"
    INFO_MAX = "info_max"
    RATE_BALANCED = "rate_balanced"
    SPARSE_OPTIMIZED = "sparse_optimized"


@dataclass
class MarkerQualityScore:
    marker_id: str = ""
    hmm_score: float = 0.0
    occupancy: float = 0.0
    length_ratio: float = 0.0
    copy_number_cv: float = 0.0
    # ``None`` means "not measured" (no gene tree yet);
    # Phase 3 back-fills it from the cached gene trees.
    phylogenetic_informativeness: Optional[float] = None
    hgt_risk_score: Optional[float] = None
    functional_category: str = ""
    # Sequence-informativeness axis. Distinct from
    # ``phylogenetic_informativeness`` (the support proxy back-filled from
    # Gene trees) — the two quantities are physically different and must
    # Never share a name.
    pis: Optional[int] = None
    effective_columns: Optional[int] = None

    @property
    def overall_score(self) -> Optional[float]:
        """Five-factor weighted ordering score.

        this score is an ordering key only (never a filter).
        When a component is unmeasured (``None``) the composite is ``None``
        — substituting 0 would silently read as the worst possible value.
        """
        weights = {
            "hmm": 0.25,
            "occupancy": 0.25,
            "length": 0.10,
            "informativeness": 0.20,
            "hgt_cleanliness": 0.20,
        }
        if self.phylogenetic_informativeness is None or self.hgt_risk_score is None:
            return None
        hmm_norm = min(self.hmm_score / HMM_SCORE_NORMALIZER, 1.0)
        occ_norm = self.occupancy
        len_norm = min(self.length_ratio, 1.0) if self.length_ratio <= LENGTH_RATIO_SATURATION else LENGTH_RATIO_LOW_FLOOR
        pi_norm = self.phylogenetic_informativeness
        hgt_norm = 1.0 - self.hgt_risk_score
        score = (
            weights["hmm"] * hmm_norm
            + weights["occupancy"] * occ_norm
            + weights["length"] * len_norm
            + weights["informativeness"] * pi_norm
            + weights["hgt_cleanliness"] * hgt_norm
        )
        return min(max(score, 0.0), 1.0)

    @property
    def level(self) -> MarkerLevel:
        score = self.overall_score
        if score is None:
            # Unmeasured components → no quality level can be asserted.
            return MarkerLevel.UNKNOWN
        if score >= LEVEL_1_MIN_SCORE:
            return MarkerLevel.LEVEL_1
        elif score >= LEVEL_2_MIN_SCORE:
            return MarkerLevel.LEVEL_2
        else:
            return MarkerLevel.LEVEL_3


@dataclass
class MarkerGene:
    id: str = ""
    sequences: Dict = field(default_factory=dict)
    alignment: Optional[object] = None
    trimmed_alignment: Optional[object] = None
    gene_tree: Optional[object] = None
    quality_score: Optional[MarkerQualityScore] = None
    hgt_evaluation: Optional[object] = None


@dataclass
class SelectedMarkerSet:
    markers: List[str] = field(default_factory=list)
    occupancy_scores: Dict[str, float] = field(default_factory=dict)
    quality_scores: Dict[str, MarkerQualityScore] = field(default_factory=dict)
    strategy: SelectionStrategy = SelectionStrategy.GREEDY
    total_candidates: int = 0
    mean_occupancy: float = 0.0


class ResolutionPreset(Enum):
    CONSERVATIVE = "conservative"
    STANDARD = "standard"
    EXPANDED = "expanded"
    MAG_ADAPTIVE = "mag_adaptive"


@dataclass
class ResolutionConfig:
    min_occupancy: float = 0.75
    max_markers: int = 60
    allow_multi_copy: bool = False
    quality_weighted: bool = True
    selection_strategy: SelectionStrategy = SelectionStrategy.GREEDY


# 各分辨率档位的占用率下限 (``min_occupancy``) 与 ``SelectionConfig.min_occupancy``
# 是同一语义的阈值(标记集允许的最低占用率)。两者应保持一致: 默认档位 STANDARD
# 的 0.75 与 ``SelectionConfig`` 默认值对齐; 更保守的 CONSERVATIVE 抬高到 0.90,
# 更宽松的 EXPANDED/MAG_ADAPTIVE 降到 0.50/0.40。如需修改全局占用率默认值, 记得
# 同步更新 ``SelectionConfig.min_occupancy`` 与本处的对应档位, 避免两处语义漂移。
# 注意: 此处方引用 ``SelectionConfig`` 会造成循环导入(marker.py <-> config.py),
# 故此处仅以注释登记其关系, 不实际 import。
RESOLUTION_PRESETS: Dict[ResolutionPreset, ResolutionConfig] = {
    ResolutionPreset.CONSERVATIVE: ResolutionConfig(
        min_occupancy=0.90, max_markers=30,
        allow_multi_copy=False, quality_weighted=True,
        selection_strategy=SelectionStrategy.INFO_MAX,
    ),
    ResolutionPreset.STANDARD: ResolutionConfig(
        min_occupancy=0.75, max_markers=60,
        allow_multi_copy=False, quality_weighted=True,
        selection_strategy=SelectionStrategy.GREEDY,
    ),
    ResolutionPreset.EXPANDED: ResolutionConfig(
        min_occupancy=0.50, max_markers=120,
        allow_multi_copy=True, quality_weighted=False,
        selection_strategy=SelectionStrategy.RATE_BALANCED,
    ),
    ResolutionPreset.MAG_ADAPTIVE: ResolutionConfig(
        min_occupancy=0.40, max_markers=150,
        allow_multi_copy=True, quality_weighted=True,
        selection_strategy=SelectionStrategy.SPARSE_OPTIMIZED,
    ),
}
