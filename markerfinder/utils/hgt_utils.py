"""HGT detection utility functions."""

import logging
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)


def combine_hgt_scores(
    scores: Dict[str, float],
    weights: Dict[str, float],
) -> Tuple[Optional[float], int]:
    """整合多个 HGT 指标为单一风险分（加权平均），返回 ``(risk, n_used)``。

    聚合方式：``weighted_sum = Σ weights[k] * scores[k]`` 仅对 ``scores`` 中
    存在的键累加，再除以这些键的权重之和（即缺失信号不参与归一化，避免稀释
    已有信号）。

    ：返回值附带实际参与合成的信号个数 ``n_used``。没有证据时
    返回 ``(None, 0)`` 而非 ``0.0`` —— 让"没有证据"与"风险为 0"在类型上
    不可混同。调用方在 ``n_used < 2`` 时不得自动分档（需置 UNKNOWN 并记
    ``single_signal_only``），除非显式允许单信号。

    权重来源：默认等权 ``{"rf": 0.5, "quartet": 0.5}``（见
    ``HGTConfig.hgt_score_weights``），在 HGT 风险评分尚无大规模敏感性分析
    结论前，等权是最保守且不对任一信号（远缘 RF / quartet 一致性）偏置的中性
    选择。调用方（如 ``PhylogeneticHGTDetector``）应从配置读取权重传入，而非
    硬编码。
    """
    total_weight = 0.0
    weighted_sum = 0.0
    n_used = 0
    for key, weight in weights.items():
        if key in scores:
            weighted_sum += scores[key] * weight
            total_weight += weight
            n_used += 1
    if total_weight > 0:
        return weighted_sum / total_weight, n_used
    return None, 0
