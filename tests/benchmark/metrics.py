"""Benchmark metrics — precision / recall / F1 + block
bootstrap confidence intervals (block = genome). Pure stdlib."""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple


@dataclass
class PRResult:
    precision: float
    recall: float
    f1: float
    tp: int
    fp: int
    fn: int
    ci: Tuple[float, float] = (0.0, 1.0)


def confusion(positives: Sequence[str], predicted: Sequence[str]) -> Tuple[int, int, int]:
    """TP/FP/FN over genome ids (block units)."""
    pos_set = set(positives)
    pred_set = set(predicted)
    tp = len(pos_set & pred_set)
    fp = len(pred_set - pos_set)
    fn = len(pos_set - pred_set)
    return tp, fp, fn


def precision_recall_f1(positives: Sequence[str], predicted: Sequence[str]) -> PRResult:
    tp, fp, fn = confusion(positives, predicted)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall)
        else 0.0
    )
    return PRResult(precision, recall, f1, tp, fp, fn)


def block_bootstrap_ci(
    positives: Sequence[str],
    predicted: Sequence[str],
    *,
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
) -> Tuple[float, float]:
    """Percentile CI for recall via genome-level block bootstrap.

    Resampling units are genomes (the independent blocks), never markers.
    Fixed seed -> deterministic.
    """
    blocks: Dict[str, int] = {}
    for genome in list(positives) + list(predicted):
        blocks.setdefault(genome, len(blocks))
    names = sorted(blocks)
    if not names:
        return (0.0, 1.0)
    rng = random.Random(seed)
    recalls: List[float] = []
    for _ in range(n_boot):
        sample = [names[rng.randrange(len(names))] for _ in names]
        sample_set = set(sample)
        pos = [g for g in positives if g in sample_set]
        pred = [g for g in predicted if g in sample_set]
        recalls.append(precision_recall_f1(pos, pred).recall)
    recalls.sort()
    low = recalls[int((alpha / 2) * len(recalls))]
    high = recalls[min(len(recalls) - 1, int((1 - alpha / 2) * len(recalls)))]
    return (low, high)


def evidence_coverage(n_measured: int, n_markers: int) -> float:
    """Report coverage alongside every P/R number."""
    return n_measured / n_markers if n_markers else 0.0
