"""Metrics self-check — known answers, plus the must-fail
control: shuffling a label MUST change P/R (a metric that cannot fail is
not a metric)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "benchmark"))

from metrics import (  # Noqa: E402
    block_bootstrap_ci,
    confusion,
    evidence_coverage,
    precision_recall_f1,
)

POSITIVES = ["G1", "G2", "G3", "G4"]       # Planted chimera carriers
PREDICTED = ["G1", "G2", "G5", "G6"]       # What the screen flagged


class TestMetricsSelfCheck:
    def test_confusion_counts(self):
        assert confusion(POSITIVES, PREDICTED) == (2, 2, 2)

    def test_precision_recall_f1(self):
        result = precision_recall_f1(POSITIVES, PREDICTED)
        assert result.precision == pytest.approx(0.5)
        assert result.recall == pytest.approx(0.5)
        assert result.f1 == pytest.approx(0.5)

    def test_perfect_prediction(self):
        result = precision_recall_f1(POSITIVES, POSITIVES)
        assert result.precision == result.recall == result.f1 == 1.0

    def test_empty_prediction_scores_zero(self):
        result = precision_recall_f1(POSITIVES, [])
        assert result.precision == 0.0 and result.recall == 0.0

    def test_shuffled_labels_change_pr(self):
        """Must-fail control: corrupting labels MUST move P/R."""
        original = precision_recall_f1(POSITIVES, PREDICTED)
        # Deliberately drop a true positive and add a false negative
        corrupted = precision_recall_f1(POSITIVES, ["G1", "G2", "G3"])
        assert corrupted.precision != original.precision
        assert corrupted.recall != original.recall

    def test_block_bootstrap_ci_deterministic_and_sane(self):
        low1, high1 = block_bootstrap_ci(POSITIVES, PREDICTED, seed=7)
        low2, high2 = block_bootstrap_ci(POSITIVES, PREDICTED, seed=7)
        assert (low1, high1) == (low2, high2)  # Determinism
        assert 0.0 <= low1 <= high1 <= 1.0

    def test_evidence_coverage(self):
        assert evidence_coverage(8, 10) == pytest.approx(0.8)
        assert evidence_coverage(0, 0) == 0.0  # Never NaN (A-11)
