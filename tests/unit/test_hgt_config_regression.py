"""Code-review regression tests for ``config.py`` (#2 + C33) and ``utils/hgt_utils.py`` (C33).

Locks in:
  - #2: ``HGTConfig.adaptive_far_thresholds`` now defaults to ``False``
        (was ``True``); it must be explicitly enabled.
  - C33: ``HGTConfig.hgt_score_weights`` exists and defaults to
         ``{"rf": 0.5, "quartet": 0.5}`` (configurable).
  - C33: ``combine_hgt_scores`` reads the configured weights so the combined
         HGT risk score weights RF / quartet consistency accordingly; default
         weighting is symmetric.
"""

import pytest

from markerfinder.config import HGTConfig
from markerfinder.utils.hgt_utils import combine_hgt_scores


class TestHGTConfigDefaults:
    def test_adaptive_far_thresholds_default_false(self):
        assert HGTConfig().adaptive_far_thresholds is False

    def test_default_hgt_score_weights(self):
        assert HGTConfig().hgt_score_weights == {"rf": 0.5, "quartet": 0.5}

    def test_custom_hgt_score_weights(self):
        cfg = HGTConfig(hgt_score_weights={"rf": 0.8, "quartet": 0.2})
        assert cfg.hgt_score_weights == {"rf": 0.8, "quartet": 0.2}


class TestCombineHgtScores:
    def test_default_symmetric_weighting(self):
        scores = {"rf": 0.4, "quartet": 0.6}
        out, n_used = combine_hgt_scores(scores, HGTConfig().hgt_score_weights)
        assert abs(out - 0.5) < 1e-9
        assert n_used == 2

    def test_custom_weights_rf_dominant(self):
        cfg = HGTConfig(hgt_score_weights={"rf": 0.8, "quartet": 0.2})
        scores = {"rf": 1.0, "quartet": 0.0}
        out, n_used = combine_hgt_scores(scores, cfg.hgt_score_weights)
        assert abs(out - 0.8) < 1e-9
        assert n_used == 2

    def test_custom_weights_quartet_dominant(self):
        cfg = HGTConfig(hgt_score_weights={"rf": 0.8, "quartet": 0.2})
        scores = {"rf": 0.0, "quartet": 1.0}
        out, n_used = combine_hgt_scores(scores, cfg.hgt_score_weights)
        assert abs(out - 0.2) < 1e-9
        assert n_used == 2

    def test_missing_score_key_ignored(self):
        # A signal absent from ``scores`` must not dilute the present one.
        cfg = HGTConfig(hgt_score_weights={"rf": 0.5, "quartet": 0.5})
        scores = {"rf": 1.0}
        out, n_used = combine_hgt_scores(scores, cfg.hgt_score_weights)
        assert abs(out - 1.0) < 1e-9
        assert n_used == 1  # N_used is reported, caller gates on <2
