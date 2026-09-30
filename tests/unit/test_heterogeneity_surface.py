"""The heterogeneity configuration must actually reach the product.

Two defects lived in one line. ``pipeline.py`` read
``getattr(self.config, "heterogeneity_config", None) or getattr(self.config,
"heterogeneity", None)`` -- but ``PipelineConfig`` declared neither name, so the
chain ALWAYS resolved to ``None`` and the outlier threshold used was the inline
fallback, whatever anyone set. A consumer wired to an object nobody provides is
the same shape as the baseline's ("存为属性后从未被读取"), just mirrored: read
but never written.

And the shipped docs had to admit that
``HeterogeneityConfig.enable_composition_screen`` did nothing while
``hgt_steps: composition`` did -- two spellings of one intent, one of them inert.
Both now work; the defaults still keep the shipped products identical ( /
), which the control below pins.
"""

from __future__ import annotations

import dataclasses
import logging

import pytest

from markerfinder.config import HeterogeneityConfig, PipelineConfig
from markerfinder.modules.composition import write_composition_tsv
from markerfinder.pipeline import (
    MarkerFinderPipeline,
    _composition_screen_enabled,
    _warn_about_outlier_metric,
)


# ── the resolved-consumer regression ──────────────────────────────────────

def test_control_pipeline_config_now_carries_the_section():
    config = PipelineConfig()
    assert isinstance(config.heterogeneity_config, HeterogeneityConfig)
    # The exact expression pipeline.py uses; before the field existed this
    # Resolved to None and the threshold below was unreachable.
    resolved = (
        getattr(config, "heterogeneity_config", None)
        or getattr(config, "heterogeneity", None)
    )
    assert resolved is config.heterogeneity_config


def test_defaults_are_the_shipped_ones():
    """Turning the dead read into a live one must not move any default."""
    defaults = HeterogeneityConfig()
    assert defaults.composition_warn_threshold == pytest.approx(0.15)
    assert defaults.enable_composition_screen is False
    assert defaults.composition_outlier_metric == "rcv"


def test_threshold_survives_the_reproducibility_snapshot():
    """A run's threshold has to be recoverable from its own products."""
    config = PipelineConfig()
    config.heterogeneity_config.composition_warn_threshold = 0.42
    snapshot = config.to_dict()
    assert snapshot["heterogeneity_config"]["composition_warn_threshold"] == 0.42


# ── the threshold is load-bearing, not decorative ─────────────────────────

def _flags(tmp_path, rcv_value: float, threshold):
    path = write_composition_tsv(
        [{"marker_id": "m1", "rcv": rcv_value, "gc_bias": None, "n_sequences": 12}],
        str(tmp_path),
        "mf",
        warn_threshold=threshold,
    )
    lines = path.read_text(encoding="utf-8").strip().split("\n")
    assert lines[0] == "marker_id\trcv\tgc_bias\tn_sequences\toutlier_flag"
    return lines[1].split("\t")[-1]


def test_control_the_same_rcv_flags_differently_at_the_two_thresholds(tmp_path):
    assert _flags(tmp_path, 0.30, 0.15) == "warn"
    assert _flags(tmp_path, 0.30, 0.90) == "NA"


def test_control_missing_rcv_is_na_and_not_a_zero_flag(tmp_path):
    assert _flags(tmp_path, 0.0, 0.0) == "warn", "0 >= 0 is a real comparison"
    path = write_composition_tsv(
        [{"marker_id": "m1", "rcv": None, "gc_bias": None, "n_sequences": 1}],
        str(tmp_path / "none"), "mf", warn_threshold=0.15,
    )
    row = path.read_text(encoding="utf-8").strip().split("\n")[1].split("\t")
    assert row[1] == "NA" and row[2] == "NA" and row[4] == "NA", row


# ── either spelling enables; the metric is never quietly ignored ──────────

@pytest.mark.parametrize(
    "screen_flag,heterogeneity,expected",
    [
        (False, HeterogeneityConfig(), False),
        (True, HeterogeneityConfig(), True),
        (False, HeterogeneityConfig(enable_composition_screen=True), True),
        (True, HeterogeneityConfig(enable_composition_screen=True), True),
        (False, None, False),
        (True, None, True),
    ],
)
def test_composition_screen_truth_table(screen_flag, heterogeneity, expected):
    assert (
        _composition_screen_enabled(screen_flag, heterogeneity) is expected
    ), (screen_flag, heterogeneity)


def test_unimplemented_outlier_metric_is_named_when_the_screen_runs(caplog):
    config = HeterogeneityConfig(composition_outlier_metric="chi2")
    with caplog.at_level(logging.WARNING, logger="markerfinder.pipeline"):
        _warn_about_outlier_metric(config, enabled=True)
    assert "chi2" in caplog.text
    assert "only metric is 'rcv'" in caplog.text


def test_control_no_warning_for_the_implemented_metric(caplog):
    with caplog.at_level(logging.WARNING, logger="markerfinder.pipeline"):
        _warn_about_outlier_metric(HeterogeneityConfig(), enabled=True)
        _warn_about_outlier_metric(
            HeterogeneityConfig(composition_outlier_metric="chi2"), enabled=False
        )
    assert caplog.text == "", caplog.text


def test_pipeline_still_constructs_with_the_new_section(tmp_path):
    """Smoke: the added field must not break construction or the dict round trip."""
    config = PipelineConfig(input_dir=str(tmp_path), output_dir=str(tmp_path / "out"))
    MarkerFinderPipeline(config)
    restored = dataclasses.fields(PipelineConfig)
    assert any(f.name == "heterogeneity_config" for f in restored)
