"""Regression test for review item N-B13: ``--marker-preset`` wrote the strategy
to a non-existent field.

``markerfinder/cli/config_build.py`` used to apply an explicit preset with
``selection_config.selection_strategy = preset.selection_strategy``, but the real
field of ``SelectionConfig`` (``markerfinder/config.py``) is ``strategy`` — and
``markerfinder/modules/marker_selection.py`` reads ``self.config.strategy``.
Because ``SelectionConfig`` is a plain (non-slots, non-frozen) dataclass, the
bad assignment did not raise: it hung a *shadow* attribute off the instance, so
the builder logged ``strategy=info_max`` while selection still ran ``GREEDY``.

These tests pin the fixed behaviour through the real parser and the real config
builder, and assert the two properties whose absence let the bug hide:

* the preset's strategy lands on the *declared* field, and no shadow attribute
  is created;
* the strategy the builder logs is the strategy it applied.
"""

import dataclasses
import logging

import pytest

from markerfinder.cli.config_build import _build_pipeline_config
from markerfinder.cli.parser import _build_parser
from markerfinder.config import SelectionConfig
from markerfinder.models.marker import (
    RESOLUTION_PRESETS,
    ResolutionPreset,
    SelectionStrategy,
)

# preset name -> the strategy that preset declares.
PRESET_STRATEGY = {
    "conservative": SelectionStrategy.INFO_MAX,
    "expanded": SelectionStrategy.RATE_BALANCED,
    "mag_adaptive": SelectionStrategy.SPARSE_OPTIMIZED,
    "standard": SelectionStrategy.GREEDY,
}

# --mode -> the strategy the mode alone selects.
MODE_STRATEGY = {
    "conservative": SelectionStrategy.INFO_MAX,
    "standard": SelectionStrategy.GREEDY,
    "expanded": SelectionStrategy.RATE_BALANCED,
    "mag_adaptive": SelectionStrategy.SPARSE_OPTIMIZED,
}


def _build(*extra):
    args = _build_parser().parse_args(["--output", "out", *extra])
    return _build_pipeline_config(args)


def _declared_fields():
    return {f.name for f in dataclasses.fields(SelectionConfig)}


def test_selection_config_field_is_named_strategy_not_selection_strategy():
    """The single source of the bug: the field name is ``strategy``."""
    assert "strategy" in _declared_fields()
    assert "selection_strategy" not in _declared_fields()


def test_selection_config_is_a_plain_dataclass():
    """Why the typo was silent rather than an AttributeError.

    Documented here so that making the class frozen/slots-frozen — which would
    turn every such typo into a hard failure — is a deliberate, visible change.
    """
    params = SelectionConfig.__dataclass_params__
    assert not params.frozen
    assert "__slots__" not in SelectionConfig.__dict__


@pytest.mark.parametrize("preset", sorted(PRESET_STRATEGY))
def test_preset_applies_its_strategy_to_the_declared_field(preset):
    expected = PRESET_STRATEGY[preset]
    cfg = _build("--marker-preset", preset)

    # The field the selection module actually reads.
    assert cfg.selection_config.strategy is expected
    # No shadow attribute left behind by the old typo.
    assert "selection_strategy" not in vars(cfg.selection_config)
    # The preset's declared value and the applied value come from one source.
    declared = RESOLUTION_PRESETS[ResolutionPreset(preset)].selection_strategy
    assert declared is expected


def test_standard_preset_overrides_a_non_matching_mode():
    """``standard`` declares GREEDY, which is also the default — so pair it with
    a mode that selected something else; an unapplied preset would still read
    GREEDY only if the preset were ignored, which this rules out."""
    assert _build("--mode", "expanded").selection_config.strategy is (
        SelectionStrategy.RATE_BALANCED
    )
    cfg = _build("--mode", "expanded", "--marker-preset", "standard")
    assert cfg.selection_config.strategy is SelectionStrategy.GREEDY


@pytest.mark.parametrize("preset", sorted(PRESET_STRATEGY))
def test_logged_strategy_matches_applied_strategy(preset, caplog):
    """The log line may never advertise a strategy the run does not use."""
    with caplog.at_level(logging.INFO, logger="markerfinder.cli.config_build"):
        cfg = _build("--marker-preset", preset)
    applied = cfg.selection_config.strategy.value
    assert f"Applied --marker-preset {preset}" in caplog.text
    assert f"strategy={applied}" in caplog.text
    assert f"strategy={RESOLUTION_PRESETS[ResolutionPreset(preset)].selection_strategy.value}" in caplog.text


@pytest.mark.parametrize("mode", sorted(MODE_STRATEGY))
def test_preset_none_leaves_the_configured_strategy_untouched(mode):
    """``none`` (the default) must not touch the strategy the mode selected."""
    without_flag = _build("--mode", mode)
    with_none = _build("--mode", mode, "--marker-preset", "none")

    assert without_flag.selection_config.strategy is MODE_STRATEGY[mode]
    assert with_none.selection_config.strategy is MODE_STRATEGY[mode]
    # 'none' leaves every other knob alone as well.
    assert with_none.selection_config.min_occupancy == (
        without_flag.selection_config.min_occupancy
    )
    assert with_none.selection_config.max_markers == (
        without_flag.selection_config.max_markers
    )
    assert "selection_strategy" not in vars(with_none.selection_config)


def test_preset_still_overrides_occupancy_and_budget():
    """The other two members of the triple keep working after the fix."""
    cfg = _build("--marker-preset", "conservative")
    assert cfg.selection_config.min_occupancy == 0.90
    assert cfg.selection_config.max_markers == 30
    # An explicit preset also wins over the Phase-0 adaptive budget.
    assert cfg.selection_config.user_max_markers == 30
