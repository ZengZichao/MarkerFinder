""" (deprecation notice that actually fires) and (dead enums).

 asked that ``--hgt-threshold`` stay a backward-compatible alias of
``--hgt-threshold-l1l2`` and *print a deprecation notice when used*. The code
carried a comment saying it did so, but nothing emitted anything, and the flag
defaulted to ``0.25`` — which also made "did the user pass it?" unanswerable.

 asked for ``ConflictType.ILS_SIGNAL``/``METHOD_BIAS`` to be deleted
(declared, never assigned anywhere). They were still in the enum.
"""

from __future__ import annotations

import logging

import pytest

from markerfinder.models.pipeline_types import ConflictType


def _config_from(args):
    from markerfinder.cli.config_build import _build_pipeline_config
    return _build_pipeline_config(args)


class TestDeprecationNotice:
    def test_unset_alias_is_now_observable(self):
        """The alias must default to None, or 'was it used?' cannot be known."""
        from markerfinder.cli.parser import _build_parser

        args = _build_parser().parse_args(["--output", "out"])
        assert args.hgt_threshold is None

    def test_no_notice_when_the_alias_is_not_used(self, caplog):
        from markerfinder.cli.parser import _build_parser

        args = _build_parser().parse_args(["--output", "out"])
        with caplog.at_level(logging.WARNING):
            cfg = _config_from(args)
        assert cfg.hgt_config.level_thresholds["level1_max"] == 0.25
        assert not [r for r in caplog.records
                    if "--hgt-threshold is deprecated" in r.getMessage()], (
            "a deprecation notice fired even though the user never passed "
            "--hgt-threshold"
        )

    def test_using_the_alias_prints_exactly_one_notice(self, caplog):
        from markerfinder.cli.parser import _build_parser

        args = _build_parser().parse_args(
            ["--output", "out", "--hgt-threshold", "0.30"]
        )
        with caplog.at_level(logging.WARNING):
            cfg = _config_from(args)
        fired = [
            r.getMessage() for r in caplog.records
            if "--hgt-threshold is deprecated" in r.getMessage()
        ]
        assert len(fired) == 1, fired
        assert "--hgt-threshold-l1l2" in fired[0]
        # Behaviour unchanged: the alias still moves only band 1/2.
        assert cfg.hgt_config.level_thresholds["level1_max"] == 0.30
        assert cfg.hgt_config.level_thresholds["level2_max"] == 0.60

    def test_explicit_new_flag_still_wins_but_alias_is_reported(self, caplog):
        from markerfinder.cli.parser import _build_parser

        args = _build_parser().parse_args(
            ["--output", "out", "--hgt-threshold", "0.30",
             "--hgt-threshold-l1l2", "0.45"]
        )
        with caplog.at_level(logging.WARNING):
            cfg = _config_from(args)
        assert cfg.hgt_config.level_thresholds["level1_max"] == 0.45
        assert any("--hgt-threshold is deprecated" in r.getMessage()
                   for r in caplog.records)

    def test_range_validation_tolerates_the_unset_alias(self):
        """A None default must not crash the 0-1 range check."""
        import inspect

        from markerfinder.cli import validation as cli_validation

        source = inspect.getsource(cli_validation)
        assert "args.hgt_threshold is not None" in source, (
            "the range check would raise TypeError on an unset --hgt-threshold"
        )


class TestDeadConflictTypesRemoved:
    def test_only_the_assigned_categories_remain(self):
        names = {member.name for member in ConflictType}
        assert names == {"HGT_SIGNAL", "AMBIGUOUS"}, names

    @pytest.mark.parametrize("gone", ["ILS_SIGNAL", "METHOD_BIAS"])
    def test_declared_but_never_assigned_categories_are_gone(self, gone):
        assert not hasattr(ConflictType, gone), (
            f"{gone} is back: the deprecated conflict categories must stay deleted "
            f"that nothing ever assigns"
        )

    def test_default_value_still_usable(self):
        assert ConflictType.AMBIGUOUS.value == "ambiguous"
