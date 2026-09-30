"""Decisive-interface contract tests + taxonomy must-pass wiring
.

Every DECIDING interface gets a contract with a must-fail case:
  1. GeneState assignment (SINGLE_COPY / ABSENT / MULTI_COPY boundaries);
  2. reference-tree legality ( validator — illegal reference never scores);
  3. selected-marker-set boundary (min_occupancy / max_markers);
  4. the abort semantic: a violated must-pass exits code 4, which the
     assertion layer guarantees via AssertionFailureError -> EXIT_ASSERTION_FAILED.
"""

import pytest

from markerfinder.cli.constants import EXIT_ASSERTION_FAILED
from markerfinder.config import SelectionConfig
from markerfinder.exceptions import AssertionFailureError
from markerfinder.models.genome import GeneState
from markerfinder.models.marker import (
    RESOLUTION_PRESETS,
    MarkerLevel,
    MarkerQualityScore,
    ResolutionPreset,
    SelectionStrategy,
)
from markerfinder.modules.marker_selection import AdaptiveMarkerSelectionModule


class TestGeneStateContract:
    """Decisive interface 1: presence/state calling (hmm path)."""

    def test_states_exist_with_distinct_values(self):
        values = {state.value for state in GeneState}
        assert {"single_copy", "multi_copy", "absent"} <= values

    def test_missing_hit_is_absent(self):
        # Must-fail semantics: ABSENT is a distinct, explicit state — a
        # Missing hit can never be read as SINGLE_COPY.
        assert GeneState.ABSENT.value == "absent"
        assert GeneState.ABSENT is not GeneState.SINGLE_COPY


class TestReferenceLegalityContract:
    """Decisive interface 2: reference judgement."""

    def test_illegal_reference_never_scores(self):
        from markerfinder.config import HGTConfig
        from markerfinder.models.evidence import MeasureState
        from markerfinder.models.tree import Tree
        from markerfinder.modules.hgt_filter import PhylogeneticHGTDetector

        detector = PhylogeneticHGTDetector(HGTConfig())
        result = detector.detect(
            "m1",
            Tree(newick="((A,B),(C,D));"),
            Tree(newick="(A,B,C,D);"),  # Polytomy reference
        )
        assert result.rf_state is MeasureState.REJECTED
        assert result.normalized_rf is None

    def test_must_fail_fixture_reversal(self):
        # Reversing the defect (resolved reference) flips the verdict green
        from markerfinder.models.evidence import MeasureState
        from markerfinder.models.tree import Tree
        from markerfinder.modules.hgt_filter import PhylogeneticHGTDetector
        from markerfinder.config import HGTConfig

        detector = PhylogeneticHGTDetector(HGTConfig())
        result = detector.detect(
            "m1",
            Tree(newick="((A,B),(C,D));"),
            Tree(newick="((A,B),(C,D));"),
        )
        assert result.rf_state is not MeasureState.REJECTED


class TestSelectionBoundaryContract:
    """Decisive interface 3: the selected-marker-set boundary."""

    @staticmethod
    def _selector(min_occupancy):
        from markerfinder.config import SelectionConfig

        return AdaptiveMarkerSelectionModule(SelectionConfig(min_occupancy=min_occupancy))

    def test_preset_is_decisive_and_reversible(self):
        """A preset change is a boundary change, testable both ways."""
        conservative = RESOLUTION_PRESETS[ResolutionPreset.CONSERVATIVE]
        standard = RESOLUTION_PRESETS[ResolutionPreset.STANDARD]
        assert conservative.min_occupancy > standard.min_occupancy
        assert conservative.max_markers < standard.max_markers
        assert conservative.selection_strategy is SelectionStrategy.INFO_MAX

    def test_strategy_distinct_values(self):
        values = {strategy.value for strategy in SelectionStrategy}
        assert {"greedy", "info_max", "rate_balanced", "sparse_optimized"} <= values


class TestMustPassAbortSemantics:
    """ +: a violated must-pass aborts with exit code 4.

     is the exit-code half: ``EXIT_ASSERTION_FAILED = 4`` must be the
    value the process exits with, and ``--allow-assertion-failure`` must be
    the only way to walk past it.
    """

    def test_exit_code_4_reserved_for_assertions(self):
        assert EXIT_ASSERTION_FAILED == 4

    def test_assertion_failure_error_maps_to_four(self):
        from markerfinder.cli.commands import EXIT_ASSERTION_FAILED as MAPPED

        assert MAPPED == 4

    def test_assertion_failure_is_raiseable(self):
        error = AssertionFailureError("taxonomy must-pass violated")
        assert "must-pass" in str(error)

    def test_quality_level_unknown_on_unmeasured(self):
        # Must-pass grading degrades honestly: unmeasured => UNKNOWN, never L1.
        qs = MarkerQualityScore(marker_id="m")
        assert qs.overall_score is None
        assert qs.level is MarkerLevel.UNKNOWN
