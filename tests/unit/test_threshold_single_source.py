"""The conflict-band thresholds must be one configurable source, not scenery.

's wording is "``ConflictDetector.__init__`` 的 ``rf_threshold`` /
``quartet_threshold`` 要么接入决策，要么删除". Half of it was scenery:

* the module always constructed ``ConflictDetector`` with defaults, so no run
  could reach either knob;
* ``quartet_threshold`` had **zero readers** in the whole package — a documented
  parameter that turns nothing;
* the note's lower band was recomputed as ``rf_threshold / 3.0``, i.e. a second
  source of truth for the same 0.1 that ``recommend_tree`` reads from
  ``config.recommend_rf_low``. Two sources drift, and the report prose would have
  quietly disagreed with the recommendation.

Both bounds now come from config, and ``quartet_threshold`` is deleted. These
tests are boundary tests, because a threshold that is only tested at its centre
passes even when the comparison is off by one band.
"""

from __future__ import annotations

import pytest

from markerfinder.config import PhylogeneticConfig
from markerfinder.models.pipeline_types import (
    CoalescentResult,
    ConflictReport,
    PhylogeneticResult,
    SupermatrixResult,
)
from markerfinder.modules.phylogenetic_inference import (
    ConflictDetector,
    PhylogeneticInferenceModule,
)
from markerfinder.models.tree import Tree


def _note(detector: ConflictDetector, norm_rf: float) -> str:
    return detector._topology_note("((A,B),(C,D));", "((A,B),(C,D));", norm_rf, None)


def test_control_default_bands_classify_low_moderate_strong():
    detector = ConflictDetector(config=PhylogeneticConfig())
    assert _note(detector, 0.05).startswith("Low topological conflict")
    assert _note(detector, 0.20).startswith("Moderate topological conflict")
    assert _note(detector, 0.50).startswith("Strong topological conflict")
    # Exact boundaries: the band starts where the config says it starts.
    assert _note(detector, 0.099).startswith("Low")
    assert _note(detector, 0.1).startswith("Moderate")
    assert _note(detector, 0.299).startswith("Moderate")
    assert _note(detector, 0.3).startswith("Strong")


def test_control_config_actually_moves_the_note():
    """Same measurement, different config -> different sentence.

    Without this, restoring a hardcoded 0.1/0.3 would keep every other test green.
    """
    config = PhylogeneticConfig(recommend_rf_low=0.3, recommend_rf_moderate=0.5)
    detector = ConflictDetector(config=config)
    assert _note(detector, 0.35).startswith("Moderate"), (
        "0.35 is 'Strong' under the defaults; the config must be what decides"
    )
    assert _note(detector, 0.25).startswith("Low")


def test_control_explicit_override_still_wins():
    detector = ConflictDetector(config=PhylogeneticConfig(), rf_threshold=0.05)
    assert detector.rf_threshold == 0.05
    assert _note(detector, 0.2).startswith("Strong")


def test_low_bound_is_the_same_field_the_recommendation_reads():
    """One source of truth: the note and recommend_tree cannot disagree."""
    config = PhylogeneticConfig(recommend_rf_low=0.22)
    detector = ConflictDetector(config=config)
    assert detector.low_bound == pytest.approx(0.22)
    assert _note(detector, 0.2).startswith("Low")
    assert _note(detector, 0.25).startswith("Moderate")


def test_quartet_threshold_knob_is_gone():
    """'s 'either wire it or delete it': it was deleted, and stays deleted."""
    detector = ConflictDetector(config=PhylogeneticConfig())
    assert not hasattr(detector, "quartet_threshold")
    with pytest.raises(TypeError):
        ConflictDetector(quartet_threshold=0.7)


def test_module_hands_its_own_config_to_the_detector():
    config = PhylogeneticConfig(recommend_rf_moderate=0.45)
    module = PhylogeneticInferenceModule(config)
    assert module.conflict_detector is not None
    assert module.conflict_detector.rf_threshold == pytest.approx(0.45)
    assert module.conflict_detector.low_bound == pytest.approx(config.recommend_rf_low)


def _recommend(low: float, moderate: float, rf: float):
    config = PhylogeneticConfig(recommend_rf_low=low, recommend_rf_moderate=moderate)
    module = PhylogeneticInferenceModule(config)
    result = PhylogeneticResult(
        supermatrix=SupermatrixResult(tree=Tree(newick="((A,B),(C,D));")),
        coalescent=CoalescentResult(
            species_tree=Tree(newick="((A,B),(C,D));"), species_tree_source="astral",
        ),
        conflict_report=ConflictReport(
            rf_distance=4, normalized_rf=rf, quartet_agreement=0.9
        ),
    )
    return module.recommend_tree(result), config


def test_recommendation_bands_follow_the_same_config():
    """Two-sided, and measured rather than guessed.

    normalized RF 0.2 is "moderate -> medium" under the shipped bands and
    "low -> high" once ``recommend_rf_low`` moves above it: one number, two
    configs, two verdicts. (The direction here was measured first; asserting a
    remembered direction is how a test ends up locking in the opposite of the
    implementation.)
    """
    moderate_verdict, _ = _recommend(0.1, 0.3, 0.2)
    low_verdict, _ = _recommend(0.3, 0.5, 0.2)
    assert moderate_verdict.confidence == "medium"
    assert "Moderate conflict" in moderate_verdict.reason
    assert low_verdict.confidence == "high"
    assert "Low conflict" in low_verdict.reason

    strong, _ = _recommend(0.1, 0.3, 0.6)
    assert strong.confidence == "low"
    assert "Strong conflict" in strong.reason


def test_recommendation_and_topology_note_never_disagree_about_the_band():
    """The single-source property exists to guarantee.

    The recommendation and the human-readable topology note are produced by
    different classes; before this change the note recomputed its lower bound as
    ``rf_threshold / 3``, so the two could name different bands for one number.
    """
    for low, moderate in ((0.1, 0.3), (0.25, 0.6), (0.4, 0.45)):
        config = PhylogeneticConfig(
            recommend_rf_low=low, recommend_rf_moderate=moderate
        )
        detector = ConflictDetector(config=config)
        module = PhylogeneticInferenceModule(config)
        for rf in (0.05, low, 0.2, moderate, 0.8):
            note = detector._topology_note("((A,B),(C,D));", "((A,B),(C,D));", rf, None)
            recommendation, _ = _recommend(low, moderate, rf)
            band = note.split(" ")[0].rstrip(".").capitalize()  # Low/Moderate/Strong
            assert band in recommendation.reason, (
                f"low={low} moderate={moderate} rf={rf}: note says {note!r} but the "
                f"recommendation says {recommendation.reason!r}"
            )
