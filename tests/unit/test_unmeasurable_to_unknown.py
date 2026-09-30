"""Unmeasurable → UNKNOWN three-state tests.

The zero-discrimination experiment is locked here as a must-fail control: if anyone
reintroduces a numeric placeholder for an unmeasured signal, these tests
turn red.
"""

import pytest

from markerfinder.config import HGTConfig
from markerfinder.models.evidence import MeasureState, UnknownReason
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.pipeline_types import PhyloStepResult
from markerfinder.models.tree import Tree
from markerfinder.modules.hgt_filter import (
    HGTDecisionEngine,
    PhylogeneticHGTDetector,
)
from markerfinder.utils.hgt_utils import combine_hgt_scores


def _ete3_available() -> bool:
    try:
        import ete3  # Noqa: F401
        return True
    except Exception:
        return False


class TestQuartetUnmeasurable:
    def test_tips_lt4_yields_unknown(self):
        """One unmeasured critical component blocks the whole grading.

        The marker goes through the existing ``MarkerLevel.UNKNOWN`` channel with
        ``risk_basis == 'unscreened'`` instead of being scored from whichever half
        happened to measure — the rule that keeps a partially-measured marker out
        of Level 1/2/3 ( extended to partial measurability).
        """
        det = PhylogeneticHGTDetector(HGTConfig())
        gene = Tree(newick="(A,B);")
        ref = Tree(newick="((A,B),(C,D));")
        res = det.detect("m1", gene, ref)
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m1", res)
        assert ev.level == MarkerLevel.UNKNOWN
        assert ev.risk_basis == "unscreened"
        assert "quartet" in ev.notes or "rf" in ev.notes

    def test_zero_measurable_quartets_yields_none_with_reason(self):
        # A gene tree with 2 tips and a reference with disjoint tips: RF
        # Still computes (0 shared), quartet has nothing to evaluate.
        det = PhylogeneticHGTDetector(HGTConfig())
        gene = Tree(newick="(X,Y);")
        ref = Tree(newick="((A,B),(C,D));")
        q, reason = det._calculate_quartet_consistency(gene, ref)
        assert q is None
        assert reason

    def test_identical_trees_still_measure(self):
        det = PhylogeneticHGTDetector(HGTConfig())
        t = Tree(newick="((A,B),(C,D));")
        res = det.detect("m1", t, t)
        assert res.normalized_rf == 0.0
        assert res.quartet_consistency == 1.0
        assert res.risk_basis == "two_signal_weighted"


class TestCombineSignals:
    def test_empty_scores_returns_none_not_zero(self):
        risk, n_used = combine_hgt_scores({}, {"rf": 0.5, "quartet": 0.5})
        assert risk is None
        assert n_used == 0

    def test_single_signal_reports_n_used(self):
        risk, n_used = combine_hgt_scores({"rf": 0.5}, {"rf": 0.5, "quartet": 0.5})
        assert risk == 0.5
        assert n_used == 1

    def test_single_signal_not_graded(self):
        # Synthesis over one signal must not auto-grade. detect
        # With an unmeasurable quartet therefore lands in UNKNOWN, with the
        # Distinguishable SINGLE_SIGNAL_ONLY-flavoured reason attached.
        det = PhylogeneticHGTDetector(HGTConfig())
        gene = Tree(newick="(A,B);")
        ref = Tree(newick="((A,B),(C,D));")
        res = det.detect("m1", gene, ref)
        assert res.unmeasured_components()  # Quartet unmeasurable here
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m1", res)
        assert ev.level == MarkerLevel.UNKNOWN


class TestEvaluateMarkerUnknown:
    def test_rf_none_yields_unknown_not_level2(self):
        # Zero-discrimination must-fail control: an unmeasured RF must NOT produce the old
        # Placeholder flow (0.5 → level_2 with high confidence).
        res = PhyloStepResult(
            gene_id="m",
            rf_distance=None,
            normalized_rf=None,
            quartet_consistency=0.8,
            rf_state=MeasureState.NOT_MEASURABLE,
            rf_reason="synthetic: rf unmeasurable",
        )
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m", res)
        assert ev.level == MarkerLevel.UNKNOWN
        assert ev.confidence == "unknown"
        assert "rf" in ev.notes

    def test_quartet_none_yields_unknown(self):
        res = PhyloStepResult(
            gene_id="m",
            rf_distance=2,
            normalized_rf=1.0,
            quartet_consistency=None,
            quartet_state=MeasureState.NOT_MEASURABLE,
            quartet_reason="synthetic: no quartet evaluated",
        )
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m", res)
        assert ev.level == MarkerLevel.UNKNOWN
        assert "quartet" in ev.notes

    def test_unknown_still_not_excluded(self):
        # Pipeline exclusion rule is {level == LEVEL_3} — UNKNOWN is kept.
        res = PhyloStepResult(gene_id="m")
        res.rf_state = MeasureState.NOT_MEASURABLE
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m", res)
        excluded = {ev.level} == {MarkerLevel.LEVEL_3}
        assert not excluded
        assert ev.level == MarkerLevel.UNKNOWN

    def test_unknown_sources_distinguishable(self):
        # UNKNOWN reasons must be machine-distinguishable.
        engine = HGTDecisionEngine(HGTConfig())
        ev_none = engine.evaluate_marker("m1")
        assert "no_reference" in ev_none.notes
        res = PhyloStepResult(gene_id="m2")
        res.rf_state = MeasureState.NOT_MEASURABLE
        res.rf_reason = "tip cap exceeded"
        ev_rf = engine.evaluate_marker("m2", res)
        assert "rf_unmeasurable" in ev_rf.notes
        assert ev_none.notes != ev_rf.notes

    def test_confidence_not_synonym_of_level(self):
        # Confidence tracks distance to the grading boundary.
        engine = HGTDecisionEngine(HGTConfig())
        # Risk=0.26 sits just above the L1/L2 boundary (0.25) → low.
        near = PhyloStepResult(overall_risk=0.26)
        ev_near = engine.evaluate_marker("m", near)
        assert ev_near.level == MarkerLevel.LEVEL_2
        assert ev_near.confidence == "low"
        # Risk=0.0 is far from both boundaries → high.
        far = PhyloStepResult(overall_risk=0.0)
        ev_far = engine.evaluate_marker("m", far)
        assert ev_far.level == MarkerLevel.LEVEL_1
        assert ev_far.confidence == "high"

    def test_evidence_coverage_reported(self):
        engine = HGTDecisionEngine(HGTConfig())
        good = PhyloStepResult(overall_risk=0.1)
        bad = PhyloStepResult(gene_id="bad")
        bad.rf_state = MeasureState.NOT_MEASURABLE
        report = engine.generate_hgt_report(
            [
                engine.evaluate_marker("a", good),
                engine.evaluate_marker("b", good),
                engine.evaluate_marker("c", bad),
            ]
        )
        assert report.n_actually_measured == 2
        assert report.evidence_coverage == pytest.approx(2 / 3)
        assert report.unknown_reason_counts.get("rf_unmeasurable") == 1

    def test_quality_score_none_when_hgt_unmeasured(self):
        # UNKNOWN markers must not contribute a 0.0 "clean" HGT risk
        # To the quality score — pipeline back-fills None for UNKNOWN.
        from markerfinder.models.marker import MarkerQualityScore

        qs = MarkerQualityScore(marker_id="m", phylogenetic_informativeness=0.8)
        qs.hgt_risk_score = None  # What pipeline.py now does for UNKNOWN
        assert qs.overall_score is None  # Composite honestly unavailable


@pytest.mark.skipif(
    not _ete3_available(),
    reason="the full-path zero-discrimination check needs live tree comparison; the pure-Python "
    "path is covered by test_identical_trees_still_measure — NOT EXECUTED "
    "when ete3 is unavailable",
)
class TestZeroDiscriminationLivePath:
    def test_e0_discrimination_with_ete3(self):
        cfg = HGTConfig()
        det = PhylogeneticHGTDetector(cfg)
        eng = HGTDecisionEngine(cfg)
        ref = Tree(newick="((A,B),(C,D));")
        clash = Tree(newick="((A,C),(B,D));")
        a = eng.evaluate_marker("clash", det.detect("clash", clash, ref))
        b = eng.evaluate_marker("same", det.detect("same", ref, ref))
        assert a.overall_risk != b.overall_risk
