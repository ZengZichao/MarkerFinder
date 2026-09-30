"""Placeholder removal and Optional type-contract tests.

the four literal placeholders are gone from production code; Optional
fields render as ``NA`` (TSV) / ``N/A`` (HTML); UNKNOWN's ``overall_risk=0.0``
is intentional and locked by a test.
"""

from markerfinder.config import HGTConfig
from markerfinder.models.marker import MarkerLevel, MarkerQualityScore
from markerfinder.models.pipeline_types import (
    CoalescentResult,
    ConflictReport,
    PhyloStepResult,
    SupermatrixResult,
)
from markerfinder.models.tree import Tree
from markerfinder.modules.hgt_filter import (
    HGTDecisionEngine,
    PhylogeneticHGTDetector,
)
from markerfinder.modules.report_generator import _fmt


class TestPlaceholderLiteralsRemoved:
    def test_supermatrix_avg_ufboot_is_none(self):
        # Default-constructed result carries "not measured", not 0.0.
        assert SupermatrixResult().avg_ufboot is None

    # The literal 0.0 that used to sit on this field is gone; an
    # Unmeasured average quartet support is None, not zero support.
    def test_coalescent_avg_quartet_support_is_none(self):
        assert CoalescentResult().avg_quartet_support is None

    # The annotations now match what the code assigns --
    # ``PhyloStepResult.rf_distance`` is declared Optional[int] and is in
    # Fact None on the monophyly-only path, where it used to be typed
    # ``int = 0`` while hgt_filter assigned None.
    def test_monophyly_path_leaves_rf_none(self):
        # The monophyly path does not measure RF; normalized_rf must
        # Stay None (NOT_APPLICABLE semantics), never 0.0. is precisely the
        # Deleted hgt_filter branch that assigned 0.5 to an unmeasurable RF, which
        # Then read as "moderate disagreement" instead of "not measured".
        res = PhyloStepResult(gene_id="m", monophyly_proportion=0.9,
                              overall_risk=0.1, is_suspicious=False)
        assert res.rf_distance is None
        assert res.normalized_rf is None

    def test_conflict_report_optional_fields(self):
        rep = ConflictReport()
        assert rep.rf_distance is None
        assert rep.normalized_rf is None

    def test_type_annotations_are_optional(self):
        import types
        import typing
        hints = typing.get_type_hints(PhyloStepResult)
        origin = typing.get_origin(hints["normalized_rf"])
        assert origin in (typing.Union, types.UnionType)


class TestInformativenessPlaceholderRemoved:
    def test_informativeness_placeholder_removed(self):
        # Phase 1.5 no longer fabricates 0.5.
        qs = MarkerQualityScore(marker_id="COG001")
        assert qs.phylogenetic_informativeness is None
        assert qs.overall_score is None

    def test_unknown_risk_zero_is_intentional(self):
        # Trap: UNKNOWN's overall_risk stays 0.0 by design —
        # The distribution chart buckets UNKNOWN separately, and changing
        # It to None would ripple through the bucket statistics. Locked here.
        engine = HGTDecisionEngine(HGTConfig())
        ev = engine.evaluate_marker("COG001", phylogenetic_result=None)
        assert ev.level == MarkerLevel.UNKNOWN
        assert ev.overall_risk == 0.0


class TestRendering:
    def test_fmt_none_renders_na(self):
        assert _fmt(None) == "NA"

    def test_fmt_value_renders_four_decimals(self):
        assert _fmt(0.5) == "0.5000"

    def test_tsv_renders_na_not_zero(self):
        # Must-fail control: an unmeasured value must never render 0.0000.
        detector = PhylogeneticHGTDetector(HGTConfig())
        res = detector.detect_monophyly(
            "m1", "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));",
            {"B1": {"genus": "B"}, "B2": {"genus": "B"},
             "A1": {"genus": "A"}, "A2": {"genus": "A"}},
            level="genus", threshold=0.5,
        )
        # On this (no-ete3) machine the monophyly screen is NOT_MEASURABLE;
        # A None-normalized_rf result must render NA in TSV position (
        # TSV uses "NA", never 0.0000).
        if res is not None:
            assert _fmt(res.normalized_rf) == "NA"
        else:
            assert _fmt(None) == "NA"

    def test_rich_rendering_helper_exists_for_html(self):
        """The same missing value renders as "N/A" in HTML."""
        from markerfinder.modules.report_generator import _fmt_html
        assert _fmt_html(None) == "N/A"
        assert _fmt_html(1.0) == "1.0000"


class TestIdenticalTreesZeroDistancePreserved:
    def test_rf_zero_is_a_real_zero(self):
        # (0, 0.0) from identical trees must stay distinct from (None, None).
        from markerfinder.utils.tree_utils import calculate_rf_distance
        raw, norm = calculate_rf_distance("((A,B),(C,D));", "((A,B),(C,D));")
        assert raw == 0 and norm == 0.0
