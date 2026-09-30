"""Reference-tree legality tests.

The Error-A fixture (``collapsed_sponge_other.nwk``) is the retraction's
merged-reference defect: the two focal groups the competing hypotheses
require to be independent are already one clade in the reference.
"""

from pathlib import Path

import pytest

from markerfinder.config import HGTConfig
from markerfinder.models.evidence import MeasureState
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.pipeline_types import PhyloStepResult
from markerfinder.models.tree import Tree
from markerfinder.modules.hgt_filter import HGTDecisionEngine, PhylogeneticHGTDetector
from markerfinder.utils.reference import (
    CLADE_COLLAPSED,
    MALFORMED,
    POLYTOMY,
    TIPSET_MISMATCH,
    TOO_FEW_TIPS,
    validate_reference_tree,
)

FIXTURES = Path(__file__).parent.parent / "fixtures" / "reference"


def _fixture(name: str) -> str:
    return (FIXTURES / name).read_text().strip()


class TestReferenceLegality:
    def test_legal_reference_passes(self):
        verdict = validate_reference_tree(
            _fixture("resolved_ref_4tips.nwk"), "((A,B),(C,D));"
        )
        assert verdict.ok
        assert verdict.code is None

    def test_polytomy_reference_rejected(self):
        verdict = validate_reference_tree(
            _fixture("polytomy_root.nwk"), "(A,B);"
        )
        assert not verdict.ok
        assert verdict.code == POLYTOMY

    def test_unfolding_polytomy_makes_it_legal(self):
        # Must-fail reversal: resolve the star → the reference turns legal.
        verdict = validate_reference_tree(
            _fixture("resolved_ref_4tips.nwk"), "(A,B);"
        )
        assert verdict.ok

    def test_tipset_mismatch_rejected(self):
        verdict = validate_reference_tree(
            _fixture("resolved_ref_4tips.nwk"), "((A,B),(X,Y));"
        )
        assert not verdict.ok
        assert verdict.code == TIPSET_MISMATCH

    def test_three_tip_reference_rejected(self):
        verdict = validate_reference_tree(
            _fixture("three_tip_ref.nwk"), "(A,B);"
        )
        assert not verdict.ok
        assert verdict.code == TOO_FEW_TIPS

    def test_malformed_reference_rejected(self):
        verdict = validate_reference_tree("((A,B;", "(A,B);")
        assert not verdict.ok
        assert verdict.code == MALFORMED

    # Error A's mechanism -- focal groups already collapsed in the
    # Reference -- is refused outright. cover the sibling cases
    # (coverage and the per-rule reason) in this class.
    def test_collapsed_target_clade_rejected(self):
        # Error A reproduction: focal groups {S1,S2} and {O1} are merged on
        # A single reference edge, so the reference cannot arbitrate.
        verdict = validate_reference_tree(
            _fixture("collapsed_sponge_other.nwk"),
            _fixture("collapsed_sponge_other.nwk"),
            target_clades=[{"S1", "S2"}, {"O1"}],
        )
        assert not verdict.ok
        assert verdict.code == CLADE_COLLAPSED

    def test_uncollapsed_reference_with_same_clades_passes(self):
        verdict = validate_reference_tree(
            _fixture("uncollapsed_sponge_other.nwk"),
            _fixture("uncollapsed_sponge_other.nwk"),
            target_clades=[{"S1", "S2"}, {"O1", "O2"}],
        )
        assert verdict.ok

    def test_superset_reference_is_accepted(self):
        # Requires tipset(ref) ⊇ tipset(gene) — a superset is legal.
        verdict = validate_reference_tree(
            _fixture("superset_ref.nwk"), _fixture("resolved_ref_4tips.nwk")
        )
        assert verdict.ok


class TestDetectorIntegration:
    def test_illegal_reference_gives_unknown_not_risk(self):
        det = PhylogeneticHGTDetector(HGTConfig())
        gene = Tree(newick="((A,B),(C,D));")
        polytomy_ref = Tree(newick=_fixture("polytomy_root.nwk"))
        res = det.detect("m1", gene, polytomy_ref)
        assert res.rf_state is MeasureState.REJECTED
        assert res.quartet_state is MeasureState.REJECTED
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m1", res)
        assert ev.level == MarkerLevel.UNKNOWN
        assert "POLYTOMY" in ev.notes or "reference" in ev.notes

    def test_legal_reference_still_scores(self):
        det = PhylogeneticHGTDetector(HGTConfig())
        gene = Tree(newick="((A,B),(C,D));")
        ref = Tree(newick="((A,B),(C,D));")
        res = det.detect("m1", gene, ref)
        assert res.risk_basis == "two_signal_weighted"
        assert res.normalized_rf is not None


class TestMonophylyProvenance:
    def test_rank_fallback_recorded(self):
        """Test_level_used/n_total/n_mono land on the result."""
        det = PhylogeneticHGTDetector(HGTConfig())
        taxonomy_map = {
            "g1": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G1"},
            "g2": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G1"},
            "g3": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G2"},
            "g4": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G2"},
        }
        nwk = "((g1,g2),(g3,g4));"
        res = det.detect_monophyly("m1", nwk, taxonomy_map, level="genus", threshold=0.5)
        if res is None:
            pytest.skip(
                "ete3 unavailable: MAD rooting not measurable — "
                "ASSERTION NOT EXECUTED in this interpreter"
            )
        assert res.test_level_used is not None
        assert res.n_total >= 1
        assert res.n_mono is not None

    def test_cross_rank_not_compared(self):
        """A proportion measured at a fallback rank must not be graded
        against the configured threshold — it degrades to UNKNOWN."""
        det = PhylogeneticHGTDetector(HGTConfig())
        taxonomy_map = {
            "g1": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G1"},
            "g2": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G1"},
            "g3": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G2"},
            "g4": {"domain": "D", "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "G2"},
        }
        nwk = "((g1,g2),(g3,g4));"
        res = det.detect_monophyly("m1", nwk, taxonomy_map, level="species", threshold=0.5)
        if res is None:
            pytest.skip(
                "ete3 unavailable: MAD rooting not measurable — "
                "ASSERTION NOT EXECUTED in this interpreter"
            )
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker("m1", res)
        if res.cross_rank_comparison:
            assert ev.level == MarkerLevel.UNKNOWN
            assert "cross_rank_incomparable" in ev.notes
        else:
            assert ev.level == MarkerLevel.LEVEL_1

    def test_rejected_reference_result_carries_reason(self):
        res = PhyloStepResult(gene_id="m")
        res.rf_state = MeasureState.REJECTED
        res.rf_reason = "reference illegal [POLYTOMY]"
        unmet = res.unmeasured_components()
        assert unmet and unmet[0].value == "reference_illegal"
