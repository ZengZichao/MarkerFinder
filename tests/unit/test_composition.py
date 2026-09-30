"""Composition/GC diagnostics tests.

The prohibition is locked here: composition evidence is a parallel
column and must never merge into the HGT risk or the quality score.
"""

import pytest

from markerfinder.modules.composition import (
    aa_composition,
    chi2_vs_group_mean,
    gc_content,
    gc_or_codon_bias,
    rcv,
)


class TestRcv:
    def test_rcv_zero_for_homogeneous(self):
        seqs = ["ACDEFGHIKLMNPQRSTVWY"] * 3
        assert rcv(seqs) == pytest.approx(0.0, abs=1e-9)

    def test_rcv_flags_planted_outlier(self):
        base = "ACDEFGHIKLMNPQRSTVWY"
        acidic = "DEDEDEDEDEDEDEDEDEDE"  # Composition-skewed outlier
        seqs = [base, base, base, acidic]
        with_outlier = rcv(seqs)
        without = rcv([base, base, base, base])
        assert with_outlier > without
        assert with_outlier > 0.2

    def test_group_labels_worst_deviation(self):
        base = "ACDEFGHIKLMNPQRSTVWY"
        skewed = "DEDEDEDEDEDEDEDEDEDE"
        seqs = [base, base, skewed, skewed]
        labels = ["g1", "g1", "g2", "g2"]
        assert rcv(seqs, labels) > 0.0


class TestGc:
    def test_protein_input_gc_not_applicable(self):
        #.faa input has no GC: explicit None (NOT_APPLICABLE), never 0.
        assert gc_content("ACDEFGHIKLMNPQRSTVWY", protein=True) is None
        assert gc_or_codon_bias(["ACDEF", "GHIKL"], protein=True) is None

    def test_nucleotide_gc_spread(self):
        at_rich = ["AAAA", "AAAA", "AAAA"]
        mixed = ["GGGG", "AAAA", "GGGG", "AAAA"]
        assert gc_or_codon_bias(at_rich, protein=False) == pytest.approx(0.0)
        assert gc_or_codon_bias(mixed, protein=False) == pytest.approx(0.5)


class TestChi2:
    def test_outlier_scores_higher(self):
        base = "ACDEFGHIKLMNPQRSTVWY"
        group = [base, base, base]
        assert chi2_vs_group_mean("DEDEDEDEDEDE", group) > chi2_vs_group_mean(base, group)

    def test_empty_inputs_return_none(self):
        assert chi2_vs_group_mean("", ["ACDEF"]) is None
        assert chi2_vs_group_mean("ACDEF", []) is None


class TestNfr05NoMerging:
    def test_composition_not_merged_into_overall_score(self):
        """Composition metrics live outside the quality/risk pipeline."""
        from dataclasses import fields
        from markerfinder.models.marker import MarkerQualityScore
        from markerfinder.models.pipeline_types import HGTEvaluation

        score_fields = {f.name for f in fields(MarkerQualityScore)}
        eval_fields = {f.name for f in fields(HGTEvaluation)}
        assert "composition_rcv" not in score_fields
        assert "composition_rcv" not in eval_fields
        # Overall_score weights (0.25/0.25/0.10/0.20/0.20) sum to 1.0 without
        # Any composition component:
        qs = MarkerQualityScore(
            hmm_score=100.0, occupancy=1.0, length_ratio=1.0,
            phylogenetic_informativeness=1.0, hgt_risk_score=0.0,
        )
        # Weights sum: 0.25+0.25+0.10+0.20+0.20 = 1.00
        assert qs.overall_score == pytest.approx(1.0)

    def test_aa_composition_sums_to_one(self):
        comp = aa_composition("ACDEFGHIKLMNPQRSTVWY")
        assert sum(comp.values()) == pytest.approx(1.0)
        assert aa_composition("---XXX") == {}
