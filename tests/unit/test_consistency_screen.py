"""Marker-level cross-consistency screen tests."""

import pytest

from markerfinder.exceptions import UnsupportedCriterion
from markerfinder.models.evidence import Measurement
from markerfinder.modules.consistency_screen import (
    BENCHMARK_EXPECTED_TRUTH,
    BOUNDS,
    ConsistencyGrade,
    MarkerConsistency,
    guard,
    screen,
    write_excluded_profile,
)

# Same 4-taxon set, two reference topologies.
GENE_AB = "((A,B),(C,D));"          # Gene tree "supports" the AB|CD pairing
CONCAT_AB = "((A,B),(C,D));"        # Concat agrees with the gene tree
ASTRAL_CD = "((A,C),(B,D));"        # Astral opposes it


class TestScreen:
    def test_both_sides_same_orientation_above_threshold_is_consistent(self):
        results = screen(
            {"M1": GENE_AB}, CONCAT_AB, CONCAT_AB, stringency=5,
        )
        assert results[0].grade is ConsistencyGrade.CONSISTENT

    def test_opposite_sides_is_inconsistent(self):
        results = screen(
            {"M1": GENE_AB}, CONCAT_AB, ASTRAL_CD, stringency=5,
        )
        # Gene tree agrees with concat, disagrees with astral at max stringency
        assert results[0].grade is ConsistencyGrade.INCONSISTENT

    def test_one_side_weak_is_inconclusive(self):
        # A star gene tree carries no orientation: both sides unresolved.
        results = screen(
            {"M1": "(A,B,C,D);"}, CONCAT_AB, CONCAT_AB, stringency=5,
        )
        assert results[0].grade is ConsistencyGrade.INCONCLUSIVE

    def test_illegal_reference_never_scores(self):
        results = screen(
            {"M1": GENE_AB}, "(A,B);", CONCAT_AB, stringency=5,
        )
        assert results[0].grade is ConsistencyGrade.INCONCLUSIVE
        assert any("reference_illegal" in r for r in results[0].reasons)

    def test_insufficient_sites_is_inconclusive(self):
        results = screen(
            {"M1": GENE_AB}, CONCAT_AB, CONCAT_AB,
            stringency=5, min_sites=100, pis_map={"M1": 3},
        )
        assert results[0].grade is ConsistencyGrade.INCONCLUSIVE
        assert "insufficient_sites" in results[0].reasons

    def test_stringency_ladder_monotone(self):
        """ — stricter bounds never turn inconsistent into
        consistent when a weaker band already refused it."""
        gene = "((A,B),(C,D));"   # Agrees but only moderately strong
        concat = "(((A,B),C),D);" # Moderately different
        astral = "(((A,B),D),C);"
        grades = []
        for s in sorted(BOUNDS):
            results = screen({"M1": gene}, concat, astral, stringency=s)
            grades.append(results[0].grade)
        inconsistent_index = [
            i for i, g in enumerate(grades) if g is ConsistencyGrade.INCONSISTENT
        ]
        consistent_index = [
            i for i, g in enumerate(grades) if g is ConsistencyGrade.CONSISTENT
        ]
        if consistent_index:
            assert min(consistent_index) <= min(inconsistent_index or [len(grades)])


class TestGate:
    def test_gate_blocks_when_prereq_missing(self, monkeypatch):
        """With truth data absent, consistency/hybrid refuse to run.
        Simulated by stubbing out one probe: the benchmark truth YAML now exists
        on disk, so the gate unlocked on its own, which is the documented
        behaviour."""
        import markerfinder.modules.consistency_screen as cs

        original_probe = cs._probe
        monkeypatch.setattr(
            cs, "_probe",
            lambda key: False if key == BENCHMARK_EXPECTED_TRUTH else original_probe(key),
        )
        with pytest.raises(UnsupportedCriterion) as excinfo:
            guard("consistency")
        assert BENCHMARK_EXPECTED_TRUTH in str(excinfo.value)

    def test_gate_unlocks_with_expected_yaml_present(self):
        """With the truth YAML in place the gate passes —
        automatic unlock, no human switch."""
        guard("consistency")  # Must NOT raise now that expected/*.yaml exists

    def test_gate_risk_always_passes(self):
        guard("risk")  # Must never raise

    def test_gate_cannot_be_faked_by_flag(self, monkeypatch):
        """AD: a completion flag cannot unlock the gate when the
        artifact probe fails."""
        import markerfinder.modules.consistency_screen as cs

        class FakeCfg:
            wp11_done = True  # The flag an impostor would add

        monkeypatch.setattr(cs, "_probe", lambda key: False)
        with pytest.raises(UnsupportedCriterion):
            guard("consistency", FakeCfg())


class TestStatisticsDiscipline:
    def test_no_chisquare_on_markers(self):
        """No marker-as-independent-unit chi-square exists in the
        consistency screen — the module exposes bootstrap/jackknife-free
        descriptive grading only."""
        import inspect
        from markerfinder.modules import consistency_screen

        source = inspect.getsource(consistency_screen)
        assert "chi2" not in source.lower().replace("chi2_vs_group_mean", "")
        assert "chisquare" not in source.lower()

    def test_excluded_profile_written(self, tmp_path):
        results = screen(
            {"M1": GENE_AB, "M2": "(A,B,C,D);"}, CONCAT_AB, ASTRAL_CD, stringency=5,
        )
        path = write_excluded_profile(
            results, str(tmp_path), "mf", pis_map={"M1": 10, "M2": 2},
        )
        content = path.read_text(encoding="utf-8")
        lines = content.strip().splitlines()
        assert lines[0].startswith("marker_id\tgrade\tstringency\treasons\tpis")
        rows = {line.split("\t")[0] for line in lines[1:]}
        assert "M1" in rows  # Inconsistent => profiled
        assert rows  # Excluded markers are never silently dropped

    def test_risk_mode_default_output_unchanged(self):
        """Default mode 'risk' means the screen never runs — guard('risk')
        passes trivially and the grading path is the two-signal engine."""
        from markerfinder.config import HGTConfig
        from markerfinder.modules.hgt_filter import HGTDecisionEngine
        from markerfinder.models.pipeline_types import PhyloStepResult

        guard("risk")
        ev = HGTDecisionEngine(HGTConfig()).evaluate_marker(
            "m", PhyloStepResult(overall_risk=0.1)
        )
        assert ev.level.value == "level_1"
