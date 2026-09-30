""" (decision card is the judges' source) and (byte determinism).

 asks that the assertions, the threshold scan and the consistency
criterion take their data from the decision card and never re-derive it
elsewhere. The assertion state was reading ``ev.overall_risk``
directly instead — a parallel view of the same decision, which lets a check
product and the self-check drift apart silently if a back-fill touches one side.

 asks that two runs over the same input produce byte-identical decision
cards, assertion results and banding. These are the parts that can be proven
without an external tool chain; where a run product needs mafft/FastTree, the
gap is recorded rather than assumed closed.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from markerfinder.assertions import (
    AssertionReport,
    AssertionResult,
    run_assertions,
)
from markerfinder.modules.composition import write_composition_tsv
from markerfinder.modules.hgt_scan import run_threshold_scan
from markerfinder.modules.taxonomy_mustpass import (
    MustPassReport,
    MustPassViolation,
    write_mustpass_tsv,
)
from markerfinder.pipeline import MarkerFinderPipeline
from markerfinder.models.marker import MarkerLevel
from markerfinder.utils.informative_sites import apply_pis_floor

_view = MarkerFinderPipeline._assertion_view


def _ev(marker_id="M1", risk=0.1, card=None, level=MarkerLevel.LEVEL_2):
    return SimpleNamespace(
        marker_id=marker_id, overall_risk=risk, level=level,
        decision_card=card,
    )


class TestDecisionCardIsTheSource:
    def test_card_value_wins_over_the_attribute(self):
        ev = _ev(risk=0.9, card={"marker_id": "M1", "overall_risk": 0.1})
        view = _view(ev)
        assert view["overall_risk"] == 0.1
        assert view["risk_source"] == "decision_card"

    def test_card_risk_alias_is_honoured(self):
        view = _view(_ev(risk=0.9, card={"risk": 0.4}))
        assert view == {
            "marker_id": "M1", "overall_risk": 0.4,
            "risk_source": "decision_card",
        }

    def test_missing_card_falls_back_and_says_so(self):
        for card in (None, {}, {"detector": "species_tree_rf"}):
            view = _view(_ev(risk=0.25, card=card))
            assert view["overall_risk"] == 0.25
            assert view["risk_source"] == "attribute"

    def test_a_non_dict_card_cannot_crash_the_state(self):
        view = _view(_ev(risk=0.3, card="oops"))
        assert view["overall_risk"] == 0.3

    def test_the_card_governs_the_verdict_end_to_end(self):
        """MUST-FAIL CONTROL: an out-of-range risk recorded on the card must
        trip the assertion even though the attribute looks clean."""
        ev = _ev(risk=0.1, card={"marker_id": "M1", "overall_risk": -0.2})
        state = {
            "phylo_steps": [],
            "hgt_evaluations": [_view(ev)],
            "evidence_coverage": 1.0,
            "n_markers": 1,
            "discrimination_probe": None,
            "level_thresholds": {"level1_max": 0.25, "level2_max": 0.60},
            "far_active": False,
            "far_card_entries": [],
        }
        report = run_assertions(state, mode="run")
        fired = [r.assertion_id for r in report.results if r.passed is False]
        assert fired, (
            "the card said overall_risk=-0.2 and no assertion noticed — the "
            "card is not really the source of the verdict"
        )

    def test_shipped_card_actually_carries_the_risk(self):
        """Otherwise the code silently falls back to the attribute again.

        The fallback in _assertion_view is deliberately tolerant, which means a
        card that stopped carrying the key would NOT fail anything — would
        just quietly stop holding. This is the lock that prevents that.
        """
        import inspect

        from markerfinder.modules import hgt_filter

        source = inspect.getsource(hgt_filter)
        assert '"overall_risk": risk' in source, (
            "the decision card no longer carries overall_risk: the judges would "
            "read the attribute instead and the provenance rule would be unmet silently"
        )

    def test_no_component_re_reads_log_files(self):
        """'s other half: nobody may re-parse logs to recover decisions.

        Matches a string literal ending in ``.log`` specifically; a bare
        ``.log`` substring would also hit ``math.log(`` and report a false
        positive.
        """
        import inspect
        import re

        pattern = re.compile(r"[\"'][^\"']*\.log[\"']")
        # Positive control: the pattern must actually catch a violation.
        assert pattern.search('x = open("run.log")')
        assert not pattern.search("y = math.log(z)")
        for module in (
            "markerfinder.assertions", "markerfinder.modules.hgt_scan",
            "markerfinder.modules.consistency_screen",
        ):
            mod = __import__(module, fromlist=["*"])
            source = inspect.getsource(mod)
            match = pattern.search(source)
            assert match is None, (
                f"{module} names a log file ({match.group(0) if match else ''}); "
                f"decisions must come from the card, not re-read from logs"
            )


class TestByteDeterminism:
    """Same input, byte-identical products."""

    def _evaluations(self):
        return [
            _ev("M1", 0.10), _ev("M2", 0.26), _ev("M3", 0.61),
            _ev("M4", 0.05, level=MarkerLevel.UNKNOWN),
        ]

    def test_threshold_scan_is_identical_across_runs(self, tmp_path):
        left, right = tmp_path / "a", tmp_path / "b"
        for out in (left, right):
            run_threshold_scan(self._evaluations(), str(out), "mf")
        a = (left / "Phase5_reports" / "mf.threshold_scan.tsv").read_bytes()
        b = (right / "Phase5_reports" / "mf.threshold_scan.tsv").read_bytes()
        assert a == b and a

    def test_assertions_tsv_is_identical_across_runs(self, tmp_path):
        from markerfinder.modules.report_generator import write_assertions_tsv

        def build():
            return AssertionReport(
                results=[
                    AssertionResult(assertion_id="A-02", severity="fail",
                                    passed=False, detail="rf=1.7"),
                    AssertionResult(assertion_id="A-11", severity="warn",
                                    passed=True, detail="coverage 0.40"),
                ],
                mode="run",
            )

        first, second = tmp_path / "1.tsv", tmp_path / "2.tsv"
        write_assertions_tsv(build(), first)
        write_assertions_tsv(build(), second)
        assert first.read_bytes() == second.read_bytes()

    def test_composition_file_is_order_independent(self, tmp_path):
        rows = [
            {"marker_id": "M2", "rcv": 0.2, "gc_bias": None, "n_sequences": 5},
            {"marker_id": "M1", "rcv": 0.2, "gc_bias": None, "n_sequences": 5},
            {"marker_id": "M3", "rcv": 0.9, "gc_bias": None, "n_sequences": 5},
        ]
        left = tmp_path / "l"
        right = tmp_path / "r"
        write_composition_tsv(rows, str(left), "mf", warn_threshold=0.15)
        write_composition_tsv(list(reversed(rows)), str(right), "mf",
                              warn_threshold=0.15)
        body_l = (left / "Phase5_reports" / "mf.composition.tsv").read_bytes()
        body_r = (right / "Phase5_reports" / "mf.composition.tsv").read_bytes()
        assert body_l == body_r
        # Ties broken by marker id, worst outlier first.
        lines = body_l.decode("utf-8").strip().splitlines()
        assert [ln.split("\t")[0] for ln in lines[1:]] == ["M3", "M1", "M2"]

    def test_mustpass_tsv_is_identical_across_runs(self, tmp_path):
        def build():
            report = MustPassReport(requirements_file="x.yaml", markers_checked=2,
                                    groups_tested=3)
            report.violations.append(MustPassViolation(
                marker_id="M1", relation="monophyletic", level="genus",
                taxon="Alpha", members=("A1", "A2"),
            ))
            report.not_checked.append("M2: unsupported relation 'sister_of'")
            return report

        left, right = tmp_path / "l", tmp_path / "r"
        a = write_mustpass_tsv(build(), str(left), "mf")
        b = write_mustpass_tsv(build(), str(right), "mf")
        assert a.read_bytes() == b.read_bytes()

    def test_pis_floor_result_is_repeatable(self):
        def once():
            ev = SimpleNamespace(marker_id="M1", level=MarkerLevel.LEVEL_1,
                                 notes="", decision_card={})
            demoted = apply_pis_floor([ev], {"M1": 2}, 5)
            return demoted, ev.notes, json.dumps(ev.decision_card, sort_keys=True)

        assert once() == once()

    def test_run_mode_assertions_repeat_exactly(self):
        state = {
            "phylo_steps": [],
            "hgt_evaluations": [{"marker_id": "M1", "overall_risk": 0.1,
                                 "risk_source": "attribute"}],
            "evidence_coverage": 1.0,
            "n_markers": 1,
            "discrimination_probe": None,
            "level_thresholds": {"level1_max": 0.25, "level2_max": 0.60},
            "far_active": False,
            "far_card_entries": [],
        }
        first = run_assertions(state, mode="run")
        second = run_assertions(state, mode="run")
        assert [(r.assertion_id, r.passed, r.detail) for r in first.results] == [
            (r.assertion_id, r.passed, r.detail) for r in second.results
        ]
        # Sanity that the repetition proved something: at least one assertion
        # Actually reported, otherwise two identical empty lists also "pass".
        assert first.results
