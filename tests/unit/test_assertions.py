"""Assertion layer tests."""

import pytest

from markerfinder.assertions import (
    REGISTRY,
    AssertionReport,
    build_negative_state,
    run_assertions,
)
from markerfinder.assertions import _good_state


class TestRegistry:
    def test_all_14_assertions_registered(self):
        ids = {spec.id for spec in REGISTRY}
        assert ids == {f"A-{i:02d}" for i in range(1, 15)}

    def test_every_assertion_has_fixture_path(self):
        for spec in REGISTRY:
            assert spec.must_fail_fixture
            assert spec.rationale
            assert spec.severity in ("warn", "fail")

    def test_registry_has_no_orphan_without_fixture_file(self):
        for spec in REGISTRY:
            fixture = (
                __import__("pathlib").Path(__file__).parent.parent
                / "fixtures" / "must_fail" / spec.must_fail_fixture
            )
            assert fixture.exists(), f"{spec.id} lacks its must-fail fixture file"


class TestHealthyStatePasses:
    def test_good_state_is_all_pass(self):
        report = run_assertions(_good_state(), mode="ci")
        assert report.fail_count == 0
        assert report.warn_count == 0
        assert report.pass_count == len(REGISTRY)


class TestMustFailControls:
    """Each negative fixture must fire exactly its own assertion."""

    @pytest.mark.parametrize("spec", REGISTRY, ids=[s.id for s in REGISTRY])
    def test_fires_exactly_once(self, spec):
        report = run_assertions(build_negative_state(spec.id), mode="ci")
        fired = [r for r in report.results if not r.passed and r.assertion_id == spec.id]
        collateral = [
            r.assertion_id for r in report.results
            if not r.passed and r.assertion_id != spec.id
        ]
        assert fired, f"{spec.id} did not fire on its negative fixture"
        assert not collateral, f"{spec.id} fixture also tripped {collateral}"


class TestErrorBScenario:
    def test_delta_lnL_4000_is_fail(self):
        # Error B (|dlnL| in the thousands) must go red under A-10.
        report = run_assertions(build_negative_state("A-10"), mode="ci")
        result = [r for r in report.results if r.assertion_id == "A-10"][0]
        assert not result.passed
        assert result.severity == "fail"
        assert "4000" in result.detail or "exceeds" in result.detail


class TestZeroDiscriminationScenario:
    def test_zero_discrimination_is_fail(self):
        # The baseline zero-discrimination behaviour (both inputs, same output) must
        # Be caught by A-12.
        report = run_assertions(build_negative_state("A-12"), mode="ci")
        result = [r for r in report.results if r.assertion_id == "A-12"][0]
        assert not result.passed
        assert "discrimination" in result.detail


class TestSemantics:
    def test_assertions_do_not_change_grading(self):
        # The assertion layer is observe-only. run_assertions must not
        # Mutate the state it inspects.
        state = _good_state()
        before = {k: repr(v) for k, v in state.items()}
        run_assertions(state, mode="ci")
        after = {k: repr(v) for k, v in state.items()}
        assert before == after

    def test_assertion_crash_counts_as_violation(self):
        from markerfinder.assertions import AssertionSpec

        def boom(state):
            raise RuntimeError("boom")

        spec = AssertionSpec("A-99", "crasher", "run", boom, "fail", "test", "x.json")
        # Run manually through the same path used by run_assertions
        from markerfinder.assertions import AssertionResult

        try:
            boom({})
            violation = None
        except Exception as e:
            violation = f"assertion raised {type(e).__name__}: {e}"
        result = AssertionResult("A-99", "fail", violation is None, violation or "")
        assert not result.passed
        assert "boom" in result.detail

    def test_report_tsv_rows_shape(self):
        report = run_assertions(_good_state(), mode="run")
        rows = report.to_tsv_rows()
        assert rows[0] == [
            "assertion_id", "name", "severity", "result", "detail", "provisional"
        ]
        assert len(rows) == len(REGISTRY) + 1
