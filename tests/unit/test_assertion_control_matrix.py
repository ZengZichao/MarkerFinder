""" Control experiment: does --check's verdict come from the check bodies?

The acceptance list demands "delete any assertion's check body and --check must
go red". A one-off manual experiment proves nothing tomorrow, so it is encoded
here as a repeatable control matrix: for EVERY registered assertion, its
negative fixture must (a) fire while the check body is intact, and (b) stop
firing once the body is replaced by a no-op.

Case (b) is the point. A fixture that still "passes the check" against a neutered
body would mean the verdict is hard-coded somewhere upstream of the check —
exactly the failure mode that let the retracted paper ship 785 green tests.
"""

from __future__ import annotations

import dataclasses
import pytest

import markerfinder.assertions as assertions_module
from markerfinder.assertions import REGISTRY, run_assertions


def _fired(report, assertion_id: str) -> bool:
    """A negative fixture 'fires' when that assertion reports not-passed."""
    return any(
        r.assertion_id == assertion_id and r.passed is False
        for r in report.results
    )


def test_every_registered_assertion_has_a_negative_fixture_path():
    assert len(REGISTRY) >= 14, "G2 requires at least 14 registered assertions"
    for spec in REGISTRY:
        assert spec.must_fail_fixture, f"{spec.id} has no must-fail fixture"


@pytest.fixture(autouse=True)
def _registry_is_restored():
    """Never let a neutered assertion leak into another test."""
    snapshot = list(assertions_module.REGISTRY)
    yield
    assertions_module.REGISTRY[:] = snapshot


@pytest.mark.parametrize("index", range(len(REGISTRY)), ids=[
    spec.id for spec in REGISTRY
])
# (all fourteen assertions implemented) and (they must take effect
# At run end, in --check, and in CI): this matrix is the CI carrier, and it
# Exercises each check body rather than its mere registration.
def test_control_matrix_neutering_the_check_body_turns_check_green(index):
    spec = assertions_module.REGISTRY[index]
    state = assertions_module.build_negative_state(spec.id)

    # (a) With the real body, the fixture must trip THIS assertion.
    intact = run_assertions(state, mode="check")
    assert _fired(intact, spec.id), (
        f"{spec.id}: its own negative fixture does not fire the check — the "
        f"assertion is untestable as written"
    )

    # (b) Replace the body with a no-op: the verdict must flip. If it stays red,
    # The red came from somewhere other than this check body.
    try:
        assertions_module.REGISTRY[index] = dataclasses.replace(
            spec, check=lambda state, *a, **k: None
        )
        after = run_assertions(
            assertions_module.build_negative_state(spec.id), mode="check"
        )
    finally:
        assertions_module.REGISTRY[index] = spec
    assert not _fired(after, spec.id), (
        f"{spec.id} still reports failure with an EMPTY check body: the "
        f"--check verdict does not depend on this assertion's logic"
    )


def test_control_matrix_actually_ran_on_all_assertions():
    # Guard against a silently empty parametrisation, which would make the
    # Control above report green for zero cases.
    assert len(REGISTRY) == len({spec.id for spec in REGISTRY})
    collected = [
        spec.id for spec in REGISTRY
    ]
    assert len(collected) >= 14
    assert collected[0] == "A-01" and collected[-1] == "A-14"


def test_fr30_assertion_report_is_persisted(tmp_path):
    """The verdicts must exist as a file a reviewer can read."""
    from markerfinder.assertions import AssertionReport, AssertionResult
    from markerfinder.modules.report_generator import write_assertions_tsv

    report = AssertionReport(
        results=[
            AssertionResult(assertion_id="A-02", severity="fail", passed=False,
                            detail="normalized_rf=1.7 out of [0,1]"),
            AssertionResult(assertion_id="A-07", severity="warn", passed=True,
                            detail="empty alignment absent"),
        ],
        mode="run",
    )
    path = tmp_path / "mf.assertions.tsv"
    write_assertions_tsv(report, path)
    body = path.read_text(encoding="utf-8")
    header = body.splitlines()[0].split("\t")
    for column in ("assertion_id", "severity", "result", "detail"):
        assert column in header, header
    assert "A-02" in body and "A-07" in body


def test_missing_assertion_product_is_not_swallowed_at_debug_level():
    """A run with no audit product must say so where a reader will see it."""
    import inspect

    import markerfinder.pipeline as pipeline_module

    source = inspect.getsource(pipeline_module)
    guard = 'assertions.tsv NOT WRITTEN'
    assert guard in source, (
        "the assertions.tsv failure was downgraded again: the product could "
        "go missing with only a debug line to show for it"
    )
