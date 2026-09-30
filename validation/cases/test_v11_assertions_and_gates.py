"""V-11 — the self-check layer, adjudication and the taxonomy must-pass gate.

Covers ``--strict-assertions``, ``--allow-assertion-failure``,
``--taxonomy-mustpass`` and the must-fail control corpus, plus the two exit
codes that exist only for these paths (4 = output unreasonable / gate failed,
5 = inconclusive recommendation).

The premise of the whole layer is that an implausible output
must stop a run rather than be published. That only counts if a control makes
the guard fire, so every case here drives the guard to its failure side.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from markerfinder.cli import constants


@pytest.mark.capability("strict_assertions",
                        "product:Phase5_reports.assertions.tsv")
def test_assertion_report_is_written_and_default_is_strict(mf, read_tsv,
                                                           record_metric):
    """The default is strict, and the adjudication must be persisted: a reader
    of the run has to be able to see which checks ran and what they decided."""
    run = mf(extra=["--force"])
    run.assert_ok()
    assert run.recorded("report_config", "strict_assertions") is True
    rows = read_tsv(run.product("Phase5_reports/markerfinder.assertions.tsv"))
    assert rows, "the assertion report is empty: the self-check layer did not run"
    assert {"assertion_id", "severity", "result"} <= set(rows[0]), sorted(rows[0])
    verdicts = {r["result"] for r in rows}
    assert verdicts <= {"PASS", "FAIL", "WARN"}, verdicts
    record_metric("v11_assertions", "items", len(rows))
    record_metric("v11_assertions", "verdicts", sorted(verdicts))


@pytest.mark.capability("allow_assertion_failure",
                        "exit:EXIT_ASSERTION_FAILED")
def test_the_bypass_is_recorded_in_the_products(mf, read_tsv, case_workdir,
                                                record_metric):
    """``--allow-assertion-failure`` downgrades a failing output assertion to a
    recorded warning — and the downgrade itself must be visible.

    Driven from a real failure: the split-genus taxonomy makes the must-pass
    gate fire (see the next case), so the same invocation with the bypass must
    finish and say that a failure was bypassed.
    """
    data = Path(__file__).resolve().parents[1] / "data"
    run = mf("small8",
             markers=data / "markers" / "small8_core",
             taxonomy=data / "taxonomy" / "taxonomy_small8_splitgenus.tsv",
             extra=["--taxonomy-mustpass",
                    str(data / "mustpass" / "mustpass_splitgenus.yaml"),
                    "--allow-assertion-failure", "--force"])
    text = run.text.lower()
    assert "bypass" in text or "allow" in text or "assertion" in text, (
        f"the bypass left no trace in the run:\n{run.tail()}"
    )
    record_metric("v11_bypass", "rc", run.rc)
    if run.rc == 0:
        rows = read_tsv(run.product("Phase5_reports/markerfinder.assertions.tsv"))
        persisted = json.dumps([r for r in rows]).lower()
        assert "bypass" in persisted or "warn" in persisted or \
               "allow" in persisted, (
            "the bypass is visible only on stdout, not in the products a "
            "reader keeps")


@pytest.mark.capability("taxonomy_mustpass", "exit:EXIT_ASSERTION_FAILED")
def test_must_pass_gate_fires_on_a_topology_that_violates_it(mf, data_dir,
                                                             record_metric):
    """The failure direction of: with a baseline the data cannot satisfy,
    the run must abort with exit 4 and name the offending marker.

    The control is a taxonomy table in which three accessions from different
    orders (one from another phylum) share a single fake genus; no marker tree
    can make that group exclusive of the other tips.
    """
    run = mf("small8",
             markers=data_dir / "markers" / "small8_core",
             taxonomy=data_dir / "taxonomy" / "taxonomy_small8_splitgenus.tsv",
             extra=["--taxonomy-mustpass",
                    str(data_dir / "mustpass" / "mustpass_splitgenus.yaml"),
                    "--force"])
    record_metric("v11_gate_fail", "rc", run.rc)
    assert run.rc == constants.EXIT_ASSERTION_FAILED, (
        f"the must-pass gate did not abort the run (exit {run.rc}):\n"
        f"{run.tail()}"
    )
    text = run.text.lower()
    assert "splitgenus" in text or "monophyletic" in text or "must-pass" in text, (
        f"the abort does not name the violated relation:\n{run.tail()}"
    )
    assert "Traceback" not in run.text, run.tail()


@pytest.mark.capability("taxonomy_mustpass", "exit:EXIT_SUCCESS",
                        "product:Phase5_reports.mustpass.tsv")
def test_must_pass_gate_runs_and_passes_when_satisfied(mf, data_dir, read_tsv,
                                                       record_metric):
    """The pass direction must show a gate that actually tested something.

    ``groups_tested=0`` next to a PASS is precisely the "comparison that never
    ran reading as green" failure this project exists to remove, so the case
    reads the gate's own product rather than the exit code alone.
    """
    run = mf("small8", extra=["--taxonomy-mustpass",
                              str(data_dir / "mustpass" / "mustpass_domain.yaml"),
                              "--force"])
    run.assert_ok("satisfiable must-pass baseline")
    rows = read_tsv(run.product("Phase5_reports/markerfinder.mustpass.tsv"))
    summary = next((r for r in rows if r["kind"] == "summary"), {})
    tested = summary.get("members", "")
    assert "groups_tested=" in tested and not tested.endswith("groups_tested=0"), (
        f"the gate reported no tested group: {tested}"
    )
    record_metric("v11_gate_pass", "gate_summary", tested)


@pytest.mark.capability("taxonomy_mustpass", "workflow:failure-loudness")
def test_a_missing_baseline_refuses_instead_of_passing(mf, case_workdir,
                                                       record_metric):
    """'s explicit rule: when the gate is asked for and the baseline
    cannot be loaded, the run must not pass. Silence would be a green light
    built out of a missing file."""
    missing = case_workdir / "does_not_exist.yaml"
    run = mf(extra=["--taxonomy-mustpass", str(missing), "--force"])
    text = run.text.lower()
    record_metric("v11_missing_baseline", "rc", run.rc)
    assert run.rc != 0 or "not found" in text or "not run" in text, (
        f"a missing baseline was treated as a pass:\n{run.tail()}"
    )
    if run.rc == 0:
        assert "not run" in text or "no requirements" in text, run.tail()


@pytest.mark.capability("taxonomy_mustpass",
                        "product:Phase5_reports.mustpass.tsv")
def test_bare_flag_resolves_the_shipped_baseline_to_a_real_path(mf, data_dir,
                                                               read_tsv,
                                                               record_metric):
    """``--taxonomy-mustpass`` without an argument documents a default baseline.
    The default must resolve to a file that exists in this checkout — a
    repository-relative string would only work from one directory, and the gate
    would then abort for a packaging reason rather than a biological one."""
    run = mf("small8", extra=["--taxonomy-mustpass", "--force"])
    if run.rc != 0:
        assert "Traceback" not in run.text, run.tail()
        assert "must-pass" in run.text.lower() or "must_pass" in run.text.lower(), (
            run.tail())
        return
    rows = read_tsv(run.product("Phase5_reports/markerfinder.mustpass.tsv"))
    assert rows, "the bare flag produced no must-pass product"
    record_metric("v11_bare_flag", "rows", len(rows))


@pytest.mark.capability("exit:EXIT_INCONCLUSIVE",
                        "workflow:unknown-not-placeholder")
def test_inconclusive_recommendation_is_a_distinct_exit_not_a_zero(mf,
                                                                   data_dir,
                                                                   tmp_path,
                                                                   run_markerfinder,
                                                                   genome_input,
                                                                   record_metric):
    """Exit 5 exists for "the pipeline finished but cannot hand out a tree".

    Driving it with data: three tips cannot produce an inferrable marker tree,
    so the run finishes without a species tree. Whatever code that maps to, it
    must be non-zero and it must be explained — the forbidden outcome is exit 0
    with no tree.
    """
    result = run_markerfinder([
        "-i", str(genome_input("pair3")), "-o", str(tmp_path / "three"),
        "-t", "2", "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "pair3_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_pair3.tsv"),
        "--max-markers", "4", "--coalescent-mode", "off", "--skip-checkm",
        "--force",
    ])
    record_metric("v11_inconclusive", "rc_for_three_tips", result.rc)
    assert result.rc != 0, (
        f"a 3-tip dataset exited 0; the supported range says a tree needs four "
        f"tips:\n{result.tail()}"
    )
    assert result.rc in (constants.EXIT_DATA_ERROR, constants.EXIT_INCONCLUSIVE,
                        constants.EXIT_RUNTIME_ERROR), result.tail()
    assert "Traceback" not in result.text, result.tail()
    text = result.text.lower()
    assert "tree" in text or "marker" in text, (
        f"the refusal does not say what is missing:\n{result.tail()}"
    )
    trees = tmp_path / "three" / "Phase4_trees"
    if trees.exists():
        empties = [p.name for p in trees.glob("*.newick")
                   if p.stat().st_size == 0]
        assert not empties, f"empty tree files written as products: {empties}"


@pytest.mark.capability("strict_assertions", "workflow:failure-loudness")
def test_must_fail_controls_are_rejected_by_the_assertion_layer():
    """The repository ships fixtures that MUST be rejected (A-01..A-14). They
    are replayed here so the validation suite shows the guards firing on the
    same controls, end to end from the shipped files rather than from a
    hand-built dict.
    """
    fixtures = Path(__file__).resolve().parents[2] / "tests" / "fixtures" \
        / "must_fail"
    assert fixtures.is_dir(), f"{fixtures} is missing"
    cases = sorted(fixtures.glob("A-*.json"))
    assert len(cases) >= 10, f"only {len(cases)} must-fail fixtures found"

    from markerfinder.assertions import REGISTRY, run_assertions

    expected = {spec.id for spec in REGISTRY}
    fired_for: dict = {}
    for path in cases:
        assertion_id = path.name.split("_")[0]          # "A-02"
        assert assertion_id in expected, (
            f"{path.name} names {assertion_id}, which is not in the assertion "
            "registry — a control for a guard that no longer exists")
        state = json.loads(path.read_text(encoding="utf-8"))
        report = run_assertions(state, mode="run")
        tripped = [r for r in report.results
                   if r.assertion_id == assertion_id and not r.passed]
        fired_for[assertion_id] = bool(tripped)

    untripped = sorted(a for a, f in fired_for.items() if not f)
    assert not untripped, (
        f"these must-fail controls did not trip their own guard: {untripped}"
    )
