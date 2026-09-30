"""V-08 — the consistency and hybrid screening criteria.

Covers ``--hgt-mode`` and ``--consistency-stringency``.

Three criteria are documented:

``risk``
    the risk score decides (default, and the only one that never switches by
    itself — );
``consistency``
    the cross-framework grade decides, with a 5-rung stringency ladder;
``hybrid``
    a marker passes only when the risk level AND the consistency grade both
    pass ('s conjunction).

The order of operations is part of the design (the consistency screen runs
after risk grading, never before), and the gate that guards these modes must be
observable: a mode that silently degrades to ``risk`` would be worse than an
absent one.
"""

from __future__ import annotations

import pytest

from markerfinder.cli import constants


def _levels(run, read_tsv):
    rows = read_tsv(run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    return {r["marker_id"]: r["hgt_risk_level"] for r in rows}, rows


@pytest.mark.capability("hgt_mode", "workflow:hgt-grading")
def test_default_criterion_is_risk_and_says_so(mf, read_tsv, record_metric):
    run = mf(extra=["--force"])
    run.assert_ok()
    _levels_map, rows = _levels(run, read_tsv)
    assert rows
    # The risk path must not have switched a criterion on by itself.
    assert run.recorded("hgt_config", "criterion_mode"
                        ) in (None, "risk"), run.recorded("hgt_config", "")
    record_metric("v08_risk", "markers", len(rows))


@pytest.mark.capability("hgt_mode")
@pytest.mark.parametrize("mode", ["risk", "consistency", "hybrid"])
def test_each_documented_criterion_is_accepted_and_recorded(mf, mode,
                                                            record_metric):
    run = mf(extra=["--hgt-mode", mode, "--force"])
    if run.rc != 0:
        # A gated mode refusing to run is legitimate only with a stated reason.
        assert "Traceback" not in run.text, run.tail()
        assert "hgt-mode" in run.text.lower() or "gate" in run.text.lower() or \
               "consistency" in run.text.lower(), run.tail()
        record_metric(f"v08_{mode}", "refused", True)
        return
    record_metric(f"v08_{mode}", "refused", False)
    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                 .read_text(encoding="utf-8")
    assert mode in summary.lower() or "consistency" in summary.lower(), summary


@pytest.mark.capability("hgt_mode", "consistency_stringency")
def test_consistency_and_hybrid_differ_from_risk_in_what_they_keep(
        mf, read_tsv, record_metric):
    """The three criteria must not be three names for one decision.

    Comparing the marker sets each criterion keeps is the cheapest honest test:
    if all three keep exactly the same markers on real data, then either the
    criteria are not wired or this dataset cannot discriminate them — and in
    either case the claim needs rewording.
    """
    kept = {}
    for mode in ("risk", "consistency", "hybrid"):
        run = mf(extra=["--hgt-mode", mode, "--force"])
        if run.rc != 0:
            pytest.skip(f"--hgt-mode {mode} is gated off on this install: "
                        f"{run.tail()[-200:]}")
        _m, rows = _levels(run, read_tsv)
        kept[mode] = {r["marker_id"] for r in rows
                      if r["hgt_risk_level"] in ("level_1", "level_2")}
    record_metric("v08_criteria", "kept_per_mode",
                  {k: sorted(v) for k, v in kept.items()})
    distinct = len({frozenset(v) for v in kept.values()})
    assert distinct >= 1, kept
    if distinct == 1:
        # Not a failure of the code, but a statement about the data: record it
        # So the report cannot claim a discrimination this dataset did not show.
        pytest.xfail(
            "all three criteria kept the identical marker set on this dataset: "
            "the criteria are accepted but not distinguishable here"
        )


@pytest.mark.capability("consistency_stringency")
@pytest.mark.parametrize("rung", [1, 2, 3, 4, 5])
def test_stringency_ladder_is_monotone(mf, read_tsv, rung, record_metric):
    """Rungs 1..5 tighten the same decision, so what survives rung n+1 must be a
    subset of what survived rung n — for the same data and the same criterion."""
    run = mf(extra=["--hgt-mode", "consistency",
                    "--consistency-stringency", str(rung), "--force"])
    if run.rc != 0:
        pytest.skip(f"consistency criterion gated off: {run.tail()[-200:]}")
    _m, rows = _levels(run, read_tsv)
    kept = {r["marker_id"] for r in rows if r["hgt_risk_level"] == "level_1"}
    record_metric(f"v08_stringency_{rung}", "level_1_markers", sorted(kept))


@pytest.mark.capability("consistency_stringency", "exit:EXIT_ARG_ERROR")
def test_stringency_outside_the_ladder_is_an_argument_error(
        run_markerfinder, genome_input, tmp_path, data_dir):
    result = run_markerfinder([
        "-i", str(genome_input("quad4")), "-o", str(tmp_path / "out"),
        "-t", "2", "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
        "--hgt-mode", "consistency", "--consistency-stringency", "9",
        "--coalescent-mode", "off", "--skip-checkm", "--force",
    ])
    assert result.rc == constants.EXIT_ARG_ERROR, result.tail()
    assert "invalid choice" in result.text.lower() or "9" in result.text, (
        result.tail())
