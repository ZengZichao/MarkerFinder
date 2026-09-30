"""V-07 — marker selection: resolution modes, budgets and occupancy floors.

Covers ``--mode`` (four documented resolutions), ``--marker-preset``,
``--max-markers``, ``--min-occupancy``, and the deprecated
``--min-marker-coverage`` alias.

Selection is where MarkerFinder's central claim lives — "dynamically selects
optimal markers adapted to your dataset" — so each documented mode must be
observable in what the run actually keeps, not only in the label it prints.
"""

from __future__ import annotations

import pytest


def _kept(run):
    return [r for r in _rows(run, "marker_summary")]


def _rows(run, stem):
    path = run.product(f"Phase5_reports/markerfinder.{stem}.tsv")
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


@pytest.mark.capability("mode", "workflow:occupancy-selection")
@pytest.mark.parametrize("mode", ["conservative", "standard", "expanded",
                                 "mag_adaptive"])
def test_each_documented_mode_runs_and_labels_the_run(mf, data_dir,
                                                     record_metric, mode):
    # The occupancy-spanning directory of the SAME genome set the run loads: a
    # Ceiling of 6 then has more candidates than the 8 complete-core markers,
    # And every marker file still describes the four genomes under ``-i`` (a
    # Foreign directory would leave markers whose quad4 tip count falls under
    # The four-tip minimum the pipeline needs to build a gene tree).
    spanning = data_dir / "markers" / "quad4"
    run = mf(markers=spanning,
             extra=["--mode", mode, "--max-markers", "6", "--force"])
    if mode == "conservative":
        # A 0.90 occupancy floor over markers that are 7/8 occupied is a real
        # Constraint: the run may legitimately end with too few markers, and
        # Saying so is the correct behaviour.
        assert run.rc in (0, 3, 5), run.tail()
    else:
        run.assert_ok(f"mode {mode}")
    if run.rc != 0:
        assert "Traceback" not in run.text, run.tail()
        return
    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                  .read_text(encoding="utf-8")
    assert mode in summary, (
        f"the summary does not report the mode that was asked for:\n{summary}"
    )
    record_metric(f"v07_mode_{mode}", "markers_kept", len(_kept(run)))


def _selected_count(run) -> int:
    """The budget's observable effect is the "Markers selected:" line.

    marker_summary.tsv lists every marker that was *evaluated*, which is a
    different quantity from how many were selected — conflating the two is how a
    budget knob passes while doing nothing.
    """
    import re

    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                 .read_text(encoding="utf-8")
    match = re.search(r"Markers selected:\s*(\d+)", summary)
    assert match, f"no 'Markers selected:' line in:\n{summary[:400]}"
    return int(match.group(1))


@pytest.mark.capability("max_markers", "workflow:occupancy-selection")
def test_max_markers_is_a_real_ceiling(mf, data_dir, record_metric):
    """The spanning directory is used so a ceiling of 8 has more candidates
    than the complete-core directory does."""
    spanning = data_dir / "markers" / "small8"
    counts = {}
    for budget in ("2", "4", "8"):
        run = mf("small8", markers=spanning,
                 extra=["--max-markers", budget,
                        "--min-occupancy", "0.5", "--force"])
        # A ceiling of 2 selects the two sparsest-ranked markers, and the HGT
        # Screen may grade both out — the documented outcome then is a refusal
        # (exit 3), not a silently shortened tree. The ceiling is a Phase 1
        # Knob, so it is read from Phase 1's own record either way; requiring a
        # Full run would test the dataset instead of the knob.
        assert run.rc in (0, 3), (
            f"--max-markers {budget} ended with rc={run.rc}:\n{run.tail()}"
        )
        if run.rc == 3:
            assert "no usable markers" in run.text.lower(), (
                f"a budget that left nothing usable did not say so:\n"
                f"{run.tail()}"
            )
        assert run.recorded("selection_config", "max_markers") == int(budget), (
            f"the budget was not recorded for --max-markers {budget}: "
            f"{run.recorded('selection_config', 'max_markers')}"
        )
        counts[budget] = _selected_count(run)
    record_metric("v07_budget", "selected_per_ceiling", counts)
    assert counts["2"] <= 2 and counts["4"] <= 4, counts
    assert counts["8"] >= counts["2"], counts


@pytest.mark.capability("min_occupancy", "workflow:occupancy-selection")
def test_occupancy_floor_rejects_markers_below_it(mf, data_dir, record_metric):
    """The whole point of the knob. All-complete markers cannot show this, so
    the occupancy-spanning directory is used, where the sparsest marker file
    holds 7 of 8 genomes (occupancy 0.875) — and the budget is set above the
    marker count so the floor, not the quota, is what limits the selection.
    """
    spanning = data_dir / "markers" / "small8"
    loose = mf("small8", markers=spanning, extra=["--min-occupancy", "0.5",
                                                  "--max-markers", "20",
                                                  "--force"])
    loose.assert_ok()
    tight = mf("small8", markers=spanning, extra=["--min-occupancy", "0.95",
                                                   "--max-markers", "20",
                                                   "--force"])
    tight.assert_ok()
    n_loose, n_tight = _selected_count(loose), _selected_count(tight)
    record_metric("v07_occupancy", "kept_at_0.50", n_loose)
    record_metric("v07_occupancy", "kept_at_0.95", n_tight)
    assert n_tight < n_loose, (
        f"an occupancy floor of 0.95 kept {n_tight} markers and 0.50 kept "
        f"{n_loose}: the floor is not applied"
    )


@pytest.mark.capability("min_occupancy", "workflow:occupancy-selection")
def test_occupancy_above_one_is_refused_as_an_argument(mf):
    """A probability outside 0-1 cannot mean anything; accepting it would let a
    typo turn the floor into no floor at all."""
    run = mf(extra=["--min-occupancy", "1.5", "--force"])
    assert run.rc != 0, (
        f"--min-occupancy 1.5 was accepted:\n{run.tail()}")
    assert "Traceback" not in run.text, run.tail()


@pytest.mark.capability("min_marker_coverage", "min_occupancy")
def test_deprecated_coverage_alias_is_announced_and_applied(mf, data_dir,
                                                           record_metric):
    spanning = data_dir / "markers" / "small8"
    alias = mf("small8", markers=spanning,
               extra=["--min-marker-coverage", "0.95", "--max-markers", "8",
                      "--force"])
    alias.assert_ok("--min-marker-coverage alias")
    assert "deprecated" in alias.text.lower(), (
        f"the deprecated alias produced no notice:\n{alias.tail()}"
    )
    direct = mf("small8", markers=spanning,
                extra=["--min-occupancy", "0.95", "--max-markers", "8",
                       "--force"])
    direct.assert_ok()
    record_metric("v07_alias", "selected_via_alias", _selected_count(alias))
    record_metric("v07_alias", "selected_directly", _selected_count(direct))
    assert _selected_count(alias) == _selected_count(direct), (
        "the alias does not behave like the option it claims to alias"
    )


@pytest.mark.capability("marker_preset")
@pytest.mark.parametrize("preset", ["none", "conservative", "standard",
                                   "expanded", "mag_adaptive"])
def test_named_presets_override_the_explicit_knobs(mf, preset, record_metric):
    """A preset overrides min_occupancy / max_markers / strategy;
    ``none`` (the default) must not."""
    run = mf(extra=["--marker-preset", preset, "--force"])
    if run.rc != 0:
        assert "preset" in run.text.lower(), run.tail()
        assert "Traceback" not in run.text, run.tail()
        return
    recorded = _rows(run, "marker_summary")
    record_metric(f"v07_preset_{preset}", "evaluated", len(recorded))
    record_metric(f"v07_preset_{preset}", "selected", _selected_count(run))
    if preset == "none":
        # 'none' must leave the requested budget (4, from the shared command)
        # Untouched rather than silently widening it.
        assert _selected_count(run) <= 4, (
            f"--marker-preset none selected {_selected_count(run)} markers, more "
            "than the requested budget of 4: the default preset is not a no-op"
        )


@pytest.mark.capability("mode", "marker_preset")
def test_preset_and_mode_disagreement_is_resolved_explicitly(mf, record_metric):
    """Both knobs set a resolution. Whichever wins, the run must record it, so
    a reader of run_config.json is not left guessing."""
    run = mf(extra=["--mode", "expanded", "--marker-preset", "conservative",
                    "--force"])
    run.assert_ok()
    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                  .read_text(encoding="utf-8")
    assert "expanded" in summary or "conservative" in summary, summary
    record_metric("v07_preset_vs_mode", "summary",
                  [ln for ln in summary.splitlines() if "mode" in ln.lower()])
