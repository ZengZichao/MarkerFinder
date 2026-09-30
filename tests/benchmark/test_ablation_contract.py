#!/usr/bin/env python3
"""The reverse-ablation comparator refuses unattributable runs.

The reverse-ablation contract is the one requirement that is pure logic, so unlike
the rest of the benchmark it can be proven on this machine without
mafft/IQ-TREE/genomes. What is
proven here is the *refusal* semantics, because that is the requirement:
"同一份数据只改一个因子". A comparator that happily returns a correlation for two
runs differing in three things manufactures exactly the kind of number the
retracted sponge paper's failure mode was built on (Steenwyk & King 2025,
Science, doi:10.1126/science.adw9456; retracted 2026-02-05).

Controls are planted on both sides: identical rankings must measure tau = 1.0,
a single swap must measure less, and a refusal must be a refusal rather than a
plausible coefficient.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1] / "benchmark"
if str(BENCH) not in sys.path:
    sys.path.insert(0, str(BENCH))

import ablation  # Noqa: E402
from ablation import (  # Noqa: E402
    AblationRefusal,
    compare_runs,
    differing_factors,
    kendall_tau_b,
    main,
    rank_by_score,
    read_ranking,
    write_ablation_tsv,
)


# ── the one-factor rule ───────────────────────────────────────────────────

def test_control_differing_factors_counts_both_directions():
    assert differing_factors({"a": "1", "b": "2"}, {"a": "1", "b": "2"}) == []
    assert differing_factors({"a": "1"}, {"a": "2"}) == ["a"]
    # A key present on one side only IS a difference
    assert differing_factors({"a": "1", "x": "on"}, {"a": "2"}) == ["a", "x"]


def test_two_changed_factors_are_refused_by_name():
    with pytest.raises(AblationRefusal) as excinfo:
        compare_runs(
            {"m1": 1, "m2": 2, "m3": 3},
            {"m1": 2, "m2": 1, "m3": 3},
            {"hgt_filter": "on", "concatenation": "yes"},
            {"hgt_filter": "off", "concatenation": "no"},
        )
    message = str(excinfo.value)
    assert "hgt_filter" in message and "concatenation" in message, message
    assert "2 factors differ" in message, message


def test_zero_changed_factors_is_refused_as_a_dead_switch():
    """A comparison that changes zero factors is a dead switch: refuse it."""
    with pytest.raises(AblationRefusal) as excinfo:
        compare_runs(
            {"m1": 1, "m2": 2, "m3": 3},
            {"m1": 1, "m2": 2, "m3": 3},
            {"hgt_filter": "on"},
            {"hgt_filter": "on"},
        )
    assert "nothing to" in str(excinfo.value)


def test_too_few_shared_markers_is_refused():
    with pytest.raises(AblationRefusal) as excinfo:
        compare_runs({"m1": 1, "m2": 2}, {"m1": 2, "m2": 1}, {"f": "a"}, {"f": "b"})
    assert "not a ranking" in str(excinfo.value)


# ── the measurement ───────────────────────────────────────────────────────

def test_control_identical_rankings_measure_perfect_correlation():
    report = compare_runs(
        {"m1": 1, "m2": 2, "m3": 3, "m4": 4},
        {"m1": 1, "m2": 2, "m3": 3, "m4": 4},
        {"hgt_filter": "on"},
        {"hgt_filter": "off"},
    )
    assert report.tau == pytest.approx(1.0)
    assert report.moved == []
    assert report.max_displacement == 0
    assert report.factor == "hgt_filter"


def test_single_adjacent_swap_measures_the_hand_computed_tau():
    """4 markers, one transposed pair: 6 pairs, 1 discordant -> tau = 4/6."""
    report = compare_runs(
        {"m1": 1, "m2": 2, "m3": 3, "m4": 4},
        {"m1": 1, "m2": 2, "m3": 4, "m4": 3},
        {"hgt_filter": "on"},
        {"hgt_filter": "off"},
    )
    assert report.tau == pytest.approx(4 / 6)
    assert [(m, d) for m, _, _, d in report.moved] == [("m3", 1), ("m4", -1)]
    assert report.max_displacement == 1


def test_reverse_ranking_is_negatively_correlated_not_one():
    ranks = {"m1": 1, "m2": 2, "m3": 3, "m4": 4}
    flipped = {"m1": 4, "m2": 3, "m3": 2, "m4": 1}
    report = compare_runs(ranks, flipped, {"f": "a"}, {"f": "b"})
    assert report.tau == pytest.approx(-1.0)
    assert report.max_displacement == 3


def test_all_tied_ranks_are_undefined_rather_than_perfect():
    """A tie on every pair carries no ordering; reporting 1.0 would be a lie."""
    with pytest.raises(ValueError):
        kendall_tau_b([2.0, 2.0, 2.0], [5.0, 5.0, 5.0])


def test_ties_are_shared_not_sequential():
    ranks = rank_by_score({"m1": 90.0, "m2": 90.0, "m3": 50.0})
    assert ranks == {"m1": 1, "m2": 1, "m3": 3}


def test_control_ranking_is_independent_of_input_order():
    """Tie-breaking by dict order would make re-runs disagree."""
    scores = {"b": 1.0, "a": 1.0, "c": 2.0}
    assert rank_by_score(scores) == rank_by_score(dict(reversed(list(scores.items()))))


# ── products and the CLI contract ─────────────────────────────────────────

def _write_tsv(path: Path, rows) -> Path:
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_read_ranking_names_unmeasurable_rows_instead_of_guessing(tmp_path, capsys):
    tsv = _write_tsv(
        tmp_path / "a.tsv",
        ["marker_id\toverall_score", "m1\t0.9", "m2\tNA", "m3\t0.1"],
    )
    scores = read_ranking(tsv, "marker_id", "overall_score")
    assert scores == {"m1": 0.9, "m3": 0.1}
    err = capsys.readouterr().err
    assert "NOT_MEASURABLE" in err and "m2" in err


def test_read_ranking_refuses_a_missing_column(tmp_path):
    tsv = _write_tsv(tmp_path / "a.tsv", ["marker_id\tscore", "m1\t0.5"])
    with pytest.raises(AblationRefusal) as excinfo:
        read_ranking(tsv, "marker_id", "overall_score")
    assert "overall_score" in str(excinfo.value)


def test_read_ranking_refuses_a_file_without_numbers(tmp_path):
    tsv = _write_tsv(tmp_path / "a.tsv", ["marker_id\toverall_score", "m1\tNA"])
    with pytest.raises(AblationRefusal):
        read_ranking(tsv, "marker_id", "overall_score")


def test_written_product_is_deterministic(tmp_path):
    report = compare_runs(
        {"m1": 1, "m2": 2, "m3": 3},
        {"m1": 2, "m2": 1, "m3": 3},
        {"hgt_filter": "on"},
        {"hgt_filter": "off"},
    )
    first = write_ablation_tsv(report, tmp_path, "ab")
    second = write_ablation_tsv(report, tmp_path / "other", "ab")
    assert first.read_bytes() == second.read_bytes()
    lines = first.read_text(encoding="utf-8").strip().split("\n")
    assert lines[0] == "marker_id\trank_baseline\trank_ablated\tdisplacement"
    assert sorted(lines[1:]) == ["m1\t1\t2\t1", "m2\t2\t1\t-1", "m3\t3\t3\t0"]


def test_cli_refuses_to_ablate_without_a_named_factor(tmp_path, capsys):
    a = _write_tsv(tmp_path / "a.tsv", ["marker_id\toverall_score", "m1\t3", "m2\t2", "m3\t1"])
    b = _write_tsv(tmp_path / "b.tsv", ["marker_id\toverall_score", "m1\t1", "m2\t2", "m3\t3"])
    assert main(["--run-a", str(a), "--run-b", str(b)]) == 2
    assert "NOT EXECUTED" in capsys.readouterr().err


def test_cli_reports_a_single_factor_ablation(tmp_path, capsys):
    a = _write_tsv(tmp_path / "a.tsv", ["marker_id\toverall_score", "m1\t3", "m2\t2", "m3\t1"])
    b = _write_tsv(tmp_path / "b.tsv", ["marker_id\toverall_score", "m1\t1", "m2\t2", "m3\t3"])
    assert main([
        "--run-a", str(a), "--run-b", str(b),
        "--factor", "hgt_filter=on:off", "--out-dir", str(tmp_path / "out"),
    ]) == 0
    out = capsys.readouterr().out
    assert "factor hgt_filter" in out
    assert "tau_b=-1.0000" in out
    assert (tmp_path / "out" / "ablation.ablation.tsv").exists()


def test_cli_refuses_two_factors(tmp_path, capsys):
    a = _write_tsv(tmp_path / "a.tsv", ["marker_id\toverall_score", "m1\t3", "m2\t2", "m3\t1"])
    b = _write_tsv(tmp_path / "b.tsv", ["marker_id\toverall_score", "m1\t1", "m2\t2", "m3\t3"])
    assert main([
        "--run-a", str(a), "--run-b", str(b),
        "--factor", "hgt_filter=on:off", "--factor", "tree_method=fasttree:iqtree",
    ]) == 2
    assert "2 factors differ" in capsys.readouterr().err


# ── the benchmark driver exposes the stage (wiring lock) ──────────────────

RANK_A = ["marker_id\toverall_score", "m1\t3", "m2\t2", "m3\t1"]
RANK_B = ["marker_id\toverall_score", "m1\t1", "m2\t2", "m3\t3"]


def test_benchmark_driver_runs_the_ablation_stage_without_external_tools(
    tmp_path, capsys, monkeypatch
):
    """ Must be reachable from run_benchmark, not only from its own CLI.

    The stage deliberately sits before the tool probe: two finished runs are
    enough evidence to compare, and demanding mafft to read two TSVs would push
    users into comparing rankings by hand.

    The tool absence is SIMULATED, not assumed. Reading this test as "the
    machine must have no aligners" made it pass only on a bare environment and
    fail on a fully installed one, where the driver legitimately gets further
    and refuses on the missing downloads instead.
    """
    import run_benchmark

    monkeypatch.setattr(run_benchmark.shutil, "which", lambda _name: None)

    a = _write_tsv(tmp_path / "a.tsv", RANK_A)
    b = _write_tsv(tmp_path / "b.tsv", RANK_B)
    rc = run_benchmark.main([
        "--ablation-a", str(a), "--ablation-b", str(b),
        "--ablation-factor", "hgt_filter=on:off",
    ])
    captured = capsys.readouterr()
    # The stage itself must have produced its verdict...
    assert "factor hgt_filter" in captured.out, captured.out
    assert "tau_b=-1.0000" in captured.out, captured.out
    #...and the driver then still refuses the rest of the benchmark, because the
    # External tools are absent. rc 2 here is the CORRECT outcome: it proves the
    # Ablation ran first and that nothing downstream was allowed to pass on
    # Borrowed credit.
    assert rc == 2
    assert "external tools missing" in captured.err, captured.err


def test_benchmark_driver_refuses_a_one_sided_ablation(tmp_path, capsys):
    import run_benchmark

    a = _write_tsv(tmp_path / "a.tsv", RANK_A)
    assert run_benchmark.main(["--ablation-a", str(a)]) == 2
    assert "--ablation-b" in capsys.readouterr().err


def test_benchmark_driver_propagates_the_two_factor_refusal(tmp_path, capsys):
    import run_benchmark

    a = _write_tsv(tmp_path / "a.tsv", RANK_A)
    b = _write_tsv(tmp_path / "b.tsv", RANK_B)
    rc = run_benchmark.main([
        "--ablation-a", str(a), "--ablation-b", str(b),
        "--ablation-factor", "hgt_filter=on:off",
        "--ablation-factor", "tree_method=fasttree:iqtree",
    ])
    assert rc == 2
    assert "2 factors differ" in capsys.readouterr().err
