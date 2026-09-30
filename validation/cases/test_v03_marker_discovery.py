"""V-03 — marker discovery: the ``hmm`` and ``gtdb_tk`` modes.

Covers ``--marker-mode``, ``--marker-hmm-dir``, ``--marker-db-source``,
``--min-hmm-score``, ``--gtdb-markers-dir`` and ``--db-dir``.

Both documented discovery routes are run end to end. In ``hmm`` mode the
software identifies markers itself by scanning each input proteome against
TIGRFAM/Pfam profiles; in ``gtdb_tk`` mode it consumes a per-marker FASTA
directory produced upstream. The two routes must produce the same KIND of
evidence (occupancy, marker summary, HGT grading), and each must say which
marker source it used in ``run_config.json``.
"""

from __future__ import annotations

import pytest

from markerfinder.cli import constants


def _summary(run):
    return run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
               .read_text(encoding="utf-8")


@pytest.mark.capability("marker_mode", "marker_hmm_dir", "min_hmm_score",
                        "product:Phase5_reports.marker_summary.tsv",
                        "exit:EXIT_SUCCESS")
def test_hmm_mode_runs_end_to_end_and_records_its_profile_source(
        mf, data_dir, hmm_dir):
    run = mf(extra=["--marker-mode", "hmm", "--marker-hmm-dir",
                    str(hmm_dir("core")), "--min-hmm-score", "20",
                    "--force"])
    run.assert_ok("hmm mode")
    assert run.recorded("selection_config", "marker_mode") == "hmm", \
        run.recorded("selection_config", "")
    assert "core" in str(run.recorded("selection_config", "marker_hmm_dir")), \
        run.recorded("selection_config", "marker_hmm_dir")
    assert run.recorded("selection_config", "min_hmm_score") == 20.0
    summary = _summary(run)
    assert "Markers selected:" in summary, summary
    assert (run.out_dir / "Phase5_reports"
                       / "markerfinder.marker_summary.tsv").exists()


@pytest.mark.capability("min_hmm_score", "marker_mode")
def test_raising_the_bitscore_floor_cannot_add_markers(mf, hmm_dir,
                                                       record_metric):
    """A stricter hit threshold is a stricter filter, never a wider one.

    The assertion is differential (two runs, one knob), which is what makes it
    a test of the knob rather than of the dataset.
    """
    loose = mf(extra=["--marker-mode", "hmm", "--marker-hmm-dir",
                      str(hmm_dir("core")), "--min-hmm-score", "5.0",
                      "--max-markers", "8", "--force"])
    loose.assert_ok("hmm mode, loose bitscore floor")
    tight = mf(extra=["--marker-mode", "hmm", "--marker-hmm-dir",
                      str(hmm_dir("core")), "--min-hmm-score", "2000.0",
                      "--max-markers", "8", "--force"])
    n_loose = len(_summary_rows(loose))
    n_tight = len(_summary_rows(tight)) if tight.out_dir else 0
    record_metric("v03_bitscore", "markers_at_score_5", n_loose)
    record_metric("v03_bitscore", "markers_at_score_2000", n_tight)
    assert n_loose >= n_tight, (
        f"--min-hmm-score 2000 selected {n_tight} markers while 5.0 selected "
        f"{n_loose}: the floor does not filter"
    )


def _summary_rows(run):
    path = run.out_dir / "Phase5_reports" / "markerfinder.marker_summary.tsv"
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    return lines[1:] if len(lines) > 1 else []


@pytest.mark.capability("marker_db_source", "db_dir", "marker_mode",
                        "workflow:provenance-recorded")
def test_hmm_mode_auto_discovers_the_bundled_library(mf, record_metric):
    """Omitting ``--marker-hmm-dir`` is the documented no-GTDB-Tk route: the
    bundled TIGRFAM/Pfam library under ``--db-dir`` must be found by domain
    (``--marker-db-source auto``, the default) and the directory actually used
    must be recorded, because a run that silently scanned nothing would still
    print a report.``
    """
    run = mf(omit=("--gtdb-markers-dir",),
             extra=["--marker-mode", "hmm", "--db-dir", "db", "--force"])
    run.assert_ok(
        "hmm mode with auto-discovery. If this fails because discovery is not "
        "reachable when --marker-hmm-dir is absent, the claim in the README "
        "('when omitted the bundled library under --db-dir is auto-discovered') "
        "is false and the documentation must change"
    )
    assert run.recorded("selection_config", "marker_mode") == "hmm"
    discovered = str(run.recorded("selection_config", "marker_hmm_dir"))
    assert discovered.strip(), (
        "the run recorded an empty marker_hmm_dir: auto-discovery found a "
        "library but did not say which one, so the result cannot be reproduced"
    )
    assert "Markers selected:" in _summary(run), _summary(run)
    record_metric("v03_autodiscovery", "discovered_dir", discovered)


@pytest.mark.capability("marker_db_source")
def test_db_source_none_requires_an_explicit_profile_directory(
        run_markerfinder, genome_input, repo_root, tmp_path):
    """``--marker-db-source none`` disables auto-discovery, so hmm mode without
    ``--marker-hmm-dir`` must be a named configuration error, not a run with
    zero markers."""
    result = run_markerfinder([
        "-i", str(genome_input("pair3")), "-o", str(tmp_path / "out"),
        "-t", "2", "--marker-mode", "hmm", "--marker-db-source", "none",
        "--skip-checkm", "--force",
    ])
    assert result.rc != 0, result.tail()
    assert "--marker-hmm-dir" in result.text, (
        f"the refusal does not name the missing option:\n{result.tail()}")


@pytest.mark.capability("marker_db_source")
def test_db_source_cog_is_reported_not_silently_ignored(run_markerfinder,
                                                        genome_input, tmp_path):
    """``cog`` names a library layout this checkout does not ship. Accepting
    the value and then scanning nothing would be the silent-zero failure mode."""
    result = run_markerfinder([
        "-i", str(genome_input("pair3")), "-o", str(tmp_path / "cog"),
        "-t", "2", "--marker-mode", "hmm", "--marker-db-source", "cog",
        "--db-dir", "db", "--skip-checkm", "--force",
    ])
    text = result.text
    if result.rc == 0:
        assert "cog" in text.lower(), (
            "the run claims to have used the cog library but never says so")
    else:
        assert "cog" in text.lower() or "hmm" in text.lower(), result.tail()
        assert "Traceback" not in text, result.tail()


@pytest.mark.capability("gtdb_markers_dir", "marker_mode")
def test_gtdb_tk_mode_without_a_marker_directory_is_a_named_refusal(
        run_markerfinder, genome_input, tmp_path):
    result = run_markerfinder([
        "-i", str(genome_input("pair3")), "-o", str(tmp_path / "nout"),
        "-t", "2", "--marker-mode", "gtdb_tk", "--skip-checkm", "--force",
    ])
    assert result.rc in (constants.EXIT_ARG_ERROR, constants.EXIT_DATA_ERROR,
                         constants.EXIT_RUNTIME_ERROR), result.tail()
    assert "gtdb" in result.text.lower(), result.tail()
    assert "Traceback" not in result.text, result.tail()


@pytest.mark.capability("gtdb_markers_dir", "marker_mode")
def test_gtdb_tk_mode_rejects_a_directory_with_no_marker_files(
        run_markerfinder, genome_input, tmp_path, empty_input_dir):
    empty = empty_input_dir("empty_markers")
    result = run_markerfinder([
        "-i", str(genome_input("pair3")), "-o", str(tmp_path / "out"),
        "-t", "2", "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(empty), "--skip-checkm", "--force",
    ])
    assert result.rc != 0, (
        "an empty marker directory produced a 'successful' run with no markers")
    assert "marker" in result.text.lower(), result.tail()


@pytest.mark.capability("marker_mode", "marker_hmm_dir", "gtdb_markers_dir")
def test_both_modes_leave_the_mode_they_ran_in_the_record(mf, data_dir,
                                                          hmm_dir, marker_set,
                                                          taxonomy_table):
    """``run_config.json`` is the reproducibility record: a reader must be able
    to tell which discovery route produced the numbers."""
    hmm_run = mf(extra=["--marker-mode", "hmm", "--marker-hmm-dir",
                        str(hmm_dir("core")), "--force"])
    hmm_run.assert_ok("hmm route")
    gtdb_run = mf(markers=marker_set("small8"),
                  taxonomy=taxonomy_table("small8"),
                  extra=["--marker-mode", "gtdb_tk", "--force"])
    gtdb_run.assert_ok("gtdb_tk route")
    for run, expected in ((hmm_run, "hmm"), (gtdb_run, "gtdb_tk")):
        assert run.recorded("selection_config", "marker_mode") == expected, (
            expected, run.recorded("selection_config", "marker_mode")
        )
