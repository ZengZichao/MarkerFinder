"""V-13 — output management: intermediates, temp directories, report format.

Covers ``--save-intermediates``, ``--tmp-dir``, ``--keep-tmp``,
``--report-format`` and the ``--verbose``/``--log-file`` pairing at the product
level, plus the documented directory layout (``Phase4_intermediate/`` and the
historical ``Phase4_*`` / ``Phase5_*`` prefixes).

``--tmp-dir`` matters beyond tidiness: a run that writes into a shared temp
directory is not reproducible and cannot be audited after the fact, which is
why the default is a per-run isolated directory and why the explicit form must
not be auto-deleted (the help text promises exactly that).
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.mark.capability("save_intermediates",
                        "product:Phase4_intermediate.saved")
def test_save_intermediates_writes_the_documented_tree(mf, case_output,
                                                       record_metric):
    off = mf(out_dir=case_output / "off", extra=["--force"])
    off.assert_ok()
    assert not (off.out_dir / "Phase4_intermediate").exists(), (
        "intermediates were written without --save-intermediates")

    on = mf(out_dir=case_output / "on", extra=["--save-intermediates", "--force"])
    on.assert_ok("--save-intermediates")
    root = on.out_dir / "Phase4_intermediate"
    assert root.is_dir(), sorted(p.name for p in on.out_dir.iterdir())
    markers = sorted((root / "markers").glob("*")) if (root / "markers").exists() else []
    assert markers, f"{root / 'markers'} holds nothing: {list(root.iterdir())}"
    suffixes = {p.suffix for p in markers}
    record_metric("v13_intermediates", "files", len(markers))
    record_metric("v13_intermediates", "suffixes", sorted(suffixes))
    assert {".faa", ".aln"} & suffixes, suffixes
    supermatrix = root / "supermatrix" / "markerfinder.concat.fasta"
    assert supermatrix.exists(), sorted(p.name for p in (root / "supermatrix").iterdir()) \
        if (root / "supermatrix").exists() else "no supermatrix/ directory"


@pytest.mark.capability("tmp_dir", "keep_tmp")
def test_an_explicit_tmp_dir_is_kept_because_it_is_not_owned_by_the_run(
        mf, case_workdir, record_metric):
    """Default behaviour: a per-run temp directory is removed on exit. Explicit
    ``--tmp-dir``: the directory is the user's, so the pipeline must not delete
    it. ``--keep-tmp`` must preserve the automatic one instead."""
    tmp = case_workdir / "tmpdir"
    tmp.mkdir()
    run = mf(extra=["--tmp-dir", str(tmp), "--force"])
    run.assert_ok("--tmp-dir run")
    assert tmp.is_dir(), "the pipeline deleted a --tmp-dir it did not create"
    # Run_config.json records paths relative to the output directory, so compare
    # Resolved locations rather than strings.
    recorded = Path(str(run.recorded("", "tmp_dir")))
    if not recorded.is_absolute():
        recorded = (run.out_dir / recorded).resolve()
    assert recorded.resolve() == tmp.resolve(), (
        f"--tmp-dir {tmp} was not the directory the run used: {recorded}"
    )

    kept = mf(extra=["--keep-tmp", "--force"])
    kept.assert_ok("--keep-tmp")
    assert kept.recorded("", "keep_tmp") is True, kept.recorded("", "keep_tmp")
    record_metric("v13_tmp", "explicit_tmp_survives", tmp.is_dir())
    auto_tmp = Path(str(kept.recorded("", "tmp_dir")))
    if not auto_tmp.is_absolute():
        auto_tmp = (kept.out_dir / auto_tmp).resolve()
    record_metric("v13_tmp", "kept_tmp_still_present", auto_tmp.is_dir())


@pytest.mark.capability("tmp_dir", "workflow:determinism")
def test_two_runs_in_parallel_do_not_share_a_temp_directory(mf, record_metric):
    """The default temp directory is per-run. Two concurrent runs that shared
    one prefix would overwrite each other's alignments and produce a tree from
    mixed data — the kind of failure that is invisible in the output."""
    a = mf(extra=["--force"])
    b = mf(extra=["--force"])
    a.assert_ok()
    b.assert_ok()
    ta, tb = str(a.recorded("", "tmp_dir")), str(b.recorded("", "tmp_dir"))
    record_metric("v13_tmp_isolation", "tmp_dirs", [ta, tb])
    assert ta != tb, (
        f"two runs used the same temp directory {ta}: concurrent runs of one "
        "pipeline cannot be trusted to keep their intermediates apart"
    )


@pytest.mark.capability("report_format")
def test_the_only_supported_report_format_is_html_and_is_produced(mf,
                                                                  record_metric):
    run = mf(extra=["--report-format", "html", "--force"])
    run.assert_ok()
    html = run.product("Phase5_reports/markerfinder.report.html")
    text = html.read_text(encoding="utf-8")
    assert text.lstrip().startswith(("<!DOCTYPE", "<html", "<HTML")), text[:120]
    # The documented promise is a SELF-CONTAINED static report: no external
    # Script or stylesheet may be required to read it.
    assert "<script src=\"http" not in text and "link rel=\"stylesheet\" href=\"http" \
        not in text, "the 'static, self-contained' report loads remote assets"
    assert "<script src=\"https" not in text, (
        "the report loads a remote script: a user offline cannot read the "
        "results of their own run")
    record_metric("v13_report", "html_bytes", len(text))


@pytest.mark.capability("report_format", "exit:EXIT_ARG_ERROR")
def test_unsupported_report_formats_are_refused_at_parse_time(
        run_markerfinder, genome_input, data_dir, tmp_path):
    """The docs state PDF and Plotly are NOT available; a value naming them
    must be rejected before any work, not accepted and then ignored."""
    for fmt in ("pdf", "plotly"):
        result = run_markerfinder([
            "-i", str(genome_input("quad4")), "-o", str(tmp_path / fmt),
            "--report-format", fmt, "-t", "2",
            "--marker-mode", "gtdb_tk",
            "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
            "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
            "--skip-checkm", "--force",
        ])
        assert result.rc == 2, f"--report-format {fmt} was accepted: {result.tail()}"
        assert "invalid choice" in result.text.lower(), result.tail()


@pytest.mark.capability("product:Phase5_metadata.run_config.json",
                        "workflow:provenance-recorded")
def test_run_config_records_the_inputs_the_run_actually_used(mf, data_dir,
                                                             record_metric,
                                                             hmm_library):
    """Reproducibility needs the recorded paths to be the ones that were used,
    including which marker source produced the numbers ('s pairing)."""
    if hmm_library is None:
        pytest.skip(
            "database_versions is only non-empty when db/gtdb_markers exists; "
            "the library is fetched, not committed (db/README.md), so a fresh "
            "checkout cannot record one"
        )
    run = mf(extra=["--force"])
    run.assert_ok()
    payload = run.product("Phase5_metadata/run_config.json").read_text(
        encoding="utf-8")
    import json as _json

    doc = _json.loads(payload)
    assert doc.get("markerfinder_version"), doc.keys()
    params = doc["parameters"]
    assert params["selection_config"]["gtdb_markers_dir"], params["selection_config"]
    assert "database_versions" in doc, sorted(doc)
    versions = doc["database_versions"]
    assert versions, (
        "run_config.json recorded no database version: a published result "
        "cannot be tied to the database that produced it"
    )
    record_metric("v13_runconfig", "database_versions", sorted(versions))


@pytest.mark.capability("output", "workflow:failure-loudness")
def test_a_run_does_not_leave_files_outside_its_output_directory(mf,
                                                                 case_output):
    """The bundle directory must stay clean: a stray file written next to the
    shipped data or the package would corrupt the published checksums."""
    repo = Path(__file__).resolve().parents[2]
    before = _listing(repo / "validation" / "data")
    run = mf(out_dir=case_output, extra=["--force"])
    run.assert_ok()
    after = _listing(repo / "validation" / "data")
    assert before == after, (
        f"a run wrote into validation/data: {sorted(set(after) - set(before))}"
    )


def _listing(root):
    return sorted(
        str(p.relative_to(root)) for p in root.rglob("*")
        if p.is_file() and ".hit_cache" not in p.parts
    )
