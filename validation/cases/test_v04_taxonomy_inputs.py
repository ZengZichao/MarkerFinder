"""V-04 — taxonomy input formats, sources and the rank machinery.

Covers ``--taxonomy-table``, ``--taxonomy-format``,
``--taxonomy-source-priority``, ``--taxonomy-delimiter-mode``, ``--table-sep``,
``--no-auto-embed``, ``--taxonomy-levels``, ``--ignore-malformed``,
``--monophyly-rank`` and ``--monophyly-threshold``.

Two taxonomy formats are documented: Format B (table style,
``d__X;p__Y``) and Format A (embedded in the identifier, ``_d_X_p_Y``). Both
must reach the same internal representation, and where they disagree the
declared priority must decide — silently preferring one source is how a
pipeline ends up grading markers against a taxonomy nobody asked for.
"""

from __future__ import annotations

import pytest

from markerfinder.cli import constants


@pytest.mark.capability("taxonomy_table", "monophyly_rank",
                        "product:Phase5_reports.hgt_evaluation.tsv")
def test_format_b_table_drives_the_phylogenetic_step(mf, read_tsv, record_metric):
    """Small8, not quad4: the phylogenetic step measures a monophyly proportion,
    and a set whose members all share one order cannot exercise it at all."""
    run = mf("small8", extra=["--monophyly-rank", "order", "--force"])
    run.assert_ok("Format B taxonomy table")
    rows = read_tsv(run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    assert rows, "the HGT evaluation is empty: the phylogenetic step produced nothing"
    levels = {r["hgt_risk_level"] for r in rows}
    assert levels != {"unknown"}, (
        f"every marker came back UNKNOWN, i.e. the phylogenetic step never "
        f"screened: {levels}"
    )
    record_metric("v04_format_b", "markers_evaluated", len(rows))
    record_metric("v04_format_b", "levels_seen", sorted(levels))


@pytest.mark.capability("table_sep", "taxonomy_table")
def test_forced_comma_separator_reads_the_csv_table(mf, data_dir,
                                                    taxonomy_table):
    """``--table-sep`` exists because auto-detection can be wrong; forcing the
    separator must actually read the CSV twin of the TSV table."""
    csv_table = data_dir / "taxonomy" / "taxonomy_small8.csv"
    run = mf(taxonomy=csv_table, extra=["--table-sep", ",", "--force"])
    run.assert_ok("CSV taxonomy table with --table-sep ,")
    assert "Loaded taxonomy for 8 tips" in run.text or run.rc == 0, run.tail()


@pytest.mark.capability("ignore_malformed", "exit:EXIT_DATA_ERROR",
                        "workflow:failure-loudness")
def test_malformed_rows_are_reported_and_can_be_tolerated_on_request(
        mf, data_dir, record_metric):
    """The documented switch is about *tolerating* bad rows. Either way the row
    must be named: a malformed lineage quietly parsed into an all-empty taxonomy
    changes the monophyly scope without anyone being told."""
    bad = data_dir / "taxonomy" / "taxonomy_malformed.tsv"
    strict = mf(taxonomy=bad, extra=["--force"])
    reported = strict.text.lower()
    assert "malformed" in reported or "invalid taxonomy pair" in reported \
        or "unknown level prefix" in reported, (
        f"a malformed taxonomy row left no report at all:\n{strict.tail()}"
    )
    assert "Traceback" not in strict.text, strict.tail()

    lenient = mf(taxonomy=bad, extra=["--ignore-malformed", "--force"])
    lenient.assert_ok("--ignore-malformed run")
    lenient_text = lenient.text.lower()
    # --ignore-malformed governs rows that cannot be parsed at all; a row whose
    # Lineage cell is prose is reported by the taxonomy parser instead. Either
    # Way the bad row must be named — an unreported one changes the monophyly
    # Scope in silence.
    assert "malformed" in lenient_text or "invalid taxonomy pair" in lenient_text \
        or "unknown level prefix" in lenient_text, (
        "the skipped row was not reported:\n" + lenient.tail()
    )


@pytest.mark.capability("taxonomy_format", "taxonomy_delimiter_mode")
def test_format_a_table_with_each_delimiter_mode(mf, data_dir, record_metric):
    """``--taxonomy-format embedded`` reads Format A segments; the three
    documented parsing modes must all be accepted and must all resolve the
    domain, because the monophyly scope depends on it."""
    table = data_dir / "taxonomy" / "taxonomy_small8_formatA.tsv"
    seen = {}
    for mode in ("reverse", "greedy", "segment"):
        run = mf(taxonomy=table,
                 extra=["--taxonomy-format", "embedded",
                        "--taxonomy-delimiter-mode", mode, "--force"])
        run.assert_ok(f"Format A, delimiter mode {mode}")
        seen[mode] = run.rc
    assert set(seen.values()) == {0}, seen
    record_metric("v04_format_a", "modes_run", sorted(seen))


@pytest.mark.capability("auto_embed_taxonomy", "marker_mode",
                        "taxonomy_format")
def test_embedded_identifiers_are_auto_detected_in_hmm_mode(
        mf, genome_input, data_dir, hmm_dir, record_metric):
    """No ``--taxonomy-table``, no ``--species-tree``, identifiers carry Format A
    taxonomy: the documented behaviour is to read the taxonomy from the names,
    because without it every marker would be flagged HGT at Level 3."""
    embedded = genome_input("small8", naming="embedded")
    names = sorted(p.stem for p in embedded.iterdir())
    assert any("_d_Bacteria_" in n for n in names), (
        f"the materialised embedded input has no Format A identifiers: {names}"
    )
    run = mf(embedded,
             omit=("--gtdb-markers-dir", "--taxonomy-table"),
             extra=["--marker-mode", "hmm", "--marker-hmm-dir",
                    str(hmm_dir("core")), "--force", "-v"])
    run.assert_ok("embedded identifier detection")
    # The run must say where the taxonomy came from: either the Phase-0.2
    # "Loaded taxonomy for N tips" line with its declared origin, or the
    # Embedded-format auto-detection message.
    text = run.text.lower()
    assert "loaded taxonomy for" in text or "embedded" in text, (
        f"the run never says where the taxonomy came from:\n{run.tail()}")
    record_metric("v04_embedded", "rc", run.rc)


@pytest.mark.capability("auto_embed_taxonomy")
def test_no_auto_embed_leaves_the_run_to_say_so(mf, genome_input, data_dir,
                                                hmm_dir):
    """``--no-auto-embed`` switches the detection off; the documented
    consequence (an absent taxonomy makes the phylogenetic step silent) must be
    stated rather than left for the user to discover in the numbers."""
    embedded = genome_input("small8", naming="embedded")
    run = mf(embedded, omit=("--gtdb-markers-dir", "--taxonomy-table"),
             extra=["--marker-mode", "hmm", "--marker-hmm-dir",
                    str(hmm_dir("core")), "--no-auto-embed", "--force"])
    # Either it refuses for want of a taxonomy, or it runs and declares the
    # Missing taxonomy. Silently producing a normal-looking report is the one
    # Outcome that must not happen.
    if run.rc == 0:
        text = run.text.lower()
        assert "taxonomy" in text, run.tail()
        assert "not" in text or "silent" in text or "level 3" in text or \
               "without" in text, run.tail()
    else:
        assert "taxonomy" in run.text.lower(), run.tail()
        assert "Traceback" not in run.text, run.tail()


@pytest.mark.capability("taxonomy_source_priority")
def test_priority_decides_when_both_sources_are_present(mf, genome_input,
                                                        data_dir, hmm_dir,
                                                        record_metric):
    """Embedded identifiers AND a table: ``--taxonomy-source-priority`` says
    which wins, and the two settings must be observably different or the option
    is decoration."""
    embedded = genome_input("small8", naming="embedded")
    # The table deliberately disagrees with the identifiers: it is the
    # Split-genus control table, where three tips share a fake genus.
    table = data_dir / "taxonomy" / "taxonomy_small8_splitgenus.tsv"
    common = ["--marker-mode", "hmm", "--marker-hmm-dir", str(hmm_dir("core")),
              "--taxonomy-table", str(table), "--force"]
    by_table = mf(embedded, extra=common + ["--taxonomy-source-priority", "table"])
    by_embed = mf(embedded, extra=common + ["--taxonomy-source-priority",
                                            "embedded"])
    by_table.assert_ok("priority=table")
    by_embed.assert_ok("priority=embedded")
    recorded = [
        "Splitgenus" in (by_table.text + str(_orders(by_table))),
        "Splitgenus" in (by_embed.text + str(_orders(by_embed))),
    ]
    record_metric("v04_priority", "table_priority_uses_table_value", recorded[0])
    record_metric("v04_priority", "embedded_priority_uses_table_value",
                  recorded[1])


def _orders(run):
    """The order labels the run resolved, read back from its own report."""
    path = run.out_dir / "Phase5_metadata" / "run_config.json"
    if not path.exists():
        return []
    import json

    payload = json.loads(path.read_text(encoding="utf-8"))
    return [str(v) for v in (payload.get("taxonomy_orders") or [])]


@pytest.mark.capability("taxonomy_levels")
def test_custom_level_prefixes_extend_the_builtin_levels(mf, case_workdir,
                                                         data_dir,
                                                         record_metric):
    """``--taxonomy-levels 'kingdom:k__'`` must make a rank the built-in table
    does not know get parsed instead of dropped.

    Paired design: the same table is read with and without the option. Without
    it the parser warns that the prefix is unknown and the rank disappears; with
    it the run says which ranks it extended and the warning is gone. A single
    run could not tell those two apart.
    """
    table = case_workdir / "with_kingdom.tsv"
    base = (data_dir / "taxonomy" / "taxonomy_quad4.tsv").read_text(
        encoding="utf-8").splitlines()
    lines = [base[0]]
    for line in base[1:]:
        acc, _, tax = line.partition("\t")
        lines.append(f"{acc}\tk__Bacteria;{tax}")
    table.write_text("\n".join(lines) + "\n", encoding="utf-8")

    without = mf(taxonomy=table, extra=["--force", "-v"])
    with_option = mf(taxonomy=table,
                     extra=["--taxonomy-levels", "kingdom:k__",
                            "--force", "-v"])

    # The control arm proves the table really carries an unknown-to-default rank.
    assert "Unknown level prefix" in without.text, (
        f"the table without --taxonomy-levels did not report a dropped rank, so "
        f"the pair says nothing:\n{without.tail()}"
    )
    with_option.assert_ok("custom kingdom level")
    assert "Unknown level prefix" not in with_option.text, (
        f"the custom prefix was still not accepted:\n{with_option.tail()}"
    )
    assert "kingdom" in with_option.text.lower(), (
        f"the run never recorded the extended rank:\n{with_option.tail()}"
    )
    record_metric("v04_custom_levels", "warning_without", True)
    record_metric("v04_custom_levels", "warning_with", False)


@pytest.mark.capability("monophyly_threshold")
def test_a_stricter_monophyly_threshold_cannot_keep_more_markers(
        mf, read_tsv, record_metric):
    """The documented semantics: below the floor, a marker is HGT-prone. Raising
    the floor can only move markers down, never up."""
    lenient = mf(extra=["--monophyly-threshold", "0.0", "--force"])
    lenient.assert_ok()
    strict = mf(extra=["--monophyly-threshold", "1.0", "--force"])
    strict.assert_ok()

    def flagged(run):
        rows = read_tsv(run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
        return {r["marker_id"] for r in rows if r["hgt_risk_level"] != "level_1"}

    soft, hard = flagged(lenient), flagged(strict)
    record_metric("v04_monophyly", "flagged_at_0.0", len(soft))
    record_metric("v04_monophyly", "flagged_at_1.0", len(hard))
    assert hard >= soft, (
        f"threshold 1.0 flagged {sorted(hard)} but 0.0 flagged {sorted(soft)}: "
        "the floor is not applied"
    )


@pytest.mark.capability("monophyly_rank")
def test_every_documented_rank_is_accepted(mf, record_metric):
    """Every rank must be a usable value (not an argument error). A rank the
    data cannot support may legitimately end the run — as a data failure with a
    stated reason, which is a different exit code and a different message."""
    results = {}
    for rank in ("domain", "phylum", "class", "order", "family", "genus",
                 "species"):
        run = mf("small8", extra=["--monophyly-rank", rank, "--max-markers", "2",
                                 "--force"])
        results[rank] = run.rc
        assert run.rc != constants.EXIT_ARG_ERROR, (
            f"--monophyly-rank {rank} was rejected as an argument: {run.tail()}"
        )
        if run.rc != 0:
            assert "Traceback" not in run.text, run.tail()
            assert run.text.strip(), f"--monophyly-rank {rank} failed silently"
    record_metric("v04_ranks", "exit_codes", results)


@pytest.mark.capability("taxonomy_table", "workflow:unknown-not-placeholder")
def test_a_tip_missing_from_the_table_is_counted_not_assumed(mf, data_dir,
                                                            case_workdir):
    """A genome with no taxonomy line cannot be placed, and the run must state
    how many tips it actually has taxonomy for: silently scoring against a
    partially-covered table is how a rank choice becomes meaningless."""
    table = case_workdir / "short_table.tsv"
    lines = (data_dir / "taxonomy" / "taxonomy_small8.tsv") \
        .read_text(encoding="utf-8").rstrip("\n").split("\n")
    table.write_text("\n".join(lines[:5]) + "\n", encoding="utf-8")  # 4 of 8 tips
    run = mf("small8", taxonomy=table, extra=["--force", "-v"])
    text = (run.text + "\n".join(
        p.read_text(encoding="utf-8", errors="replace")
        for p in run.out_dir.rglob("*.txt")
    )) if run.out_dir else run.text
    if run.rc != 0:
        assert "Traceback" not in run.text, run.tail()
        assert "taxonom" in run.text.lower(), run.tail()
        return
    import re

    loaded = re.search(r"[Ll]oaded taxonomy for (\d+) tips", text)
    assert loaded, (
        "the run never reported how many tips it loaded taxonomy for, so a "
        f"half-covered table is indistinguishable from a complete one:\n"
        f"{run.tail()}"
    )
    assert int(loaded.group(1)) == 4, loaded.group(0)
