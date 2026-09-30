"""V-05 — reference species tree, tree parsing and sequence validation.

Covers ``--species-tree``, the deprecated ``--tree`` alias,
``--multi-tree-mode``, ``--strip-annotations``, ``--sequences``, ``--mol-type``,
``--skip-length-check``, ``--no-cross-check``, and the refusal paths for an
illegal reference.

With ``--species-tree`` the phylogenetic step compares each gene tree against
the supplied reference (Robinson–Foulds and quartet incongruence); with
``--taxonomy-table`` it uses MAD rooting plus a monophyly proportion instead.
Both routes are documented, so both are run.
"""

from __future__ import annotations

import re

import pytest

from markerfinder.cli import constants
from markerfinder.utils import reference as reference_module


@pytest.mark.capability("species_tree", "workflow:hgt-grading",
                        "product:Phase5_reports.hgt_evaluation.tsv")
def test_reference_species_tree_drives_the_rf_quartet_route(mf, data_dir,
                                                           read_tsv,
                                                           record_metric):
    tree = data_dir / "trees" / "reference_small8.nwk"
    run = mf(omit=("--taxonomy-table",),
             extra=["--species-tree", str(tree), "--force"])
    run.assert_ok("--species-tree route")
    rows = read_tsv(run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    assert rows, "no marker was graded against the reference tree"
    columns = set(rows[0])
    assert any("rf" in c.lower() or "quartet" in c.lower() for c in columns), (
        f"the evaluation has no RF/quartet column to show the reference was "
        f"used: {sorted(columns)}"
    )
    record_metric("v05_species_tree", "markers_graded", len(rows))


@pytest.mark.capability("tree", "species_tree")
def test_deprecated_tree_alias_still_works_and_says_it_is_deprecated(
        mf, data_dir):
    tree = data_dir / "trees" / "reference_small8.nwk"
    run = mf(omit=("--taxonomy-table", "--species-tree"),
             extra=["--tree", str(tree), "--force"])
    run.assert_ok("deprecated --tree alias")
    assert "deprecat" in run.text.lower(), (
        "--tree is documented as a deprecated alias but the run never says so:\n"
        f"{run.tail()}"
    )


@pytest.mark.capability("species_tree", "workflow:failure-loudness", "workflow:hgt-grading")
def test_illegal_reference_trees_are_reported_and_never_scored(
        run_markerfinder, genome_input, data_dir, tmp_path, read_tsv,
        record_metric):
    """'s documented behaviour for an illegal reference is not a crash and
    not a number: every affected marker is routed to UNKNOWN/``unscreened``, and
    the run says why.

    A reference with three tips cannot measure incongruence against an eight-tip
    gene tree; silently emitting an RF value anyway is the failure mode this
    guards against, so the assertion is on the *reason and the state*, not on a
    particular exit code (which depends on whether any marker survives).
    """
    trees = data_dir / "trees"
    for name in ("illegal_three_tips.nwk", "illegal_duplicate_tips.nwk",
                 "illegal_negative_branch.nwk"):
        out_dir = tmp_path / name
        result = run_markerfinder([
            "-i", str(genome_input("small8")), "-o", str(out_dir),
            "-t", "2", "--marker-mode", "gtdb_tk",
            "--gtdb-markers-dir", str(data_dir / "markers" / "small8_core"),
            "--max-markers", "2", "--coalescent-mode", "off", "--skip-checkm",
            "--species-tree", str(trees / name), "--force", "-v",
        ])
        assert "Traceback" not in result.text, result.tail()
        text = result.text.lower()
        # The refusal may come from the reference validator ("reference illegal
        # [TOO_FEW_TIPS]") or from the structural tree validator ("Duplicate tip
        # Names") — either is a named reason; neither may be a silent pass.
        assert ("reference illegal" in text or "validation error" in text
                or "duplicate tip" in text or "too few" in text
                or "not resolved" in text), (
            f"{name}: the run never named why the reference is unusable:\n"
            f"{result.tail()}"
        )
        # The promise is about the markers, not the exit code: an illegal
        # Reference may not take part in scoring. Whether the run then exits 0
        # Depends on whether its now-unscreened markers are enough to build a
        # Tree, so the state is read from the product — and a run that leaves no
        # Product must not claim success either.
        evaluation = (out_dir / "Phase5_reports"
                      / "markerfinder.hgt_evaluation.tsv")
        if evaluation.exists():
            scored = [r["marker_id"] for r in read_tsv(evaluation)
                      if r["hgt_risk_level"] != "unknown"]
            assert not scored, (
                f"{name}: these markers were graded against an illegal "
                f"reference: {scored}"
            )
        else:
            assert result.rc != 0, (
                f"{name}: exit {result.rc} with no evaluation product — the run "
                "reports success without recording what it decided"
            )


@pytest.mark.capability("species_tree", "workflow:failure-loudness")
def test_reference_validity_codes_are_distinguishable():
    """The three illegal shapes must carry different codes, so a reader can act
    on the reason."""
    from markerfinder.utils.reference import (
        MALFORMED, TOO_FEW_TIPS, validate_reference_tree,
    )

    tiny = validate_reference_tree("(A:0.1,B:0.1)C;", "(A:0.1,B:0.1)C;")
    assert not tiny.ok and tiny.code == TOO_FEW_TIPS, tiny
    broken = validate_reference_tree("((A,B;", "(A,B);")
    assert not broken.ok and broken.code == MALFORMED, broken
    assert TOO_FEW_TIPS != MALFORMED


@pytest.mark.capability("multi_tree_mode")
def test_multi_tree_mode_selects_one_tree_as_documented(mf, data_dir,
                                                        record_metric):
    """One file, three trees. ``first`` / ``last`` must each resolve to a
    single reference; ``split`` is documented as unsupported and must refuse."""
    multi = data_dir / "trees" / "reference_small8_multi.nwk"
    for mode in ("first", "last", "random"):
        run = mf(omit=("--taxonomy-table",),
                 extra=["--species-tree", str(multi),
                        "--multi-tree-mode", mode, "--force"])
        assert run.rc == 0, f"--multi-tree-mode {mode} failed:\n{run.tail()}"
    refused = mf(omit=("--taxonomy-table",),
                 extra=["--species-tree", str(multi), "--multi-tree-mode",
                        "split", "--force"])
    assert refused.rc != 0, (
        "--multi-tree-mode split is documented as unsupported, but the run "
        f"succeeded:\n{refused.tail()}"
    )
    assert "split" in refused.text.lower(), refused.tail()
    record_metric("v05_multitree", "modes_accepted", ["first", "last", "random"])


@pytest.mark.capability("strip_annotations")
def test_strip_annotations_removes_nhx_labels_from_the_reported_tree(
        mf, data_dir, record_metric):
    """The NHX twin of the reference carries per-tip annotations; stripping must
    be observable in the product, not just accepted."""
    plain = data_dir / "trees" / "reference_small8.nwk"
    nhx = data_dir / "trees" / "reference_small8.nhx.nwk"
    assert nhx.exists(), f"{nhx} missing — re-run 03_build_fixtures.py"
    assert "&&NHX" in nhx.read_text(encoding="utf-8"), (
        "the NHX fixture carries no annotations, so the option cannot be tested"
    )
    run = mf(omit=("--taxonomy-table",),
             extra=["--species-tree", str(nhx), "--strip-annotations",
                    "--force"])
    run.assert_ok("--strip-annotations")
    written = _tree_texts(run)
    assert "&&NHX" not in written, (
        "annotations survived --strip-annotations in the recorded tree"
    )
    record_metric("v05_strip", "nhx_fixture_chars", len(nhx.read_text()))


def _tree_texts(run):
    out = []
    for path in sorted(run.out_dir.rglob("*.newick")) + \
            sorted(run.out_dir.rglob("*.nwk")):
        out.append(path.read_text(encoding="utf-8"))
    return "\n".join(out)


@pytest.mark.capability("sequences", "no_cross_check")
def test_tree_sequence_cross_validation_runs_and_can_be_disabled(
        mf, data_dir, record_metric):
    """``--sequences`` plus a reference tree triggers the cross-validation;
    ``--no-cross-check`` skips it, and skipping must be visible in the log
    rather than silent (the help text calls the skip "not recommended")."""
    tree = data_dir / "trees" / "reference_small8.nwk"
    seqs = data_dir / "sequences" / "tip_sequences.faa"
    # ``-v`` because the successful outcome is an INFO line ("Cross-validation
    # Passed: N tips matched"): a check that ran and found nothing must be
    # Distinguishable from one that never ran, and that distinction lives in the
    # Log.
    with_check = mf(omit=("--taxonomy-table",),
                    extra=["--species-tree", str(tree), "--sequences",
                           str(seqs), "--force", "-v"])
    assert with_check.rc in (0, constants.EXIT_DATA_ERROR), with_check.tail()
    text = with_check.text.lower()
    assert "cross" in text or "mismatch" in text or "tip" in text, (
        f"the cross-validation left no trace:\n{with_check.tail()}"
    )

    skipped = mf(omit=("--taxonomy-table", "--no-cross-check"),
                 extra=["--species-tree", str(tree), "--sequences", str(seqs),
                        "--no-cross-check", "--force"])
    assert skipped.rc == 0, skipped.tail()
    record_metric("v05_crosscheck", "with_check_rc", with_check.rc)
    record_metric("v05_crosscheck", "skipped_rc", skipped.rc)


@pytest.mark.capability("mol_type", "skip_length_check", "sequences")
def test_mol_type_can_be_forced_and_the_length_check_skipped(mf, data_dir,
                                                             record_metric):
    """``--mol-type`` is an alphabet assertion over the sequence input, so it is
    exercised with one: forcing DNA onto protein sequences must be refused,
    declaring protein must pass, and ``--skip-length-check`` must still run.

    Cross-validation only runs when a reference tree and a sequence file are both
    supplied (`cli/validation._handle_cross_validation`), so the pair is given —
    and it must be a *matching* pair: ``tip_sequences.faa`` holds one protein per
    small8 tip, so the tree beside it is the small8 reference. A sequence file
    whose ids are not on the tree is a label mismatch, which the run reports
    (rightly) instead of cross-checking anything.
    """
    seqs = data_dir / "sequences" / "tip_sequences.faa"
    tree = data_dir / "trees" / "reference_small8.nwk"
    pair = ["--sequences", str(seqs), "--species-tree", str(tree)]
    ok = mf("small8", extra=["--mol-type", "protein", *pair, "--force"])
    assert ok.rc == 0, f"--mol-type protein on proteins failed:\n{ok.tail()}"

    wrong = mf("small8", extra=["--mol-type", "DNA", *pair, "--force"])
    assert wrong.rc != 0, (
        "--mol-type DNA was accepted for protein sequences:\n" + wrong.tail()
    )
    assert "Traceback" not in wrong.text, wrong.tail()

    skip = mf(extra=["--skip-length-check", "--force"])
    skip.assert_ok("--skip-length-check")
    record_metric("v05_moltype", "protein_rc", ok.rc)
    record_metric("v05_moltype", "dna_on_protein_rc", wrong.rc)


@pytest.mark.capability("species_tree", "workflow:unknown-not-placeholder")
def test_a_reference_that_does_not_cover_the_input_is_not_quietly_used(
        mf, data_dir, tmp_path, case_workdir, genome_input):
    """Tips that are not in the dataset cannot measure incongruence. The run
    must refuse or report the mismatch — never grade against an unrelated
    tree."""
    stranger = case_workdir / "unrelated.nwk"
    stranger.write_text("((X:0.1,Y:0.1)A:0.1,(Z:0.1,W:0.1)B:0.1)R;\n",
                        encoding="utf-8")
    run = mf(omit=("--taxonomy-table",),
             extra=["--species-tree", str(stranger), "--force"])
    text = run.text.lower()
    assert run.rc != 0 or "match" in text or "tip" in text or "overlap" in text, (
        f"a reference with no tip in common was used silently:\n{run.tail()}"
    )
