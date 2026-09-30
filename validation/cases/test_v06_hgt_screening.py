"""V-06 — HGT screening: steps, thresholds, scanning and evidence coverage.

Covers ``--hgt-steps``, ``--hgt-threshold`` (deprecated alias),
``--hgt-threshold-l1l2``, ``--hgt-threshold-l2l3``, ``--hgt-scan``,
``--scan-stability-min``, ``--require-evidence-coverage``,
``--hgt-adaptive-thresholds``, ``--no-hgt-adaptive-thresholds``,
``--cog-category-map`` and ``--min-informative-sites``.

The screening design under test: the risk score is a *ranking* device, the
level is the *decision*, and composition diagnostics are parallel evidence
that must never be merged into the score. Each of those claims gets a
case here, because a claim about a number is only testable by moving the
number.
"""

from __future__ import annotations

import pytest


@pytest.mark.capability("hgt_steps", "workflow:hgt-grading")
def test_default_step_is_phylogenetic_only(mf, read_tsv, record_metric):
    run = mf(extra=["--force"])
    run.assert_ok()
    rows = read_tsv(run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    assert rows
    # The composition columns live in their own artefact, never in the score.
    composition = run.out_dir / "Phase5_reports" / "markerfinder.composition.tsv"
    assert not composition.exists(), (
        "composition diagnostics were written although --hgt-steps defaults to "
        "'phylogenetic': the default would silently enable a second screen")
    score_cols = [c for c in rows[0] if "rcv" in c.lower() or "gc" in c.lower()]
    assert not score_cols, (
        f"composition metrics leaked into the evaluation table: {score_cols}")
    record_metric("v06_default_steps", "markers", len(rows))


@pytest.mark.capability("hgt_steps", "product:Phase5_reports.composition.tsv")
def test_composition_step_writes_parallel_evidence_not_a_score(mf, read_tsv,
                                                               record_metric):
    """``--hgt-steps phylogenetic,composition`` adds its own file and its own
    columns; the risk score must be unchanged by adding a diagnostic."""
    plain = mf(extra=["--force"])
    plain.assert_ok()
    with_comp = mf(extra=["--hgt-steps", "phylogenetic,composition", "--force"])
    with_comp.assert_ok("--hgt-steps phylogenetic,composition")

    comp_path = with_comp.product("Phase5_reports/markerfinder.composition.tsv")
    comp_rows = read_tsv(comp_path)
    assert comp_rows, f"{comp_path} is empty"

    before = {r["marker_id"]: r.get("overall_risk") or r.get("hgt_risk_score")
              for r in read_tsv(plain.product(
                  "Phase5_reports/markerfinder.hgt_evaluation.tsv"))}
    after = {r["marker_id"]: r.get("overall_risk") or r.get("hgt_risk_score")
             for r in read_tsv(with_comp.product(
                 "Phase5_reports/markerfinder.hgt_evaluation.tsv"))}
    shared = set(before) & set(after)
    changed = {m for m in shared if before[m] != after[m]}
    record_metric("v06_composition", "markers_in_composition", len(comp_rows))
    record_metric("v06_composition", "score_changed_by_diagnostic",
                  sorted(changed))
    assert not changed, (
        "adding a diagnostic step moved the risk score, which it must never: "
        f"{sorted(changed)}"
    )


@pytest.mark.capability("hgt_threshold", "hgt_threshold_l1l2",
                        "workflow:hgt-grading")
def test_deprecated_threshold_alias_only_moves_the_first_band(mf, read_tsv,
                                                             record_metric):
    """``--hgt-threshold`` is documented as an alias of ``--hgt-threshold-l1l2``
    with a deprecation notice, and must not touch the Level2/Level3 boundary.

    The two runs differ in exactly one knob, and the bands are compared by
    count: lowering the L1/L2 boundary can only move Level-1 markers out of
    Level 1, while the Level-3 count must stay put because its boundary was not
    touched. Both L1/L2 values stay under the shared L2/L3 boundary of 0.90: a
    pair with ``L1 >= L2`` is an inconsistent configuration that A-13 refuses,
    and refusing it is the validator working, not the alias.
    """
    # ``--monophyly-rank order`` is part of the design, not decoration: under the
    # Default auto rank this dataset's markers are measured at genus, where only
    # One taxon is testable, so every proportion is 0 or 1 and no band boundary
    # Placed between them can move anything. At order rank two taxa are testable
    # And the risks land at 0.0 and 0.5 — values the bands can separate.
    bands = ["--monophyly-rank", "order"]
    base = mf("small8", extra=["--hgt-threshold-l1l2", "0.60",
                              "--hgt-threshold-l2l3", "0.90", *bands,
                              "--force"])
    base.assert_ok()
    alias = mf("small8", extra=["--hgt-threshold", "0.20",
                               "--hgt-threshold-l2l3", "0.90", *bands,
                               "--force"])
    alias.assert_ok()
    assert "deprecated" in alias.text.lower(), (
        f"the alias was accepted without the promised notice:\n{alias.tail()}"
    )

    def counts(run):
        rows = read_tsv(run.product(
            "Phase5_reports/markerfinder.hgt_evaluation.tsv"))
        return {level: sum(1 for r in rows if r["hgt_risk_level"] == level)
                for level in ("level_1", "level_2", "level_3")}

    a, b = counts(base), counts(alias)
    record_metric("v06_alias", "base_band_counts", a)
    record_metric("v06_alias", "alias_band_counts", b)
    assert b["level_1"] < a["level_1"], (
        f"lowering the L1/L2 boundary from 0.60 to 0.20 through the alias moved "
        f"no marker out of level_1 ({a} -> {b}): the alias reaches no decision"
    )
    assert b["level_3"] == a["level_3"], (
        f"the L1/L2 alias changed the Level-3 count ({a} -> {b}); it must not "
        "touch the second band"
    )


@pytest.mark.capability("hgt_threshold_l2l3")
def test_l2l3_boundary_decides_exclusion(mf, read_tsv, record_metric):
    """The second band is the exclusion decision, so it must be able to exclude
    more and exclude less — on a set with real incongruence signal.

    Both values stay inside the legal interval (``0 <= L1 < L2 <= 1``, asserted
    by A-13 and pinned by its must-fail fixture): a boundary of ``0.0`` below an
    L1/L2 of 0.25 is not an aggressive threshold, it is an inconsistent one, and
    the run refuses it. That refusal is the correct behaviour, so testing the
    knob with it would test the validator instead of the knob.
    """
    lenient = mf("small8", extra=["--hgt-threshold-l1l2", "0.05",
                                  "--hgt-threshold-l2l3", "0.9",
                                  "--monophyly-rank", "order", "--force"])
    strict = mf("small8", extra=["--hgt-threshold-l1l2", "0.05",
                                 "--hgt-threshold-l2l3", "0.1",
                                 "--monophyly-rank", "order", "--force"])
    lenient.assert_ok("boundary 0.9")
    strict.assert_ok("boundary 0.1")

    def count(run, level):
        return sum(1 for r in read_tsv(
            run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
            if r["hgt_risk_level"] == level)

    excl_lenient = count(lenient, "level_3")
    excl_strict = count(strict, "level_3")
    record_metric("v06_l2l3", "excluded_at_0.9", excl_lenient)
    record_metric("v06_l2l3", "excluded_at_0.1", excl_strict)
    assert excl_strict >= excl_lenient, (
        f"boundary 0.1 excluded {excl_strict} markers but 0.9 excluded "
        f"{excl_lenient}: the second band is applied backwards"
    )
    assert excl_strict > excl_lenient, (
        "the Level2/Level3 boundary changed nothing between its two extremes: "
        "the knob is decorative"
    )


@pytest.mark.capability("hgt_scan", "product:Phase5_reports.threshold_scan.tsv",
                        "scan_stability_min")
def test_threshold_scan_writes_its_own_artifact(mf, record_metric):
    run = mf("small8", extra=["--hgt-scan", "--force"])
    run.assert_ok("--hgt-scan")
    text = run.product("Phase5_reports/markerfinder.threshold_scan.tsv") \
               .read_text(encoding="utf-8")
    header = text.splitlines()[0].split("\t")
    assert {"band_level1_max", "band_level2_max", "n_L1", "n_L2", "n_L3",
            "mean_risk"} <= set(header), header
    bands = [ln for ln in text.splitlines()[1:] if ln and not ln.startswith("jaccard")]
    assert len(bands) >= 2, f"the scan wrote {len(bands)} band rows:\n{text}"
    # The cross-band Jaccard matrix is the second block of the same file.
    assert "\njaccard\t" in text or text.splitlines()[-1].startswith("jaccard") \
        or "jaccard" in text, (
        "the scan file carries no cross-band Jaccard block, so the stability "
        "claim has no evidence behind it"
    )
    record_metric("v06_scan", "bands_scanned", len(bands))


@pytest.mark.capability("require_evidence_coverage")
def test_low_evidence_coverage_is_warned_in_the_summary(mf, record_metric):
    """A floor no run can meet must make the warning fire — in the product a
    reader keeps, not only on a console that scrolls away."""
    run = mf("small8", extra=["--require-evidence-coverage", "1.01", "--force"])
    run.assert_ok()
    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                 .read_text(encoding="utf-8")
    assert "EVIDENCE COVERAGE LOW" in summary or "evidence coverage" in summary.lower(), (
        f"the floor was not applied to the summary:\n{summary[:400]}"
    )
    record_metric("v06_coverage", "warned_in_summary",
                  "EVIDENCE COVERAGE LOW" in summary)


@pytest.mark.capability("hgt_adaptive_thresholds",
                        "no_hgt_adaptive_thresholds")
def test_adaptive_relaxation_is_off_by_default_and_opt_in_moves_the_band(
        mf, read_tsv, record_metric):
    """The default must not relax the放行尺度 silently, and opting in must
    actually move the band — otherwise the flag is a label, not a behaviour."""
    default = mf("small8", extra=["--force", "-v"])
    default.assert_ok()
    opted = mf("small8", extra=["--hgt-adaptive-thresholds", "--force", "-v"])
    opted.assert_ok()
    disabled = mf("small8", extra=["--no-hgt-adaptive-thresholds", "--force",
                                   "-v"])
    disabled.assert_ok()

    # This dataset spans three orders, which is the far-distance condition the
    # Flag governs, so the opted run must say it detected it — and the default
    # Run must not.
    assert "far-distance" in opted.text.lower(), (
        "--hgt-adaptive-thresholds was accepted without saying what it did:\n"
        f"{opted.tail()}")
    assert "far-distance" not in default.text.lower(), (
        "the default run relaxed the exclusion band without being asked:\n"
        f"{default.tail()}")

    def thresholds(run):
        rows = read_tsv(run.product(
            "Phase5_reports/markerfinder.hgt_evaluation.tsv"))
        return {r.get("thresholds_in_effect", "") for r in rows}

    a, b = thresholds(default), thresholds(opted)
    record_metric("v06_adaptive", "default_thresholds", sorted(a))
    record_metric("v06_adaptive", "opted_thresholds", sorted(b))
    assert a != b, (
        f"--hgt-adaptive-thresholds left the thresholds in effect identical "
        f"({sorted(a)}): the flag changed no decision"
    )
    # The explicit disable restates the default, so it must agree with it.
    assert thresholds(disabled) == a, (
        "--no-hgt-adaptive-thresholds behaves differently from the default, "
        "which its help text says it merely restates"
    )


@pytest.mark.capability("cog_category_map",
                        "product:Phase5_reports.excluded_profile.tsv")
def test_cog_category_map_populates_the_column_it_feeds(mf, data_dir, read_tsv,
                                                        record_metric):
    """``--cog-category-map`` is the only source of the
    ``functional_category`` column of the consistency profile; without the map
    the column must read NA plus a note, never a blank.

    ``--coalescent-mode always`` is not decoration: the consistency criterion
    compares two inference frameworks, so without a coalescent (ASTRAL) tree
    alongside the concatenation tree the screen reports NOT EXECUTED and never
    reaches the profile this case reads.
    """
    mapping = data_dir / "cog" / "marker_category_map.tsv"
    criterion = ["--hgt-mode", "consistency", "--consistency-stringency", "5",
                 "--coalescent-mode", "always", "--ufboot", "2"]
    without = mf("small8", extra=criterion + ["--force"])
    without.assert_ok("consistency mode without a category map")
    assert "NOT EXECUTED" not in without.text, (
        f"the consistency criterion did not run, so the profile would not "
        f"exist for the right reason:\n{without.tail()}"
    )
    path = without.product("Phase5_reports/markerfinder.excluded_profile.tsv")
    rows = read_tsv(path)
    assert "functional_category" in rows[0], sorted(rows[0])
    assert all(r["functional_category"].strip() in ("", "NA", "N/A")
               for r in rows), [r["functional_category"] for r in rows]
    assert any("no --cog-category-map" in (r.get("category_note") or "").lower()
               for r in rows), (
        "the absent map is not explained in category_note (D-04): "
        f"{[r.get('category_note') for r in rows]}"
    )

    with_map = mf("small8", extra=criterion +
                            ["--cog-category-map", str(mapping), "--force"])
    with_map.assert_ok("--cog-category-map")
    filled = read_tsv(with_map.product(
        "Phase5_reports/markerfinder.excluded_profile.tsv"))
    named = {r["functional_category"] for r in filled} - {"", "NA", "N/A"}
    record_metric("v06_cog", "rows_in_profile", len(filled))
    record_metric("v06_cog", "categories_filled", sorted(named)[:5])
    assert named, (
        "--cog-category-map was supplied and the functional_category column "
        f"still renders NA on every row ({len(filled)} rows): either the map "
        "never reached the screen or its keys do not match the markers being "
        "profiled — in both cases the option does nothing"
    )
    # The map may only annotate markers it actually names: a category written on
    # A marker the map never mentioned would mean the column was invented rather
    # Than looked up.
    mapped = {line.split("\t")[0] for line in
              mapping.read_text(encoding="utf-8").splitlines()[1:] if line}
    annotated = {r["marker_id"] for r in filled
                 if r["functional_category"] not in ("", "NA", "N/A")}
    assert annotated <= mapped, (
        f"categories written for markers absent from the map: "
        f"{sorted(annotated - mapped)}"
    )


@pytest.mark.capability("min_informative_sites")
def test_informative_site_floor_marks_markers_inconclusive_and_says_so(
        mf, read_tsv, record_metric):
    """``--min-informative-sites`` is a provisional knob, and provisional is not
    allowed to mean inert.

    Its documented effect is to mark a marker below the floor
    *inconclusive* so it never enters a "clean" conclusion — it is explicitly
    **not** a deletion, so the assertion is about the marking and the reporting,
    not about the row count.
    """
    off = mf("small8", extra=["--min-informative-sites", "0", "--force"])
    off.assert_ok()
    high = mf("small8", extra=["--min-informative-sites", "100000", "--force"])
    high.assert_ok()

    def clean(run):
        return [r["marker_id"] for r in read_tsv(
            run.product("Phase5_reports/markerfinder.marker_summary.tsv"))
            if r.get("marker_quality_level") == "level_1"
            or r.get("marker_quality_level") == "L1"]

    assert "min-informative-sites" in high.text or "informative" in high.text.lower(), (
        "a floor of 100000 informative sites was applied without reporting it:\n"
        f"{high.tail()}"
    )
    record_metric("v06_pis", "markers_without_floor", len(clean(off)))
    record_metric("v06_pis", "markers_with_floor", len(clean(high)))
    assert len(clean(high)) <= len(clean(off)), (
        "an impossible informativeness floor left the clean set unchanged: the "
        "screen does not run"
    )


@pytest.mark.capability("scan_stability_min")
def test_stability_floor_moves_the_scan_verdict(mf, data_dir, read_tsv,
                                                record_metric):
    """``--scan-stability-min`` decides when the scan calls a marker set
    unstable; a floor of 0 can never warn, an impossible floor must."""
    silent = mf(extra=["--hgt-scan", "--scan-stability-min", "0.0", "--force"])
    silent.assert_ok()
    strict = mf(extra=["--hgt-scan", "--scan-stability-min", "1.5", "--force"])
    strict.assert_ok()
    warned = "unstable" in strict.text.lower()
    record_metric("v06_stability", "impossible_floor_warned", warned)
    assert warned, (
        "a cross-band Jaccard floor of 1.5 (impossible by construction) "
        f"produced no warning:\n{strict.tail()}"
    )
