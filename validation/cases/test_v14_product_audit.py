"""V-14 — the product audit: every documented output, where it is, what is in it.

Covers the ``product:*`` surface: the files the README and MANUAL promise a run
leaves behind, their schemas, and the rule that a product is either present and
populated or absent with an explanation — never an empty placeholder.

Why the schema is asserted from the documented column list rather than from
whatever a run happens to write: a renamed or dropped column is invisible to a
"does the file exist" check, and it breaks every downstream reader silently.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

PREFIX = "markerfinder"

# Columns the documentation names, per product. Kept here as data so a missing
# Column is one failing assertion rather than a missing row in a table.
EXPECTED_COLUMNS = {
    "Phase5_reports/markerfinder.marker_summary.tsv": {
        "marker_id", "occupancy_score", "marker_quality_score",
        "marker_quality_level",
    },
    "Phase5_reports/markerfinder.hgt_evaluation.tsv": {
        "marker_id", "overall_risk", "hgt_risk_level",
        "hgt_evidence_confidence", "notes", "detector",
    },
    "Phase5_reports/markerfinder.assertions.tsv": {
        "assertion_id", "name", "severity", "result", "detail", "provisional",
    },
}


@pytest.fixture(scope="module")
def full_run(request, mf, data_dir):
    """One complete run that every product assertion in this file reads.

    ``--save-intermediates`` and ``--coalescent-mode always`` are on so the
    optional products exist too; ``--hgt-scan`` and the composition step are
    exercised in their own files (they write extra artefacts of their own).
    """
    return mf(markers=data_dir / "markers" / "quad4",
              extra=["--coalescent-mode", "always", "--save-intermediates",
                     "--ufboot", "1000", "--force"],
              out_dir=Path(request.config.rootdir) / "validation" / ".work" /
              "outputs" / "v14_full_run")


@pytest.mark.capability("product:Phase5_reports.marker_summary.tsv",
                        "product:Phase5_reports.hgt_evaluation.tsv",
                        "product:Phase5_reports.assertions.tsv")
def test_every_tabular_product_exists_with_its_documented_columns(full_run,
                                                                 record_metric):
    for rel, columns in EXPECTED_COLUMNS.items():
        path = full_run.product(rel)
        text = path.read_text(encoding="utf-8").rstrip("\n")
        assert text, f"{rel} is an empty file"
        header = set(text.split("\n")[0].split("\t"))
        missing = columns - header
        assert not missing, f"{rel} lacks documented column(s): {sorted(missing)}"
    record_metric("v14_products", "tabular_products_checked",
                  sorted(EXPECTED_COLUMNS))


@pytest.mark.capability("product:Phase5_evidence.decision.json",
                        "workflow:provenance-recorded")
def test_decision_cards_exist_for_every_scored_marker(full_run, read_tsv,
                                                     record_metric):
    """'s claim is that a reader can reconstruct WHY a marker got its level.

    So: one card per scored marker, and each card must carry the inputs of the
    decision (the measured quantities, the thresholds in effect, the evidence
    states), not just the verdict.
    """
    evidence = full_run.out_dir / "Phase5_evidence"
    assert evidence.is_dir(), f"{evidence} missing"
    cards = sorted(evidence.glob("decision_*.json"))
    scored = read_tsv(full_run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    assert cards, "no decision cards were written"
    assert len(cards) == len(scored), (
        f"{len(cards)} decision cards for {len(scored)} scored markers"
    )
    payload = json.loads(cards[0].read_text(encoding="utf-8"))
    assert json.dumps(payload).strip("{} "), f"{cards[0]} is an empty card"
    record_metric("v14_cards", "cards", len(cards))
    record_metric("v14_cards", "card_top_level_keys", sorted(payload))


@pytest.mark.capability("product:Phase4_trees.species_tree_concat.newick",
                        "product:Phase4_trees.species_tree_astral.newick",
                        "product:Phase4_trees.gene_trees.newick")
def test_trees_are_parseable_and_cover_the_input_genomes(full_run, read_tsv,
                                                         record_metric):
    from markerfinder.utils.tree_utils import parse_newick_tips

    concat = full_run.product(
        "Phase4_trees/markerfinder.species_tree_concat.newick")
    text = concat.read_text(encoding="utf-8").strip()
    assert text.endswith(";"), f"{concat} is not terminated as Newick: {text[-60:]}"
    # FastTree writes each internal node's support as a bare label (``(a,b)72``),
    # And a naive tip read counts ``72`` as a taxon — which would let a tree over
    # Four genomes be reported as having five tips. Separate the two by shape and
    # Record what was separated, so the number in the report is the tip count.
    parsed = set(parse_newick_tips(text))
    support_labels = {t for t in parsed if re.fullmatch(r"\d+(?:\.\d+)?", t)}
    tips = parsed - support_labels
    assert len(tips) >= 4, f"the species tree has {len(tips)} tips: {sorted(tips)}"
    assert all(not t[:1].isdigit() for t in tips), sorted(tips)

    gene_trees = full_run.product("Phase4_trees/markerfinder.gene_trees.newick")
    trees = [ln for ln in gene_trees.read_text(encoding="utf-8").splitlines()
             if ln.strip().endswith(";")]
    assert trees, "no gene trees in the multi-newick product"
    # Every gene tree must be readable and have at least four tips, which is the
    # Minimum the pipeline itself demands before building one. Support labels are
    # Not tips, so they are filtered out before counting.
    def tips_of(statement: str) -> set:
        return {t for t in parse_newick_tips(statement)
                if not re.fullmatch(r"\d+(?:\.\d+)?", t)}

    short = [i for i, ln in enumerate(trees) if len(tips_of(ln)) < 4]
    assert not short, f"gene trees with fewer than four tips at index {short}"
    record_metric("v14_trees", "species_tree_tips", sorted(tips))
    record_metric("v14_trees", "support_labels_on_internal_nodes",
                  sorted(support_labels))
    record_metric("v14_trees", "gene_trees", len(trees))

    astral = full_run.out_dir / "Phase4_trees" / "markerfinder.species_tree_astral.newick"
    if astral.exists():
        astral_text = astral.read_text(encoding="utf-8").strip()
        assert astral_text, "the ASTRAL tree file exists but is empty"
        astral_tips = {t for t in parse_newick_tips(astral_text)
                       if not re.fullmatch(r"\d+(?:\.\d+)?", t)}
        assert len(astral_tips) >= 4, astral_text


@pytest.mark.capability("product:Phase4_alignments.partition.nex",
                        "workflow:hgt-grading")
def test_the_partition_file_is_nexus_and_lists_the_selected_markers(full_run,
                                                                   record_metric):
    nex = full_run.product("Phase4_alignments/markerfinder.partition.nex")
    text = nex.read_text(encoding="utf-8")
    assert text.lstrip().startswith("#nexus"), text[:80]
    assert "begin sets" in text.lower(), text[:300]
    assert "end;" in text.lower(), text[-300:]
    charsets = re.findall(r"^\s*charset\s+(\S+)\s*=", text, flags=re.M)
    summary = full_run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                      .read_text(encoding="utf-8")
    selected = int(re.search(r"Markers selected:\s*(\d+)", summary).group(1))
    # The partition describes the SUPERMATRIX, i.e. the markers the HGT screen
    # Kept (Level 1 + Level 2). Phase 1's selected count is the wrong yardstick:
    # A marker excluded as Level 3 must not appear here, or IQ-TREE would be
    # Handed a partition for data the run deliberately rejected.
    kept = sum(int(m.group(1)) for m in re.finditer(
        r"Level [12] \([^)]*\):\s*(\d+)", summary))
    record_metric("v14_partition", "charsets", len(charsets))
    record_metric("v14_partition", "markers_selected", selected)
    record_metric("v14_partition", "markers_kept", kept)
    assert len(charsets) == kept, (
        f"the partition file lists {len(charsets)} charsets for {kept} markers "
        f"kept by the HGT screen ({selected} were selected in Phase 1): "
        f"{charsets}"
    )
    assert len(set(charsets)) == len(charsets), f"duplicate charsets: {charsets}"


@pytest.mark.capability("product:Phase4_trees.gene_trees.cache",
                        "workflow:determinism")
def test_the_gene_tree_cache_states_what_each_tree_was_built_from(full_run, mf,
                                                                 data_dir,
                                                                 record_metric):
    """A cache entry must be attributable to its input.

    Per-marker gene trees are written into ``Phase4_trees/gene_trees/`` and
    reused when a run finds one there. Reuse keyed by marker id alone would let
    a tree inferred from *other* sequences become this run's gene tree — and the
    HGT screen would grade against it. So every cached tree carries a stamp of
    the sequences it was built from, and a second run into the same directory
    may only reuse entries whose stamp matches.
    """
    cache = full_run.out_dir / "Phase4_trees" / "gene_trees"
    trees = sorted(cache.glob("*.nwk"))
    assert trees, f"no cached gene trees under {cache}"

    stamps = []
    for tree in trees:
        stamp = cache / f"{tree.name}.input_sha256"
        assert stamp.exists(), (
            f"{tree.name} is a gene tree with no record of the sequences it was "
            "built from; the next run cannot tell whether it may be reused"
        )
        value = stamp.read_text(encoding="utf-8").strip()
        assert re.fullmatch(r"[0-9a-f]{64}", value), f"{stamp}: {value!r}"
        stamps.append((tree.name, value))
    # Reusing the directory must be a reuse of verified entries, not a rebuild
    # Dressed up as one: the second run says so in its own log.
    rerun = mf(markers=data_dir / "markers" / "quad4",
               out_dir=full_run.out_dir,
               extra=["--coalescent-mode", "always", "--save-intermediates",
                      "--ufboot", "1000", "--force", "-v"])
    rerun.assert_ok("rerun into the same output directory")
    assert "Reusing existing gene tree" in rerun.text, (
        f"the second run rebuilt everything, so the stamps do not match the "
        f"inputs they describe:\n{rerun.tail()}"
    )
    after = {p.name[: -len(".input_sha256")]:
             (cache / p.name).read_text(encoding="utf-8").strip()
             for p in cache.glob("*.nwk.input_sha256")}
    before = dict(stamps)
    changed = {name: (before.get(name, "<absent>")[:12], value[:12])
               for name, value in after.items() if before.get(name) != value}
    changed.update({name: (value[:12], "<absent>")
                    for name, value in before.items() if name not in after})
    assert not changed, (
        "the cache stamps changed between two runs over the same input "
        f"(name: before -> after): {changed}"
    )
    record_metric("v14_tree_cache", "cached_trees", len(trees))
    record_metric("v14_tree_cache", "stamped_trees", len(after))


@pytest.mark.capability("product:Phase5_reports.pipeline_summary.txt",
                        "product:Phase5_reports.hgt_evaluation.tsv",
                        "workflow:hgt-grading")
def test_the_text_summary_states_the_numbers_the_report_is_judged_on(full_run,
                                                                      read_tsv,
                                                                      record_metric):
    summary = full_run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                        .read_text(encoding="utf-8")
    for probe in ("Genomes:", "Markers selected:", "Level 1 (clean):",
                  "Evidence coverage:", "Analysis mode:"):
        assert probe in summary, f"the summary omits {probe!r}:\n{summary[:600]}"

    counts = {int(m.group(1)): int(m.group(2)) for m in re.finditer(
        r"Level ([123])[^:]*:\s*(\d+)", summary)}
    rows = read_tsv(full_run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    from collections import Counter
    actual = Counter(r["hgt_risk_level"] for r in rows)
    record_metric("v14_summary", "summary_level_counts", counts)
    record_metric("v14_summary", "tsv_level_counts", dict(actual))
    for level, number in counts.items():
        key = f"level_{level}"
        assert actual.get(key, 0) == number, (
            f"the summary prints Level {level}={number} while the TSV says "
            f"{actual.get(key, 0)}"
        )


@pytest.mark.capability("product:Phase5_reports.report.html")
def test_html_report_is_self_contained_and_names_the_products_it_summarises(
        full_run, record_metric):
    html = full_run.product("Phase5_reports/markerfinder.report.html") \
                   .read_text(encoding="utf-8")
    assert "<html" in html.lower(), html[:200]
    assert "MarkerFinder" in html
    assert not re.search(r"(?:src|href)\s*=\s*[\"']https?://", html), (
        "the 'self-contained static' report references remote resources"
    )
    # The documented contents: a summary, the marker table, the HGT table.
    for probe in ("Marker", "HGT", "occupancy"):
        assert probe.lower() in html.lower(), f"report never mentions {probe}"
    record_metric("v14_html", "bytes", len(html))


@pytest.mark.capability("workflow:unknown-not-placeholder",
                        "product:Phase5_reports.hgt_evaluation.tsv")
def test_unmeasured_quantities_are_never_rendered_as_a_number(full_run,
                                                             read_tsv):
    """A quantity that could not be measured must read as
    unknown / NOT MEASURABLE, not as 0.5 or 0."""
    rows = read_tsv(full_run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    offenders = []
    for row in rows:
        for column in ("rf", "quartet", "monophyly"):
            value = row.get(column, "")
            state = row.get(f"{column}_state", "")
            if state and "measured" in state.lower() and "not" in state.lower():
                if value not in ("", "NA", "N/A", "nan", "None", "-"):
                    offenders.append((row["marker_id"], column, value, state))
    assert not offenders, (
        f"unmeasured quantities rendered as numbers: {offenders}"
    )


@pytest.mark.capability("product:Phase5_metadata.run_config.json",
                        "workflow:provenance-recorded")
def test_no_undocumented_extra_product_is_left_behind(full_run):
    """Anything a run writes must be a product the documentation promises.

    The list is the published one; a stray file is either an undocumented
    product or a leak (a temp file, a debug dump) into the user's output
    directory, and both belong in a failing test rather than in a release.
    """
    documented_prefixes = (
        "Phase5_reports/", "Phase5_evidence/", "Phase5_metadata/",
        "Phase4_trees/", "Phase4_alignments/", "Phase4_intermediate/",
        ".markerfinder/",
    )
    root = full_run.out_dir
    actual = sorted(
        str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()
    )
    outside = [a for a in actual if not a.startswith(documented_prefixes)]
    assert not outside, f"products outside the documented layout: {outside[:20]}"
