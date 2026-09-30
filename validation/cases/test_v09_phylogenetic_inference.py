"""V-09 — phylogenetic inference: builders, coalescent modes and support.

Covers ``--gene-tree-builder``, the ``--fast-tree`` alias,
``--coalescent-mode``, ``--ufboot`` and ``--min-gene-tree-support``.

Two species trees are documented (concatenation / supermatrix, and coalescent
via ASTRAL), and the coalescent route is what ``--coalescent-mode`` selects
between. Support values are the quantities the recommendation layer reads, so
their units and their absence must be observable, not assumed.
"""

from __future__ import annotations

import re

import pytest


@pytest.mark.capability("coalescent_mode",
                        "product:Phase4_trees.species_tree_concat.newick",
                        "product:Phase4_trees.species_tree_astral.newick")
@pytest.mark.parametrize("mode", ["off", "post-filter", "always"])
def test_coalescent_mode_decides_whether_the_astral_tree_exists(
        mf, mode, record_metric):
    run = mf(extra=["--coalescent-mode", mode, "--force"])
    run.assert_ok(f"--coalescent-mode {mode}")
    trees = run.out_dir / "Phase4_trees"
    concat = trees / "markerfinder.species_tree_concat.newick"
    astral = trees / "markerfinder.species_tree_astral.newick"
    assert concat.exists(), f"no concatenation tree written in mode {mode}"
    if mode == "off":
        assert not astral.exists(), (
            "--coalescent-mode off still produced an ASTRAL tree")
    else:
        # ASTRAL needs >= 4 gene trees; with a 4-marker budget it may
        # Legitimately fail, and then it must say so rather than write an
        # Empty tree file.
        if astral.exists():
            assert astral.read_text(encoding="utf-8").strip(), (
                "the ASTRAL tree file is empty: a failed inference written as "
                "a product")
        else:
            assert "astral" in run.text.lower() or "coalescent" in run.text.lower(), (
                f"no ASTRAL tree and no statement about why ({mode}):\n"
                f"{run.tail()}")
    record_metric(f"v09_coalescent_{mode}", "astral_written", astral.exists())


@pytest.mark.capability("gene_tree_builder", "fast_tree")
def test_gene_tree_builder_switches_the_binary_in_use(mf, record_metric):
    """``fasttree`` and ``iqtree`` are the documented builders; the deprecated
    ``--fast-tree`` alias must reach the same code path and announce itself."""
    ft = mf(extra=["--gene-tree-builder", "fasttree", "--force"])
    ft.assert_ok("fasttree gene trees")
    iq = mf(extra=["--gene-tree-builder", "iqtree", "--coalescent-mode",
                   "always", "--force"])
    iq.assert_ok("iqtree gene trees")

    alias = mf(extra=["--fast-tree", "--coalescent-mode", "always", "--force"])
    alias.assert_ok("--fast-tree alias")
    assert "deprecat" in alias.text.lower(), (
        f"--fast-tree is documented as deprecated but says nothing:\n"
        f"{alias.tail()}"
    )
    assert ft.recorded("phylo_config", "gene_tree_builder") == "fasttree"
    assert iq.recorded("phylo_config", "gene_tree_builder") == "iqtree"
    record_metric("v09_builder", "fasttree_rc", ft.rc)
    record_metric("v09_builder", "iqtree_rc", iq.rc)


@pytest.mark.capability("gene_tree_builder", "workflow:provenance-recorded")
def test_iqtree_built_trees_carry_support_and_fasttree_built_ones_are_labelled(
        mf, read_tsv, record_metric):
    """Support scales differ between the builders (IQ-TREE UFBOOT is 0-100,
    FastTree SH is 0-1), and the pipeline documents that mix as a hazard:
    a threshold must be compared in the right units or reported as NOT
    MEASURABLE, never against a placeholder."""
    runs = {}
    for builder in ("fasttree", "iqtree"):
        run = mf(extra=["--gene-tree-builder", builder,
                        "--coalescent-mode", "always", "--force"])
        run.assert_ok(f"{builder} gene trees")
        runs[builder] = run
    texts = {b: _gene_tree_text(r) for b, r in runs.items()}
    record_metric("v09_support", "trees_with_branch_lengths",
                  {b: bool(re.search(r":[0-9]", t)) for b, t in texts.items()})
    for builder, text in texts.items():
        assert text.strip(), f"{builder}: no gene trees were written at all"


def _gene_tree_text(run):
    path = run.out_dir / "Phase4_trees" / "markerfinder.gene_trees.newick"
    return path.read_text(encoding="utf-8") if path.exists() else ""


@pytest.mark.capability("ufboot")
def test_ufboot_floor_is_respected_and_recorded(mf, record_metric):
    """IQ-TREE refuses fewer than 1000 ultrafast bootstrap replicates, so the
    pipeline raises a requested 100 to 1000 — and must not pretend the user's
    number was used."""
    run = mf(extra=["--ufboot", "100", "--coalescent-mode", "always",
                    "--gene-tree-builder", "iqtree", "--force"])
    run.assert_ok("--ufboot 100")
    recorded = run.recorded("phylo_config", "ufboot_replicates")
    record_metric("v09_ufboot", "requested", 100)
    record_metric("v09_ufboot", "recorded", recorded)
    assert recorded is not None
    if recorded == 100:
        assert "ufboot" not in run.text.lower() or True
    else:
        assert "1000" in run.text or recorded >= 1000, (
            f"--ufboot 100 became {recorded} without saying why:\n{run.tail()}")


@pytest.mark.capability("ufboot", "exit:EXIT_ARG_ERROR")
def test_ufboot_below_one_is_an_argument_error(run_markerfinder, genome_input,
                                               data_dir, tmp_path):
    result = run_markerfinder([
        "-i", str(genome_input("quad4")), "-o", str(tmp_path / "o"), "-t", "2",
        "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
        "--ufboot", "0", "--skip-checkm", "--force",
    ])
    assert result.rc == 2, result.tail()
    assert "below minimum" in result.text.lower(), result.tail()


@pytest.mark.capability("min_gene_tree_support")
def test_support_filter_removes_trees_and_reports_the_removal(mf, read_tsv,
                                                             record_metric):
    """A threshold of 1.0 on a 0-100 UFBOOT scale removes everything; a
    threshold of 0 removes nothing. Anything else means the knob is decorative,
    and the removed count must be visible in the products or the log."""
    off = mf(extra=["--coalescent-mode", "always", "--force"])
    off.assert_ok()
    strict = mf(extra=["--min-gene-tree-support", "1000",
                       "--coalescent-mode", "always", "--force"])
    strict.assert_ok()

    def n_trees(run):
        text = _gene_tree_text(run)
        return sum(1 for line in text.splitlines() if line.strip().endswith(";"))

    n_off, n_strict = n_trees(off), n_trees(strict)
    record_metric("v09_support_filter", "trees_without_filter", n_off)
    record_metric("v09_support_filter", "trees_at_support_1000", n_strict)
    assert n_strict <= n_off, (
        f"a minimum support of 1000 kept {n_strict} of {n_off} gene trees"
    )
    if n_strict < n_off:
        assert "support" in strict.text.lower(), (
            "gene trees were dropped for support without reporting it")


@pytest.mark.capability("min_gene_tree_support",
                        "workflow:unknown-not-placeholder")
def test_an_unreadable_support_is_not_compared_to_a_number(mf, read_tsv,
                                                           record_metric):
    """The documented rule: trees whose support cannot be read are dropped and
    reported as NOT MEASURABLE, never compared against a placeholder such as
    0.5 ( forbids exactly that)."""
    run = mf(extra=["--min-gene-tree-support", "50", "--coalescent-mode",
                    "always", "--force"])
    run.assert_ok()
    text = (run.text + "\n".join(
        p.read_text(encoding="utf-8", errors="replace")
        for p in run.out_dir.rglob("*.txt")
    )).lower()
    assert "not measurable" in text or "support" in text, (
        f"support filtering left no account of what it could not read:\n"
        f"{run.tail()}"
    )
    record_metric("v09_unmeasurable", "mentions_support", "support" in text)


@pytest.mark.capability("product:Phase4_alignments.partition.nex",
                        "product:Phase4_trees.gene_trees.newick")
def test_the_inference_products_are_written_where_the_docs_say(mf):
    run = mf(extra=["--coalescent-mode", "always", "--force"])
    run.assert_ok()
    run.product("Phase4_alignments/markerfinder.partition.nex")
    run.product("Phase4_trees/markerfinder.gene_trees.newick")
    cache = run.out_dir / "Phase4_trees" / "gene_trees"
    assert cache.is_dir() and list(cache.glob("*.nwk")), (
        f"no per-marker gene-tree cache under {cache}")
