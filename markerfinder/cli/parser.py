"""Argument parser definitions for the MarkerFinder CLI.

This module owns every argparse artefact (custom types, shared/common args,
subcommand metadata, and the top-level ``_build_parser``). It is a verbatim
relocation of the historical definitions from ``markerfinder.__main__``; the
public ``--help`` text, parameter order, defaults, and exit codes are
unchanged.
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Optional

from markerfinder._version import get_version_string
from markerfinder.config_loader import load_config_file, merge_config_with_args
from markerfinder.exceptions import ConfigError

logger = logging.getLogger(__name__)


class PathType:
    """Argparse type for path validation."""

    def __init__(self, exists: bool = True, is_file: bool = True, is_dir: bool = False):
        self.exists = exists
        self.is_file = is_file
        self.is_dir = is_dir

    def __call__(self, value: str) -> str:
        p = Path(value)
        if self.exists and not p.exists():
            raise argparse.ArgumentTypeError(f"Path does not exist: {value}")
        # Is_file / is_dir are not mutually exclusive in the constructor, so gate
        # The "expected a file" check on is_dir being unset — otherwise
        # PathType(exists=True, is_dir=True) (used by -i/--input) would wrongly
        # Reject every real directory.
        if self.is_file and not self.is_dir and p.exists() and p.is_dir():
            raise argparse.ArgumentTypeError(f"Expected a file, got directory: {value}")
        if self.is_dir and p.exists() and p.is_file():
            raise argparse.ArgumentTypeError(f"Expected a directory, got file: {value}")
        return value


class DeprecatedAliasAction(argparse._StoreAction):
    """Store into the canonical dest and remember that the alias was used.

    ``--tree`` shares its dest with ``--species-tree``, so after parsing there is
    no trace of which of the two the user typed — and a deprecation notice that
    checks ``args.tree`` can therefore never fire. This action records the fact
    under its own attribute so the notice is about the option actually used.
    """

    def __init__(self, option_strings, dest, alias_attr, **kwargs):
        self.alias_attr = alias_attr
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        setattr(namespace, self.alias_attr, True)


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add shared pipeline arguments to *parser*.

    Called twice: once on the top-level parser (real defaults, so the bare
    legacy form ``markerfinder -i X -o Y`` works) and once on a throwaway
    template whose defaults are then replaced with ``argparse.SUPPRESS`` so
    each subcommand inherits only the options the user actually supplied —
    everything else is inherited from the top-level parser.
    """
    parser.add_argument("-i", "--input",
                        type=PathType(exists=True, is_dir=True),
                        help="Input genome directory (protein FASTA: .faa/.fasta/.fa; "
                             "matching .fna files are optional and only used by CheckM)")
    parser.add_argument("-o", "--output",
                        help="Output directory (created if not exists)")
    parser.add_argument("-t", "--threads", type=IntRange(min_val=1),
                        default=min(8, os.cpu_count() or 1),
                        help="Number of threads (default: min(8, CPU count))")

    parser.add_argument("--mode",
                        choices=["conservative", "standard", "expanded", "mag_adaptive"],
                        default="standard",
                        help="Analysis resolution mode (default: standard)")

    tree_group = parser.add_argument_group("Tree & Taxonomy Options")
    tree_group.add_argument("--tree",
                            type=PathType(exists=True, is_file=True),
                            dest="species_tree",
                            action=DeprecatedAliasAction,
                            alias_attr="used_tree_alias",
                            default=None,
                            help="Deprecated alias for --species-tree. Input reference species tree file (Newick/Nexus/NHX)")
    tree_group.add_argument("--sequences",
                            type=PathType(exists=True, is_file=True),
                            help="Input sequence file (FASTA/FASTQ)")
    tree_group.add_argument("--taxonomy-table",
                            type=PathType(exists=True, is_file=True),
                            help="External taxonomy table file")
    tree_group.add_argument("--taxonomy-format",
                            choices=["table", "embedded"],
                            default="table",
                            help="Taxonomy table format: table (Format B d__X;p__Y) or embedded (Format A _d_X_p_Y)")
    tree_group.add_argument("--taxonomy-source-priority",
                            choices=["embedded", "table"],
                            default="table",
                            help="Priority when both embedded and table taxonomy exist (default: table)")
    tree_group.add_argument("--taxonomy-delimiter-mode",
                            choices=["reverse", "greedy", "segment"],
                            default="reverse",
                            help="Parsing strategy for Format A taxonomy (default: reverse)")
    tree_group.add_argument("--table-sep",
                            help="Force table separator (auto-detect if not set)")
    tree_group.add_argument("--no-auto-embed",
                            dest="auto_embed_taxonomy",
                            action="store_false",
                            default=True,
                            help="Disable auto-detection of Format A (e.g. _d_Archaea_p_...) "
                                 "embedded taxonomy from genome file/.faa identifiers. "
                                 "Auto-detection only applies in --marker-mode hmm with no "
                                 "explicit --taxonomy-table / --species-tree, where an absent "
                                 "taxonomy would otherwise leave the phylogenetic HGT step "
                                 "silent and flag all markers as HGT (Level 3).")
    tree_group.add_argument("--multi-tree-mode",
                            choices=["ask", "split", "first", "last", "random"],
                            default="ask",
                            help="How to handle multiple trees in one file (default: ask). "
                                 "first/last/random select one tree from the file; "
                                 "split is not supported for the reference species tree "
                                 "and fails with an explicit error")
    tree_group.add_argument("--strip-annotations", action="store_true",
                            help="Remove NHX/bootstrap annotations from tree")
    tree_group.add_argument("--mol-type",
                            choices=["DNA", "RNA", "protein"],
                            help="Force molecule type for sequence validation (auto-detect if not set)")
    tree_group.add_argument("--skip-length-check", action="store_true",
                            help="Skip sequence length consistency check")
    tree_group.add_argument("--no-cross-check", action="store_true",
                            help="Skip tree-sequence cross-validation (not recommended)")
    tree_group.add_argument("--ignore-malformed", action="store_true",
                            help="Skip malformed rows/files instead of terminating (default: terminate)")
    tree_group.add_argument("--taxonomy-levels",
                            help="Extend the built-in taxonomic ranks with custom ones, as "
                                 "'level:prefix' pairs (e.g. 'kingdom:k__' for Format B or "
                                 "'kingdom:_k_' for Format A). The rank is then parsed out of "
                                 "--taxonomy-table / embedded labels instead of being dropped "
                                 "with an 'Unknown level prefix' warning")

    hgt_group = parser.add_argument_group("HGT Filtering Options")
    hgt_group.add_argument("--hgt-steps",
                           type=str,
                           default="phylogenetic",
                           help="Which HGT screening steps to enable. 'phylogenetic' = "
                                "gene-tree vs species-tree incongruence, or MAD rooting + "
                                "monophyly proportion with a taxonomy table. "
                                "'composition' = RCV/GC composition-bias diagnostics, "
                                "written as PARALLEL evidence columns only (never merged "
                                "into the risk score). Default: 'phylogenetic'")
    hgt_group.add_argument("--hgt-threshold", type=float, default=None,
                           help="DEPRECATED alias for --hgt-threshold-l1l2 "
                                "(0.25 when neither is given). It only moves the "
                                "Level1/Level2 band, never Level2/Level3; prefer "
                                "--hgt-threshold-l1l2 / --hgt-threshold-l2l3.")
    # Split band thresholds. --hgt-threshold stays as a
    # Deprecated alias of --hgt-threshold-l1l2 (one deprecation notice).
    hgt_group.add_argument("--hgt-threshold-l1l2", type=float, default=None,
                           help="Level1/Level2 boundary (default: 0.25). "
                                "Overrides --hgt-threshold.")
    hgt_group.add_argument("--hgt-threshold-l2l3", type=float, default=None,
                           help="Level2/Level3 boundary (default: 0.60).")
    # Marker budget, previously hard-coded at 60.
    marker_group_kwargs = dict()
    hgt_group.add_argument("--max-markers", type=int, default=None,
                           help="Maximum number of markers selected (default: 60, "
                                "or the Phase-0 adaptive budget when the analysis "
                                "mode sets one). Giving this flag explicitly pins "
                                "the budget: it then overrides the adaptive value "
                                "instead of being overridden by it.")
    hgt_group.add_argument("--marker-preset", dest="marker_preset",
                              choices=["none", "conservative", "standard", "expanded", "mag_adaptive"],
                              default="none",
                              help="Named resolution preset (D-02). 'none' (default) keeps "
                                   "current behaviour; an explicit preset overrides "
                                   "min_occupancy/max_markers/selection strategy.")
    # Threshold scan.
    hgt_group.add_argument("--hgt-scan", action="store_true", default=False,
                           help="Scan threshold bands (default grid 0.15/0.25/0.40/0.60 "
                                "+ far band 0.25/0.95) and write "
                                "<output>/Phase5_reports/<prefix>.threshold_scan.tsv.")
    hgt_group.add_argument("--scan-stability-min", type=float, default=0.6,
                           help="Minimum acceptable cross-band Jaccard of selected "
                                "marker sets (default: 0.6).")
    hgt_group.add_argument("--require-evidence-coverage", type=float, default=0.5,
                           help="Warn prominently when HGT evidence coverage falls "
                                "below this fraction (default: 0.5).")
    # Criterion mode. Default 'risk' never auto-switches
    #; consistency/hybrid are gated.
    hgt_group.add_argument("--hgt-mode", choices=["risk", "consistency", "hybrid"],
                           default="risk",
                           help="Marker screening criterion (default: risk). "
                                "consistency: cross-framework grade decides. "
                                "hybrid: risk level AND consistency grade must "
                                "both pass. consistency/hybrid "
                                "require the consistency gate to pass.")
    hgt_group.add_argument("--consistency-stringency", type=int, choices=range(1, 6),
                           default=1,
                           help="Consistency stringency ladder 1..5, descending thresholds. Default: 1.")
    # Provisional; 0 disables the informativeness screen.
    hgt_group.add_argument("--min-informative-sites", type=int, default=0,
                           help="Minimum parsimony-informative sites per marker "
                                "(provisional; 0 = disabled). Default: 0.")
    # Optional COG -> functional category map feeding the
    # Functional_category column. Absent file => that column renders NA with a
    # Report note, instead of being silently blank.
    hgt_group.add_argument("--cog-category-map", dest="cog_category_map",
                           default=None, metavar="TSV",
                           help="Optional TSV of marker_id<TAB>functional_category "
                                "(D-04). Default: unset, column renders NA.")
    assertion_group = parser.add_argument_group("Self-Check (Assertion) Options")
    assertion_group.add_argument("--strict-assertions", dest="strict_assertions",
                                 action="store_true", default=True,
                                 help="Fail the run (exit 4) when a severity=fail output "
                                      "assertion fires. Default: on.")
    assertion_group.add_argument("--allow-assertion-failure",
                                 dest="allow_assertion_failure",
                                 action="store_true", default=False,
                                 help="Explicitly downgrade failing output assertions to a "
                                      "recorded warning instead of exiting with code 4. "
                                      "The bypass is written to the report. Default: off.")
    # Taxonomy must-pass gate. Opt-in, because real HGT
    # Markers legitimately break taxon monophyly; but when it IS asked for and
    # The baseline cannot be loaded, the run refuses rather than passing.
    from markerfinder.cli.constants import default_taxonomy_mustpass_path

    assertion_group.add_argument("--taxonomy-mustpass",
                                 dest="taxonomy_mustpass",
                                 nargs="?", const=default_taxonomy_mustpass_path(),
                                 default=None, metavar="YAML",
                                 help="Abort (exit 4) unless every gene tree satisfies the "
                                      "required taxonomic relationships in the baseline YAML "
                                      "(default: tests/benchmark/expected/"
                                      "taxonomy_mustpass.yaml). Unset = gate off, and the run "
                                      "says so explicitly.")

    marker_group = parser.add_argument_group("Marker Discovery Options")
    marker_group.add_argument("--marker-mode",
                              choices=["gtdb_tk", "hmm"],
                              default="gtdb_tk",
                              help="Marker discovery mode: 'gtdb_tk' consumes GTDB-TK ar53/bac120 per-marker "
                                   "raw/unaligned FASTA output (set --gtdb-markers-dir); 'hmm' uses HMMER "
                                   "to scan each input genome against a directory of TIGRFAM/Pfam HMM profiles "
                                   "(--marker-hmm-dir) — identifies markers internally and extracts the matched "
                                   "protein sequences, so no GTDB-TK pre-run is required (default: gtdb_tk)")
    marker_group.add_argument("--marker-hmm-dir",
                              help="[hmm mode] Directory of TIGRFAM/Pfam HMM profile files (one .HMM/.hmm per "
                                   "marker, filename stem = marker id). Used with --marker-mode hmm. Example: "
                                   "GTDB-Tk-214-Markers/ containing pfam/individual_hmms/ and tigrfam/individual_hmms/. "
                                   "When omitted in hmm mode the bundled library under --db-dir is auto-discovered "
                                   "according to --marker-db-source.")
    marker_group.add_argument("--marker-db-source",
                              choices=["auto", "cog", "gtdb", "none"],
                              default="auto",
                              help="[hmm mode, auto-discovery] How to locate the bundled TIGRFAM/Pfam HMM library "
                                   "when --marker-hmm-dir is not given: 'auto' pick ar53/bac120 under --db-dir "
                                   "by detected domain; 'cog' use only the cog_hmm tree; 'gtdb' use only "
                                   "gtdb_markers (ar53/bac120); 'none' disable auto-discovery so "
                                   "--marker-hmm-dir is required (default: auto)")
    marker_group.add_argument("--min-hmm-score", type=float, default=20.0,
                              help="[hmm mode] Minimum bitscore for an hmmsearch hit to count as a marker "
                                   "hit (default: 20.0). Increase to tighten marker identification.")
    marker_group.add_argument("--gtdb-markers-dir",
                              help="[gtdb_tk mode] Directory of GTDB-TK extracted per-marker raw/unaligned FASTA files.")
    marker_group.add_argument("--species-tree",
                              help="Reference species tree (Newick) used as "
                                   "the reference for HGT phylogenetic-step screening (RF/quartet vs "
                                   "this tree). Alternative to --taxonomy-table; if a taxonomy table is "
                                   "given, the phylogenetic step uses MAD rooting + monophyly proportion instead. "
                                   "--tree is accepted as a deprecated alias.")
    marker_group.add_argument("--monophyly-rank",
                              choices=["auto", "domain", "phylum", "class", "order",
                                       "family", "genus", "species"],
                              default="auto",
                              help="[--taxonomy-table] Taxonomic rank at which the gene-tree "
                                   "monophyly proportion is measured during the phylogenetic "
                                   "step with a taxonomy table (default: auto = one rank below "
                                   "the taxonomic scope detected from the tree's tips; naming a "
                                   "rank measures at exactly that rank, falling back upward "
                                   "only when no taxon there offers an informative split, which "
                                   "is then reported as a cross-rank comparison)")
    marker_group.add_argument("--monophyly-threshold", type=float, default=0.5,
                              help="[--taxonomy-table] Minimum monophyly proportion for a marker "
                                   "to be considered phylogenetically clean (default: 0.5). Below this, the "
                                   "marker is flagged as HGT-prone (phylogenetic-step risk = 1 - proportion).")
    marker_group.add_argument("--hgt-adaptive-thresholds", action="store_true", default=False,
                              help="(default: Off) Opt in to auto-relax the Level-2/Level-3 risk threshold on "
                                   "phylogenetically distant datasets (multiple genera/orders) where phylogenetic HGT "
                                   "risk saturates, so markers are not all discarded as Level 3. Disabled by default "
                                   "so the HGT放行尺度 is never silently relaxed; pass this flag to enable. "
                                   "Equivalent no-op disable: --no-hgt-adaptive-thresholds.")
    marker_group.add_argument("--no-hgt-adaptive-thresholds", action="store_true", default=False,
                              help="Explicitly disable --hgt-adaptive-thresholds (the default behaviour).")

    phylo_group = parser.add_argument_group("Phylogenetic Inference Options")
    phylo_group.add_argument("--ufboot", type=IntRange(min_val=1), default=1000,
                             help="UFBOOT replicates (default: 1000)")
    phylo_group.add_argument("--fast-tree", action="store_true",
                             help="Use FastTree2 instead of IQ-TREE3 for gene trees (coalescent mode). "
                                  "Deprecated alias for --gene-tree-builder fasttree; kept for backward compatibility.")
    phylo_group.add_argument("--gene-tree-builder",
                             choices=["fasttree", "iqtree"],
                             default=None,
                             help="Gene-tree builder for coalescent inference: 'fasttree' (FastTree2 with WAG, fast) "
                                  "or 'iqtree' (IQ-TREE3, more thorough). Default: fasttree. "
                                  "This option explicitly selects the per-marker gene-tree method; "
                                  "it overrides --fast-tree and config.fast_tree.")
    phylo_group.add_argument("--coalescent-mode",
                             choices=["off", "post-filter", "always"],
                             default="post-filter",
                             help="Coalescent mode (default: post-filter)")
    phylo_group.add_argument("--min-gene-tree-support", type=float, default=None,
                             help="Minimum average gene-tree local support. Gene trees below "
                                  "this threshold are filtered out of coalescent species-tree "
                                  "inference. Units must match the support scale of the trees: "
                                  "IQ-TREE UFBOOT/SH labels are 0-100, FastTree SH and ASTRAL "
                                  "local support are 0-1. Trees whose support cannot be read "
                                  "are dropped and reported as NOT MEASURABLE, never compared "
                                  "against a placeholder. Default: None (disabled; only the "
                                  "n_tips >= 4 filter is applied).")

    mag_group = parser.add_argument_group("MAG Optimization Options")
    mag_group.add_argument("--checkm-results",
                           help="Pre-computed CheckM results file")
    mag_group.add_argument("--skip-checkm", action="store_true",
                           help="Skip CheckM quality assessment and use default quality estimates "
                                "(useful for high-quality complete genomes or to speed up testing)")
    mag_group.add_argument("--min-occupancy", type=float, default=None,
                           help="Explicit minimum marker occupancy floor (0-1). When set, "
                                "overrides the Phase-0 adaptive occupancy threshold for "
                                "marker selection in all modes")
    mag_group.add_argument("--min-marker-coverage", type=float, default=None,
                           help="Deprecated alias of --min-occupancy (kept for compatibility; "
                                "the value is applied to marker-selection occupancy)")

    report_group = parser.add_argument_group("Report Options")
    report_group.add_argument("--report-format",
                              choices=["html"],
                              default="html",
                              help="Report format (only HTML is supported)")

    output_group = parser.add_argument_group("Output Management")
    output_group.add_argument("--force", action="store_true",
                              help="Overwrite existing output files")
    output_group.add_argument("--no-clobber", action="store_true",
                              help="Skip output if file already exists (do not overwrite)")
    output_group.add_argument("--save-intermediates", action="store_true",
                              help="Copy useful intermediate files (marker sequences, alignments, "
                                   "supermatrix, CheckM output) to <output>/Phase4_intermediate/")
    output_group.add_argument("--redo", action="store_true",
                              help="Re-run the requested step even if it was already completed")
    output_group.add_argument("--resume", action="store_true",
                              help="Automatically run any missing prerequisite steps up to the "
                                   "requested step")

    config_group = parser.add_argument_group("Configuration")
    config_group.add_argument("--config",
                              type=PathType(exists=True, is_file=True),
                              help="Path to YAML/TOML configuration file")
    config_group.add_argument("--log-file",
                              help="Path to log file (UTF-8, persistent, real-time flush)")

    adv_group = parser.add_argument_group("Advanced Options")
    adv_group.add_argument("--tmp-dir", default=None,
                           help="Temporary file directory (default: a per-run isolated "
                                "directory under the system temp dir, auto-cleaned on exit)")
    adv_group.add_argument("--keep-tmp", action="store_true",
                           help="Keep the temporary directory after pipeline completion (default: delete)")
    adv_group.add_argument("--db-dir", default="./db",
                           help="Database directory")
    adv_group.add_argument("--verbose", "-v", action="count", default=0,
                           help="Verbose output (stackable: -v, -vv)")


class IntRange:
    """Argparse type for integer range validation."""

    def __init__(self, min_val: Optional[int] = None, max_val: Optional[int] = None):
        self.min_val = min_val
        self.max_val = max_val

    def __call__(self, value: str) -> int:
        iv = int(value)
        if self.min_val is not None and iv < self.min_val:
            raise argparse.ArgumentTypeError(
                f"Value {iv} below minimum {self.min_val}"
            )
        if self.max_val is not None and iv > self.max_val:
            raise argparse.ArgumentTypeError(
                f"Value {iv} above maximum {self.max_val}"
            )
        return iv


_SUBCOMMAND_SPECS = [
    ("scan", "Step 1: quality preprocessing + marker scanning + sequence extraction (Phase 0+1+1.5)",
     "Phase 0+1+1.5: CheckM-style quality preprocessing, adaptive marker "
     "scanning, and marker sequence extraction. Requires -i/--input. "
     "Persists state to <output>/.markerfinder/.pipeline_state.json for the next step."),
    ("filter", "Step 2: HGT-aware marker filtering (Phase 2)",
     "Phase 2: runs the phylogenetic HGT screening step "
     "and partitions markers into quality levels. Requires prior 'scan' step. "
     "Gene trees built here are cached "
     "under <output>/Phase4_trees/gene_trees/ (historical directory name, "
     "D-05) for reuse by 'infer'."),
    ("infer", "Step 3: phylogenetic inference - concatenation + coalescent (Phase 3)",
     "Phase 3: supermatrix + coalescent species-tree inference. Requires a "
     "prior 'filter' step. Reuses cached gene trees from 'filter' when available."),
    ("report", "Step 4: report generation (Phase 4)",
     "Phase 4: generates HTML/text reports under <output>/Phase5_reports/, species trees "
     "under <output>/Phase4_trees/, partition files under <output>/Phase4_alignments/, plus "
     "Phase5_metadata/run_config.json (Phase5_*/Phase4_* are historical directory "
     "names kept verbatim, D-05). Requires a prior 'infer' step."),
]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="markerfinder",
        description="MarkerFinder: Adaptive HGT-Aware Phylogenomic Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # Full pipeline in one shot (legacy aggregate 'run')
  markerfinder -i genomes/ -o output/ -t 8 --mode standard --gtdb-markers-dir gtdbtk_out/align/marker_genes
  markerfinder -i mags/ -o output/ -t 8 --mode mag_adaptive --checkm-results checkm.tsv \
      --gtdb-markers-dir gtdbtk_out/align/marker_genes
  markerfinder -i genomes/ -o output/ --species-tree tree.nwk --taxonomy-table taxa.tsv \
      --gtdb-markers-dir gtdbtk_out/align/marker_genes
  markerfinder -i genomes/ -o output/ --config config.yaml

  # Step-by-step (state persists under <output>/.markerfinder/ between steps)
  markerfinder scan     -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes
  markerfinder filter             -o output/ -t 8
  markerfinder infer              -o output/ -t 8 --gene-tree-builder fasttree
  markerfinder report            -o output/

  # MarkerFinder consumes GTDB-TK ar53/bac120 per-marker raw/unaligned FASTA;
  # MAFFT alignment + trimAl trimming are performed internally before tree inference.
  markerfinder -i genomes/ -o output/ -t 8 \
      --gtdb-markers-dir gtdbtk_out/align/marker_genes \
      --species-tree gtdbtk_out/classify/gtdbtk.bac120.classify.tree

  markerfinder --version
  markerfinder --check
""",
    )

    # Top-level parser carries the common args WITH real defaults so the bare
    # Legacy form (no subcommand) still works exactly as before.
    _add_common_args(parser)

    # Top-level-only flags (not inherited by subcommands).
    adv_group = parser.add_argument_group("Advanced Options")
    adv_group.add_argument("--version", action="version",
                           version=f"%(prog)s {get_version_string()}")
    adv_group.add_argument("--check", "--self-test", action="store_true",
                           dest="self_test",
                           help="Run self-test: verify dependencies and built-in examples")

    # Shared subcommand template: identical common args but with SUPPRESS
    # Defaults so a subcommand only overrides attributes the user actually
    # Passed — everything else is inherited from the top-level parser. This
    # Makes both `markerfinder scan -i X -o Y` (subcommand first) and
    # `markerfinder -i X -o Y scan` (subcommand last) work identically.
    common = argparse.ArgumentParser(add_help=False)
    _add_common_args(common)
    for _action in common._actions:
        _action.default = argparse.SUPPRESS

    # ── Subcommands (optional; omitted → full 'run' aggregate for backward compat) ──
    subparsers = parser.add_subparsers(
        dest="command",
        title="subcommands (optional)",
        description=(
            "Run a single pipeline step. State is persisted under "
            "<output>/.markerfinder/ (.pipeline_state.json) and "
            "<output>/Phase5_metadata/context.json so steps can be chained. "
            "If no subcommand is given, the full pipeline runs as one aggregate 'run'."
        ),
    )
    for name, help_text, desc_text in _SUBCOMMAND_SPECS:
        subparsers.add_parser(
            name,
            parents=[common],
            help=help_text,
            description=desc_text,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    return parser


def _apply_config_file(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> argparse.Namespace:
    """Load config file and merge with CLI args (CLI > Config > Defaults).

    Args:
        args: Parsed CLI arguments.
        parser: The argument parser used to resolve default values.
    """
    if not args.config:
        return args

    try:
        config_data = load_config_file(args.config)
    except ConfigError:
        raise
    except Exception as e:
        raise ConfigError(f"Failed to load config file '{args.config}': {e}")

    args_dict = {k: v for k, v in vars(args).items() if v is not None}
    defaults = {k: parser.get_default(k) for k in vars(args)}
    merged = merge_config_with_args(config_data, args_dict, defaults)

    # 配置文件里出现的未知键只警告不报错: 静默吞掉会让用户以为生效了
    # (例如 MANUAL 曾推荐 adaptive_far_thresholds 这个并不存在的键).
    known_keys = set(vars(args).keys())
    unknown_keys = [k for k in merged if k not in known_keys]
    if unknown_keys:
        logger.warning(
            f"Config file '{args.config}' contains unrecognized keys "
            f"(ignored): {sorted(unknown_keys)}"
        )

    for k, v in merged.items():
        setattr(args, k, v)

    return args
