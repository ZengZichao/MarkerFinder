"""MarkerFinder command-line entry point: orchestration and dispatch.

This is the historical ``main`` function, relocated verbatim from
``markerfinder.__main__``. The only changes versus the original are:
  * the signal-handler / log-level bits are delegated to
:mod:`markerfinder.cli.logging_setup` (``resolve_log_level``);
  * the single-step subcommand branch is dispatched to the per-step
    ``execute_*`` helpers in:mod:`markerfinder.cli.commands`.

Everything else — argument parsing, validation order, taxonomy resolution,
config construction, the aggregate ``run`` branch, and all exit codes — is
unchanged.
"""

from __future__ import annotations

import logging
import os
from typing import List, Optional

import argparse

from markerfinder.phases import canonical_log_tag
from markerfinder.banner import print_banner
from markerfinder.exceptions import (
    ConfigError,
    CrossValidationError,
    MultiTreeError,
    PhyloToolError,
    TreeValidationError,
)
from markerfinder.pipeline import MarkerFinderPipeline
from markerfinder.utils.io import load_genomes_from_directory
from markerfinder.utils.logging_utils import setup_logging

from markerfinder.cli.parser import _build_parser, _apply_config_file
from markerfinder.cli.config_build import _build_pipeline_config
from markerfinder.cli.validation import (
    _validate_inputs,
    _check_output_conflicts,
    _handle_tree_input,
    _handle_cross_validation,
    _load_taxonomy,
    _auto_parse_embedded_taxonomy,
    _discover_bundled_hmm_dir,
)
from markerfinder.cli.self_test import _run_self_test
from markerfinder.cli.logging_setup import _register_signal_handlers, resolve_log_level
from markerfinder.cli.commands.scan import execute as execute_scan
from markerfinder.cli.commands.filter import execute as execute_filter
from markerfinder.cli.commands.infer import execute as execute_infer
from markerfinder.cli.commands.report import execute as execute_report
from markerfinder.cli.constants import (
    EXIT_SUCCESS,
    EXIT_RUNTIME_ERROR,
    EXIT_ARG_ERROR,
    EXIT_DATA_ERROR,
    EXIT_ASSERTION_FAILED,
    EXIT_INCONCLUSIVE,
)
from markerfinder.exceptions import AssertionFailureError

logger = logging.getLogger(__name__)


def main(argv: Optional[List[str]] = None) -> int:
    _register_signal_handlers()

    parser = _build_parser()

    try:
        args = parser.parse_args(argv)
    except argparse.ArgumentError as e:
        logger.error(f"Argument error: {e}")
        return EXIT_ARG_ERROR

    if args.self_test:
        # Honour an explicit --db-dir. The default is left to the
        # Self-test (repository db/), so this only ever narrows the check to the
        # Directory the user actually named; --check runs before config-file
        # Handling, so a db_dir coming from --config is not in play here.
        explicit_db = getattr(args, "db_dir", None)
        if explicit_db in (None, parser.get_default("db_dir")):
            explicit_db = None
        return _run_self_test(explicit_db)

    if args.config:
        try:
            args = _apply_config_file(args, parser)
        except ConfigError as e:
            logger.error(str(e))
            return EXIT_ARG_ERROR

    command = getattr(args, "command", None)  # None → full 'run' aggregate

    # Required-arg validation differs by command. Must run AFTER config file
    # Loading so that ``input``/``output`` values supplied only in the config
    # File are honored.
    if command is None:
        # Legacy aggregate 'run': both -i and -o required
        if not args.input or not args.output:
            parser.error("the following arguments are required: -i/--input, -o/--output")
    else:
        # Subcommand mode: every step needs -o (locates pipeline state);
        # Only 'scan' additionally requires -i
        if not args.output:
            parser.error("-o/--output is required (locates pipeline state)")
        if command == "scan" and not args.input:
            parser.error("-i/--input is required for the 'scan' subcommand")

    log_level = resolve_log_level(args)

    setup_logging(level=log_level, log_file=args.log_file)

    print_banner()

    try:
        _validate_inputs(args)
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else EXIT_ARG_ERROR

    _check_output_conflicts(args)

    try:
        tree = _handle_tree_input(args)
    except MultiTreeError as e:
        logger.error(str(e))
        return EXIT_ARG_ERROR
    except (TreeValidationError, FileNotFoundError) as e:
        logger.error(f"Tree validation failed: {e}")
        return EXIT_DATA_ERROR

    try:
        _handle_cross_validation(args)
    except CrossValidationError as e:
        logger.error(str(e))
        return EXIT_DATA_ERROR

    # ``table_taxa`` drives the phylogenetic HGT step; track *how* it was
    # Produced so Phase 0.2 can log the real origin instead of always saying
    # "external table". Parallel to ``table_taxa`` (None/NonNone shapes match).
    taxonomy_origin: str = "external table"

    try:
        table_taxa = _load_taxonomy(args)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Failed to load taxonomy table: {e}")
        return EXIT_DATA_ERROR

    # ── Resolve taxonomy for the phylogenetic HGT step ──
    # Order of precedence (highest wins):
    # 1. User-supplied --taxonomy-table file (already loaded into table_taxa).
    # 2. --taxonomy-format embedded → parse Format A out of genome ids/filenames.
    # 3. Auto-detect (hmm mode only): in hmm mode markers are identified from
    # TIGRFAM/Pfam HMMs per genome, so the only way the phylogenetic HGT
    # Step can run is the MAD + monophyly path, which needs a taxonomy map.
    # If the genome identifiers carry Format A taxonomy and the user has not
    # Supplied one, auto-enable embedded parsing — otherwise the phylogenetic
    # Step stays silent and every marker receives an unknown phylogenetic
    # Risk (defaulting to conservative Level-3), ending up as all-Level-3
    # With no species tree.
    if table_taxa is None and getattr(args, "taxonomy_format", "table") == "embedded":
        # Explicit opt-in via --taxonomy-format embedded.
        table_taxa = _auto_parse_embedded_taxonomy(args, force=True)
        if table_taxa:
            taxonomy_origin = "embedded Format A labels (--taxonomy-format embedded)"
    elif table_taxa is None and args.auto_embed_taxonomy and \
            getattr(args, "marker_mode", "hmm") == "hmm":
        # Silent auto-detect: only kick in if the genome ids actually look like
        # Format A labels in the majority, so we don't impose structure on
        # Arbitrary identifiers.
        table_taxa = _auto_parse_embedded_taxonomy(args, force=False)
        if table_taxa:
            taxonomy_origin = "embedded Format A labels (auto-detected from genome ids)"

    # In any marker mode, if after all of the above we still have no taxonomy map
    # And no species tree, the phylogenetic HGT step will not run. Warn the user
    # Loudly so the "all-UNKNOWN" outcome is explainable and actionable. UNKNOWN
    # Markers are KEPT for tree inference (never silently dropped as Level 3),
    # But their HGT status is unscreened and the report flags them as such.
    has_reference = bool(getattr(args, "species_tree", None) or getattr(args, "tree", None))
    if table_taxa is None and not has_reference:
        logger.warning(
            f"{canonical_log_tag('0.2')} No taxonomy available for the phylogenetic HGT step "
            "(no --taxonomy-table, --species-tree, or auto-detectable Format A labels). "
            "The HGT phylogenetic step (MAD-rooting + monophyly) will be skipped, so every "
            "marker's HGT status will be UNKNOWN: markers are kept for tree inference, "
            "but they were NOT screened for HGT (the report lists them under 'unknown'). "
            "Pass --taxonomy-table, or if your genome ids carry Format A taxonomy "
            "('_d_<domain>_p_<phylum>_...'), pass --taxonomy-format embedded (the default "
            "auto-detect is on)."
        )

    config = _build_pipeline_config(args, parser=parser)

    # In hmm mode with no explicit --marker-hmm-dir, auto-discover the bundled
    # TIGRFAM/Pfam HMM library under --db-dir. The directory is chosen by the
    # Detected domain of the input genomes (ar53 for Archaea, bac120 for Bacteria)
    # So the marker set matches the clade. Performed BEFORE constructing the
    # Pipeline so AdaptiveMarkerSelectionModule picks up the resolved path at
    # Init time. Skipped when the user passed --marker-hmm-dir, or when
    # --marker-db-source is "none".
    if (
        config.selection_config.marker_mode == "hmm"
        and not config.selection_config.marker_hmm_dir
        and (config.selection_config.marker_db_source or "auto") != "none"
    ):
        resolved = _discover_bundled_hmm_dir(args, config, table_taxa=table_taxa)
        if resolved:
            config.selection_config.marker_hmm_dir = resolved
            logger.info(
                f"{canonical_log_tag('1')} Auto-discovered bundled HMM library for hmm mode: {resolved} "
                f"(marker_db_source={config.selection_config.marker_db_source})"
            )
        else:
            logger.warning(
                f"{canonical_log_tag('1')} hmm mode: no bundled TIGRFAM/Pfam HMM library discovered under "
                f"--db-dir '{config.db_dir}' (marker_db_source={config.selection_config.marker_db_source}). "
                "Pass --marker-hmm-dir to point at the per-marker HMM directory."
            )

    pipeline = MarkerFinderPipeline(config)

    # ── Full aggregate run (no subcommand) ──
    if command is None:
        genomes = load_genomes_from_directory(args.input)
        if not genomes:
            logger.error(f"No genomes found in {args.input}")
            return EXIT_DATA_ERROR
        try:
            result = pipeline.run(
                genomes,
                reference_tree=tree,
                table_taxa=table_taxa,
                taxonomy_origin=taxonomy_origin,
            )
        except AssertionFailureError as e:
            # An implausible output or a failed
            # Taxonomy must-pass gate is exit 4, never 3. The stepwise route
            # Already mapped it correctly; the aggregate route caught it as the
            # PhyloToolError it subclasses, so the same failure answered with a
            # Different code depending on how the pipeline was invoked — and
            # "output unreasonable" became indistinguishable from "input bad".
            logger.error(f"Assertion failure: {e}")
            return EXIT_ASSERTION_FAILED
        except PhyloToolError as e:
            logger.error(f"Pipeline error: {e}")
            return EXIT_DATA_ERROR
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            if args.verbose >= 2 or "MARKERFINDER_DEBUG" in os.environ:
                import traceback
                traceback.print_exc()
            return EXIT_RUNTIME_ERROR
        if result is None:
            logger.error("Pipeline returned no result")
            return EXIT_RUNTIME_ERROR

        # ── 可用性汇总检验 ──
        # 管线可能在无真正结果的情况下"完成"(0 误): 典型情况是 HGT 过滤把所有
        # 标记判为 等级 3(组成成分步骤饱和 + 系统发育步骤未运行), 或是建树材料不足
        # 导致没有物种树. 这种情况不能当作成功; 要带强错误返回 EXIT_DATA_ERROR,
        # 让用户明确看到管线没有产出可用的物种树/标记集.
        if getattr(result, "usable_marker_count", 0) <= 0:
            # The refusal must name what actually happened. "all markers were
            # Flagged as Level 3" is wrong when Phase 1 selected nothing (the
            # Screen had no marker to grade), and a reader who follows that
            # Advice would loosen thresholds that were never the problem.
            hgt_report = getattr(result, "hgt_report", None)
            graded = getattr(hgt_report, "total_markers", 0) or 0
            if graded == 0:
                logger.error(
                    "Pipeline completed with NO MARKERS REACHING THE HGT SCREEN: "
                    "Phase 1 selected 0 markers, so nothing was graded and no "
                    "species tree / marker set exists. Check that the marker "
                    "directory holds files whose headers carry these genome ids, "
                    "and that --min-occupancy / --min-hmm-score / --max-markers "
                    "are not excluding every candidate."
                )
            else:
                logger.error(
                    f"Pipeline completed but has NO USABLE MARKERS after HGT "
                    f"filtering: {graded} marker(s) were graded and "
                    f"{getattr(hgt_report, 'level3_count', 0)} were excluded as "
                    f"Level 3 (Level 1: {getattr(hgt_report, 'level1_count', 0)}, "
                    f"Level 2: {getattr(hgt_report, 'level2_count', 0)}, "
                    f"UNKNOWN: {getattr(hgt_report, 'unknown_count', 0)}). "
                    "No species tree / marker set was produced. If the "
                    "phylogenetic step never activated (no --taxonomy-table, "
                    "--species-tree, or detectable Format A labels) the "
                    "composition step alone can saturate; check the HGT columns "
                    "of Phase5_reports/markerfinder.hgt_evaluation.tsv for the "
                    "evidence each marker was graded on."
                )
            return EXIT_DATA_ERROR
        if getattr(result, "species_tree_source", "none") == "none":
            logger.error(
                "Pipeline completed but produced NO SPECIES TREE (species_tree_source=none). "
                "Usable markers: %d. A concatenation (supermatrix) or coalescent tree could "
                "not be built — check that enough markers passed the HGT filter and a tree "
                "is inferrable from the dataset.",
                getattr(result, "usable_marker_count", 0),
            )
            return EXIT_DATA_ERROR
        # Exit 5 is for a run that finished and built trees but
        # Whose recommendation layer could not name a resolved one.
        from markerfinder.cli.constants import recommendation_is_inconclusive

        if recommendation_is_inconclusive(
            getattr(getattr(result, "phylogenetic", None),
                    "tree_recommendation", None)
        ):
            logger.error(
                "Pipeline completed, but the recommendation layer returned an "
                "INCONCLUSIVE verdict (no resolved tree). This is exit code 5, "
                "not a success: see the Tree recommendation line in "
                "Phase5_reports/{prefix}.pipeline_summary.txt."
            )
            return EXIT_INCONCLUSIVE
        return EXIT_SUCCESS

    # ── Single-step subcommand (scan/filter/infer/report) ──
    # 'scan' loads genomes from -i and may carry a reference species tree + taxonomy table.
    # 'filter'/'infer'/'report' load all prior state from
    # <output>/.markerfinder/.pipeline_state.json.
    if command == "scan":
        return execute_scan(args, pipeline, tree, table_taxa, taxonomy_origin, logger)
    if command == "filter":
        return execute_filter(args, pipeline, tree, table_taxa, taxonomy_origin, logger)
    if command == "infer":
        return execute_infer(args, pipeline, tree, table_taxa, taxonomy_origin, logger)
    if command == "report":
        return execute_report(args, pipeline, tree, table_taxa, taxonomy_origin, logger)
    return EXIT_ARG_ERROR
