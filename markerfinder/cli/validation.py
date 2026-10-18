"""Input validation and taxonomy/state resolution helpers for the CLI.

Verbatim relocation of the validation and input-resolution helpers that used to
live in ``markerfinder.__main__``. Every function body is preserved exactly so
that argument checks, exit codes, and warnings are unchanged.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

from markerfinder.phases import canonical_log_tag
from markerfinder.exceptions import (
    CrossValidationError,
    MultiTreeError,
    TaxonomyConflictError,
    TreeValidationError,
)
from markerfinder.cli.constants import EXIT_ARG_ERROR, EXIT_DATA_ERROR

logger = logging.getLogger(__name__)


# Pipeline "core" output artifacts. When any of these already exist in the
# Output directory, the run is treated as a collision and gated by --force /
# --no_clobber. Other files (user notes, logs, prior non-MarkerFinder
# Artifacts) are treated as foreign and ignored — MarkerFinder just coexists
# With them rather than refusing to start.
#
# Note: `.markerfinder` and `Phase5_metadata` are stepwise state / metadata
# Directories. They are listed here so a fresh aggregate run detects them,
# But subcommands (stepwise mode) explicitly allow them to already exist
# Because they are created by the preceding step.
_CORE_OUTPUT_PARTS = {
    ".markerfinder",
    "Phase5_reports",
    "Phase4_trees",
    "Phase4_alignments",
    "Phase5_metadata",
    "Phase4_intermediate",
    "tmp",
}

# Directories that are expected to exist between stepwise subcommands.
_STEPWISE_STATE_DIRS = {".markerfinder", "Phase5_metadata"}


def _detect_input_domain(
    args: argparse.Namespace,
    table_taxa: Optional[dict],
) -> Optional[str]:
    """Detect the clade domain of the input genomes: "archaea" / "bacteria".

    Resolution order:
      1. An already-resolved taxonomy map (external table or Format A labels).
      2. The Format A domain token (``_d_Archaea_`` / ``_d_Bacteria_``) embedded
         in the genome file stems / FASTA record ids.

    Returns the majority domain across genomes, or ``None`` when the genomes
    span multiple domains or no domain signal could be read — in which case the
    caller should fall back / fail loudly.
    """
    from collections import Counter
    from markerfinder.taxonomy import parse_taxonomy_embedded
    from markerfinder.utils.io import iter_genome_ids

    def _domain_from_label(label: str) -> Optional[str]:
        # Fast path for Format A labels such as ``..._d_Archaea_p_...``.
        low = label.lower()
        if "_d_archaea" in low:
            return "archaea"
        if "_d_bacteria" in low:
            return "bacteria"
        # Slow path: the full parser tolerates various delimiters.
        try:
            parsed = parse_taxonomy_embedded(label, mode=args.taxonomy_delimiter_mode or "reverse")
        except Exception:
            return None
        domain = parsed.get("domain")
        if isinstance(domain, str) and domain:
            d = domain.lower()
            if "archaea" in d:
                return "archaea"
            if "bacteria" in d:
                return "bacteria"
        return None

    labels: list = []
    if table_taxa:
        labels = list(table_taxa.keys())
    else:
        try:
            labels = list(iter_genome_ids(args.input))
        except (FileNotFoundError, ValueError):
            labels = []

    counts: Counter = Counter()
    for label in labels:
        d = _domain_from_label(label)
        if d:
            counts[d] += 1
    if not counts:
        return None
    top, n_top = counts.most_common(1)[0]
    # Require a clear majority so we don't silently pick a domain for a
    # Multi-domain / unclassifiable input.
    if n_top >= max(1, len(labels) / 2):
        return top
    return None


def _discover_bundled_hmm_dir(
    args: argparse.Namespace,
    config: object,
    table_taxa: Optional[dict],
) -> Optional[str]:
    """Resolve the bundled TIGRFAM/Pfam HMM library directory for hmm mode.

    Honors ``config.selection_config.marker_db_source``:
      - ``auto``: pick ar53/bac120 under ``--db-dir`` by detected domain.
      - ``gtdb``: only consider the gtdb_markers (ar53/bac120) tree.
      - ``cog``: only consider the cog_hmm tree.
      - ``none``: never called (guarded by the caller), but returns None here.

    Returns the absolute path to the per-marker HMM directory (containing
    one ``.HMM``/``.hmm`` per marker), or None when no bundled library can be
    located. Callers fall back to the existing "required --marker-hmm-dir" error.
    """
    from markerfinder.modules.marker_selection import gather_marker_ids

    source = (config.selection_config.marker_db_source or "auto").lower()
    if source == "none":
        logger.debug("  Auto HMM discovery disabled (--marker-db-source none).")
        return None

    db_dir = getattr(config, "db_dir", None) or args.db_dir or "./db"
    db_root = Path(db_dir)
    if not db_root.is_dir():
        logger.debug("  Auto HMM discovery: --db-dir '%s' does not exist.", db_dir)
        return None

    # Candidate (domain_subdir, root_label) pairs, ordered by preference for the
    # Requested source. ``(name, path)`` per configured marker-db source.
    candidates: list = []
    gtdb_root = db_root / "gtdb_markers"
    cog_root = db_root / "cog_hmm"

    def _gtdb_candidates(selected_domain: Optional[str]) -> list:
        pairs = []
        if selected_domain == "archaea":
            pairs.append(("ar53", gtdb_root / "ar53"))
            pairs.append(("bac120", gtdb_root / "bac120"))
        elif selected_domain == "bacteria":
            pairs.append(("bac120", gtdb_root / "bac120"))
            pairs.append(("ar53", gtdb_root / "ar53"))
        else:
            # Unknown domain: prefer whichever subdir exists, ar53 first.
            pairs.append(("ar53", gtdb_root / "ar53"))
            pairs.append(("bac120", gtdb_root / "bac120"))
        return pairs

    def _cog_candidates() -> list:
        return [("cog_hmm", cog_root)]

    if source == "cog":
        candidates = _cog_candidates()
    else:  # Auto / gtdb
        domain = _detect_input_domain(args, table_taxa) if source == "auto" else None
        candidates = _gtdb_candidates(domain)

    for name, path in candidates:
        if not path.is_dir():
            continue
        try:
            ids = gather_marker_ids(str(path))
        except FileNotFoundError:
            continue
        if ids:
            return str(path.resolve())

    return None


def _check_output_conflicts(args: argparse.Namespace) -> None:
    """Check if output files already exist and handle --force/--no-clobber.

    Collision is only raised when a *pipeline core* artifact (see
    ``_CORE_OUTPUT_PARTS``) is already present. Foreign files (e.g. a
    user-provided logfile) are logged at INFO and do not block startup.
    In stepwise mode, state/metadata directories created by the previous
    step are allowed to exist without --force.
    """
    output_dir = Path(args.output)
    is_stepwise = getattr(args, "command", None) is not None

    if not output_dir.exists():
        return

    # ``-o`` naming an existing FILE used to reach ``iterdir`` and die with a
    # Bare NotADirectoryError from inside os.listdir — an unhandled traceback
    # For what is simply a bad argument ( asks for a distinguishable,
    # Actionable failure instead).
    if output_dir.is_file():
        logger.error(
            f"Output path {args.output} is an existing FILE; MarkerFinder writes "
            "its results into a directory. Choose a different -o path or remove "
            "the file."
        )
        sys.exit(EXIT_ARG_ERROR)

    try:
        children = [p for p in output_dir.iterdir()]
    except PermissionError:
        logger.error(f"Cannot read output directory: {args.output}")
        sys.exit(EXIT_ARG_ERROR)

    if not children:
        return

    core = [p for p in children if p.name in _CORE_OUTPUT_PARTS]
    foreign = [p for p in children if p.name not in _CORE_OUTPUT_PARTS]

    if foreign and not core:
        logger.info(
            f"Output directory {args.output} contains non-pipeline files "
            f"(e.g. {[p.name for p in foreign][:5]}); will coexist with them."
        )
        return

    if core:
        # In stepwise mode, state/metadata dirs from a previous step are not
        # Considered a conflict on their own.
        if is_stepwise:
            non_state_core = [p for p in core if p.name not in _STEPWISE_STATE_DIRS]
            if not non_state_core:
                return
            core = non_state_core

        if args.no_clobber:
            logger.info(
                f"Output directory {args.output} already has pipeline outputs, "
                "will skip existing (--no-clobber)."
            )
            return
        if not args.force:
            logger.error(
                f"Output directory {args.output} already contains pipeline "
                f"outputs (e.g. {[p.name for p in core][:5]}). Use --force to "
                "overwrite or --no-clobber to skip."
            )
            sys.exit(EXIT_ARG_ERROR)


def _validate_inputs(args: argparse.Namespace) -> None:
    """Validate all input paths and arguments.

    For non-scan subcommands, -i/--input may be absent (genomes are loaded from
    persisted state); the input-directory check is skipped in that case.
    """
    if args.input is not None:
        input_path = Path(args.input)
        if not input_path.is_dir():
            logger.error(f"Input path is not a directory: {args.input}")
            sys.exit(EXIT_DATA_ERROR)

    if args.threads < 1:
        logger.error(f"Thread count must be >= 1, got {args.threads}")
        sys.exit(EXIT_ARG_ERROR)

    # --hgt-threshold is a deprecated alias and may be unset (None); only a
    # Value the user actually passed needs range validation.
    if args.hgt_threshold is not None and (
            args.hgt_threshold < 0 or args.hgt_threshold > 1):
        logger.error(f"HGT threshold must be 0-1, got {args.hgt_threshold}")
        sys.exit(EXIT_ARG_ERROR)

    # --min-occupancy 与其废弃别名 --min-marker-coverage 均可为 None(未指定,
    # 沿用 Phase 0 自适应阈值), 指定时必须落在 0-1.
    for _occ_name in ("min_occupancy", "min_marker_coverage"):
        _occ = getattr(args, _occ_name, None)
        if _occ is not None and (_occ < 0 or _occ > 1):
            logger.error(f"{_occ_name.replace('_', '-')} must be 0-1, got {_occ}")
            sys.exit(EXIT_ARG_ERROR)


def _handle_tree_input(args: argparse.Namespace) -> Optional[object]:
    """Handle reference species tree input with validation and multi-tree detection."""
    tree_path = getattr(args, "species_tree", None) or getattr(args, "tree", None)
    if not tree_path:
        return None

    from markerfinder.validation import (
        count_trees_in_file,
        detect_tree_format,
        load_tree,
    )

    fmt = detect_tree_format(tree_path)
    n_trees = count_trees_in_file(tree_path, fmt)

    if n_trees == 0:
        logger.error(f"No trees found in {tree_path} (format: {fmt})")
        sys.exit(EXIT_DATA_ERROR)

    if n_trees > 1:
        logger.warning(
            f"Multiple trees ({n_trees}) detected in {tree_path} (format: {fmt})"
        )

        mode = args.multi_tree_mode
        if mode == "ask":
            raise MultiTreeError(
                f"Multiple trees found in {tree_path}. "
                f"Re-run with --multi-tree-mode {{first,last,random}} to proceed."
            )

        # 各模式真正落实到 load_tree 的 tree_index: first=第 1 棵,
        # Last=最后 1 棵, random=确定性随机抽取 1 棵. split 需要对每棵树独立
        # 跑完整下游流程, 参考树入口不支持, 显式报错而非只打日志后照旧.
        if mode == "split":
            raise MultiTreeError(
                f"--multi-tree-mode split is not supported for a reference "
                f"species tree ({tree_path}); use first/last/random instead."
            )
        elif mode == "first":
            tree_index = 0
            logger.info(f"Using first tree (#1/{n_trees}) from {tree_path}")
        elif mode == "last":
            tree_index = n_trees - 1
            logger.info(f"Using last tree (#{n_trees}/{n_trees}) from {tree_path}")
        elif mode == "random":
            import random as _random
            selected = _random.randint(1, n_trees)
            tree_index = selected - 1
            logger.info(f"Random mode: selected tree #{selected}/{n_trees} from {tree_path}")
        else:
            tree_index = None
    else:
        tree_index = None

    try:
        tree = load_tree(tree_path, validate=True, strip_annotations=args.strip_annotations,
                         tree_index=tree_index)
    except TreeValidationError:
        raise
    except FileNotFoundError:
        raise

    if tree.n_tips > 10000:
        logger.info(
            f"Large tree detected ({tree.n_tips} tips). "
            f"Consider reducing thread count for memory."
        )

    return tree


def _handle_cross_validation(args: argparse.Namespace) -> None:
    """Run tree-sequence cross-validation if both inputs provided."""
    tree_path = getattr(args, "species_tree", None) or getattr(args, "tree", None)
    if not tree_path or not args.sequences:
        return

    if args.no_cross_check:
        logger.info("Cross-validation skipped (--no-cross-check)")
        return

    from markerfinder.validation import cross_validate

    try:
        report = cross_validate(
            tree_path, args.sequences, strict=True,
            # ``--mol-type`` and ``--skip-length-check`` are assertions over the
            # Sequence input; this is the only place the sequence file is read,
            # So dropping them here is what made the options inert.
            mol_type=getattr(args, "mol_type", None),
            skip_length_check=bool(getattr(args, "skip_length_check", False)),
        )
        logger.info(f"Cross-validation passed: {len(report.tree_tips)} tips matched")
    except CrossValidationError as e:
        logger.error(str(e))
        sys.exit(EXIT_DATA_ERROR)


def _load_taxonomy(args: argparse.Namespace) -> Optional[dict]:
    """Load external taxonomy table if provided."""
    if not args.taxonomy_table:
        return None

    from markerfinder.taxonomy import load_taxonomy_table, validate_table_tree_consistency

    try:
        table_taxa = load_taxonomy_table(
            args.taxonomy_table,
            sep=args.table_sep if args.table_sep else None,
            taxonomy_format=args.taxonomy_format,
            delimiter_mode=args.taxonomy_delimiter_mode,
            ignore_malformed=args.ignore_malformed,
            custom_levels=getattr(args, "taxonomy_levels", None),
        )
    except (FileNotFoundError, ValueError, TaxonomyConflictError) as e:
        # TaxonomyConflictError(控制字符/循环引用等用户输入错误)是 PhyloToolError
        # 而非 ValueError 子类, 必须显式捕获, 否则会以裸 traceback 逃逸 CLI 错误
        # 处理而非走规范的 EXIT_DATA_ERROR 退出码.
        logger.error(f"Failed to load taxonomy table: {e}")
        sys.exit(EXIT_DATA_ERROR)

    return table_taxa


def _auto_parse_embedded_taxonomy(args: argparse.Namespace, force: bool) -> Optional[dict]:
    """Parse Format A taxonomy embedded in genome identifiers/filenames.

    In ``--marker-mode hmm`` the built-in TIGRFAM/Pfam marker set is scanned
    directly per genome, so the only way the phylogenetic HGT step can run is
    the MAD + monophyly path — which needs a taxonomy map. When none is supplied
    we can recover it from Format A labels such as
    ``GB_GCA_..._d_Archaea_p_Halobacteriota_c_...`` carried either by the genome
    file stem or the FASTA record ids.

    Args:
        args: parsed CLI arguments (needs ``args.input`` and ``args.auto_embed_taxonomy``).
        force: when True (user passed ``--taxonomy-format embedded``), parse every
            genome regardless of how many look like Format A. when False
            (auto-detect), only enable when a clear majority (>= 50%%) of the
            genomes carry a resolvable ``domain`` rank, so we do not impose bogus
            taxonomy structure on arbitrary identifiers.

    Returns:
        A ``{genome_id: {rank: value}}`` taxonomy map, or ``None`` if nothing
        resolvable was found (or the auto-detect heuristic did not trigger).
    """
    try:
        from markerfinder.utils.io import iter_genome_ids
        from markerfinder.taxonomy import parse_taxonomy_embedded
    except Exception as e:  # Pragma: no cover - best-effort helper
        logger.debug("  Embedded-taxonomy helper import failed: %s", e)
        return None

    delimiter_mode = args.taxonomy_delimiter_mode or "reverse"
    try:
        genome_ids = list(iter_genome_ids(args.input))
    except (FileNotFoundError, ValueError) as e:
        logger.debug("  Cannot iterate genome ids for embedded-taxonomy detection: %s", e)
        return None
    if not genome_ids:
        return None

    parsed_map: dict = {}
    n_with_domain = 0
    for gid in genome_ids:
        parsed = parse_taxonomy_embedded(gid, mode=delimiter_mode)
        if parsed.get("domain") is not None:
            n_with_domain += 1
        if force or any(v is not None for v in parsed.values()):
            # Keep entries that resolve at least one rank; for the auto-detect
            # (non-forced) path we still collect them and decide below by majority.
            if any(v is not None for v in parsed.values()):
                parsed_map[gid] = parsed

    if not parsed_map:
        return None

    # Auto-detect gate: require a clear majority of genomes to expose a domain.
    # This avoids embedding taxonomy on identifiers that merely happen to contain
    # A few underscore-separated segments without real rank semantics.
    if not force:
        majority = n_with_domain >= max(1, len(genome_ids) / 2)
        if not majority:
            logger.debug(
                "  No embedded taxonomy auto-detected: only %d/%d genome ids expose a "
                "Format A 'domain' rank (need >= 50%%). Pass --taxonomy-format embedded to "
                "override, or --taxonomy-table to supply a mapping.",
                n_with_domain, len(genome_ids),
            )
            return None

    logger.info(
        f"{canonical_log_tag('0.2')} Parsed embedded taxonomy for %d tips from genome ids "
        "(format=%s, mode=%s, auto-detected=%s).",
        len(parsed_map),
        "embedded",
        delimiter_mode,
        "yes" if not force else "no (forced)",
    )
    return parsed_map
