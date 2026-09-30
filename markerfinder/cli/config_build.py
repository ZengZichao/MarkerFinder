"""Build a:class:`PipelineConfig` from parsed CLI arguments.

Verbatim relocation of ``_build_pipeline_config`` from ``markerfinder.__main__``.
The only structural change is that the previously *local* ``mode_strategy`` dict
is now a module-level constant (``_MODE_STRATEGY``); the resolved mapping and the
``SelectionStrategy`` chosen for a given ``--mode`` are identical.
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

from markerfinder.config import (
    AlignerConfig,
    HGTConfig,
    MAGConfig,
    OrthologConfig,
    PhylogeneticConfig,
    PipelineConfig,
    ReportConfig,
    SelectionConfig,
    TaxonomyConfig,
)
from markerfinder.config import HGT_LEVEL1_MAX, HGT_LEVEL2_MAX
from markerfinder.config import _default_tmp_dir as _config_tmp_root
from markerfinder.models.marker import SelectionStrategy
from markerfinder.cli.constants import EXIT_DATA_ERROR

logger = logging.getLogger(__name__)

# Maps the CLI ``--mode`` value onto a marker-selection strategy. Previously a
# Local dict inside ``_build_pipeline_config``; hoisted to module level so it is
# Defined once (behaviourally identical).
_MODE_STRATEGY = {
    "conservative": SelectionStrategy.INFO_MAX,
    "standard": SelectionStrategy.GREEDY,
    "expanded": SelectionStrategy.RATE_BALANCED,
    "mag_adaptive": SelectionStrategy.SPARSE_OPTIMIZED,
}


def _build_pipeline_config(args: argparse.Namespace, parser: Optional[argparse.ArgumentParser] = None) -> PipelineConfig:
    """Build a PipelineConfig from parsed CLI args.

    Shared by the full 'run' aggregate and the per-step subcommands. For non-scan
    steps (where -i/--input is typically absent), input_dir falls back to the
    output directory so config construction never fails.
    """
    # NOTE: keep ``file_prefix`` as a bare filename (no directory component).
    # It is used as the basename of tmp files (e.g. via
    # ``os.path.join(tmp_dir, f"{prefix}.concat.fasta")``); embedding the
    # Output directory here would nest tmp files under an
    # ``tmp_dir/output_basename/`` subpath that is never created (ENOENT).
    file_prefix = "markerfinder"
    # Non-scan subcommands do not require -i; fall back to output dir so the
    # Config object is always constructible (genomes are loaded from state).
    input_dir = args.input or args.output

    # MarkerFinder requires protein input; the composition HGT step has been
    # Removed. 接受.faa/.fasta/.fa —— 与 load_genomes_from_directory 的
    # 实际加载能力及 README「输入格式」一致, 之前只检查 *.faa 会把纯.fasta 输入
    # 误报为"核苷酸输入不支持".
    if args.input:
        from pathlib import Path as _Path
        input_path = _Path(args.input)
        if not any(input_path.glob(ext) for ext in ("*.faa", "*.fasta", "*.fa")):
            logger.error(
                "MarkerFinder requires protein input (.faa/.fasta/.fa files). "
                "The composition HGT step has been removed; nucleotide-only input is not supported."
            )
            sys.exit(EXIT_DATA_ERROR)

    # Validate the HGT step names. 'phylogenetic' drives the risk grading;
    # 'composition' turns on the parallel composition/GC
    # Evidence columns only — it never enters overall_risk/overall_score.
    _HGT_STEP_NAMES = ("phylogenetic", "composition")
    hgt_steps_raw = [s.strip().lower() for s in args.hgt_steps.split(",") if s.strip()]
    unknown_steps = [s for s in hgt_steps_raw if s not in _HGT_STEP_NAMES]
    if unknown_steps:
        msg = (
            f"unrecognized --hgt-step(s): {unknown_steps!r}. "
            f"Valid choices are: {list(_HGT_STEP_NAMES)}"
        )
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)
    if not hgt_steps_raw:
        msg = "--hgt-steps must contain: phylogenetic"
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)
    enable_phylogenetic = "phylogenetic" in hgt_steps_raw
    composition_screen_requested = "composition" in hgt_steps_raw

    # ``--hgt-threshold`` is a deprecated alias of
    # ``--hgt-threshold-l1l2``. Unset is None (that is how we know the user did
    # Not pass it); passing it keeps the historical behaviour — it moves ONLY
    # The Level1/Level2 band — and prints exactly one deprecation notice, which
    # Is what asks for and what used to be merely a comment claiming it
    # Happened.
    hgt_threshold_legacy = getattr(args, "hgt_threshold", None)
    if hgt_threshold_legacy is not None:
        logger.warning(
            "--hgt-threshold is deprecated: it is only an alias of "
            "--hgt-threshold-l1l2 and moves the Level1/Level2 band alone "
            "(Level2/Level3 stays under --hgt-threshold-l2l3). Use "
            "--hgt-threshold-l1l2 / --hgt-threshold-l2l3."
        )
    hgt_threshold_l1l2 = (
        args.hgt_threshold_l1l2
        if getattr(args, "hgt_threshold_l1l2", None) is not None
        else (
            hgt_threshold_legacy
            if hgt_threshold_legacy is not None else HGT_LEVEL1_MAX
        )
    )

    # A per-run scratch directory when the user did not name one; the function is
    # Called per build so two pipelines in one process do not share it.
    tmp_dir = args.tmp_dir if args.tmp_dir else _config_tmp_root()
    # 未显式指定 --tmp-dir 时为自动生成的每趟独立目录, 流程结束可安全清理;
    # 显式指定时不自动清理(与 --tmp-dir 帮助文本一致).
    tmp_dir_auto = not bool(args.tmp_dir)

    # 用户显式指定的占用率下限: --min-occupancy 优先, 废弃别名
    # --min-marker-coverage 兜底. 两者都未提供时保持原行为(Phase 0 自适应阈值).
    user_min_occupancy = getattr(args, "min_occupancy", None)
    if user_min_occupancy is None:
        user_min_occupancy = getattr(args, "min_marker_coverage", None)
        if user_min_occupancy is not None:
            logger.warning(
                "  --min-marker-coverage is a deprecated alias for "
                "--min-occupancy; the value is applied to marker-selection "
                "occupancy. Use --min-occupancy in new scripts."
            )

    # --max-markers is CLI-settable. None means "not given", in
    # Which case the Phase-0 adaptive budget applies; an explicit value must
    # Survive Phase 0 (see SelectionConfig.user_max_markers).
    user_max_markers = getattr(args, "max_markers", None)

    # Every other documented deprecated alias announces itself here. The help
    # Text and the MANUAL already label them as deprecated; a run that accepted
    # Them in silence left a user following an old pipeline script with no way
    # To learn that the flag they typed is a compatibility shim (or, for
    # --no-hgt-adaptive-thresholds, a documented no-op).
    if getattr(args, "used_tree_alias", False):
        logger.warning(
            "  --tree is a deprecated alias for --species-tree; use "
            "--species-tree in new scripts."
        )
    if getattr(args, "fast_tree", False):
        logger.warning(
            "  --fast-tree is a deprecated alias for --gene-tree-builder "
            "fasttree (which is also the default); use --gene-tree-builder."
        )
    if getattr(args, "no_hgt_adaptive_thresholds", False):
        logger.warning(
            "  --no-hgt-adaptive-thresholds restates the default (adaptive "
            "far-distance relaxation is off); it changes nothing."
        )

    # Resolve gene-tree builder selection. Explicit --gene-tree-builder takes
    # Precedence over the legacy --fast-tree flag and the config file values
    # (fast_tree / gene_tree_builder). This makes the two-option choice clear
    # For coalescent inference while preserving backward compatibility.
    gene_tree_builder = args.gene_tree_builder
    if gene_tree_builder is None:
        # Default to FastTree+WAG for speed; IQ-TREE3 remains available via --gene-tree-builder iqtree.
        gene_tree_builder = "fasttree"

    pipeline_config = PipelineConfig(
        input_dir=input_dir,
        output_dir=args.output,
        output_prefix=file_prefix,
        cpus=args.threads,
        mode=args.mode,
        tmp_dir=tmp_dir,
        tmp_dir_auto=tmp_dir_auto,
        keep_tmp=args.keep_tmp,
        save_intermediates=args.save_intermediates,
        db_dir=args.db_dir,
        force=args.force,
        no_clobber=args.no_clobber,
        sequences_path=args.sequences,
        mol_type=args.mol_type,
        skip_length_check=args.skip_length_check,
        strip_annotations=args.strip_annotations,
        mag_config=MAGConfig(
            checkm_results=args.checkm_results,
            skip_checkm=args.skip_checkm,
            cpus=args.threads,
            tmp_dir=tmp_dir,
            input_dir=input_dir,
            min_marker_coverage=user_min_occupancy if user_min_occupancy is not None else 0.3,
        ),
        selection_config=SelectionConfig(
            strategy=_MODE_STRATEGY.get(args.mode, SelectionStrategy.GREEDY),
            marker_mode=args.marker_mode,
            marker_hmm_dir=args.marker_hmm_dir or "",
            marker_db_source=args.marker_db_source or "auto",
            min_hmm_score=args.min_hmm_score,
            # 显式用户输入优先于自动推断: base 默认与 user 覆盖分别记录.
            min_occupancy=user_min_occupancy if user_min_occupancy is not None else 0.75,
            user_min_occupancy=user_min_occupancy,
            max_markers=user_max_markers if user_max_markers is not None else 60,
            user_max_markers=user_max_markers,
            cpus=args.threads,
            tmp_dir=tmp_dir,
            gtdb_markers_dir=args.gtdb_markers_dir or "",
            species_tree=args.species_tree or "",
        ),
        hgt_config=HGTConfig(
            enable_phylogenetic=enable_phylogenetic,
            cpus=args.threads,
            tmp_dir=tmp_dir,
            # Two independently settable bands; the legacy
            # --hgt-threshold aliases level1_max only (one deprecation notice).
            level_thresholds={
                "level1_max": (
                    args.hgt_threshold_l1l2
                    if getattr(args, "hgt_threshold_l1l2", None) is not None
                    else hgt_threshold_l1l2
                ),
                "level2_max": (
                    args.hgt_threshold_l2l3
                    if getattr(args, "hgt_threshold_l2l3", None) is not None
                    else HGT_LEVEL2_MAX
                ),
            },
            # Criterion mode fields (consumed by the gate).
            hgt_mode=getattr(args, "hgt_mode", "risk"),
            consistency_stringency=getattr(args, "consistency_stringency", 1),
            min_informative_sites=getattr(args, "min_informative_sites", 0),
            monophyly_rank=args.monophyly_rank,
            monophyly_threshold=args.monophyly_threshold,
            # Both CLI flags are opt-in store_true (default False), so the
            # Config default (adaptive_far_thresholds=False) is preserved unless
            # The user explicitly passes --hgt-adaptive-thresholds.
            adaptive_far_thresholds=bool(
                getattr(args, "hgt_adaptive_thresholds", False)
                and not getattr(args, "no_hgt_adaptive_thresholds", False)
            ),
        ),
        ortholog_config=OrthologConfig(
            cpus=args.threads,
            tmp_dir=tmp_dir,
        ),
        aligner_config=AlignerConfig(
            cpus=args.threads,
            tmp_dir=tmp_dir,
            output_prefix=file_prefix,
        ),
        phylo_config=PhylogeneticConfig(
            cpus=args.threads,
            tmp_dir=tmp_dir,
            output_prefix=file_prefix,
            ufboot_replicates=args.ufboot,
            use_fasttree=(gene_tree_builder == "fasttree") if gene_tree_builder else args.fast_tree,
            gene_tree_builder=gene_tree_builder,
            coalescent_mode=args.coalescent_mode,
            # Optional support gate. ``None`` (default, when --min-gene-tree-support
            # Is not passed) keeps the historical n_tips >= 4-only filter.
            min_gene_tree_support=getattr(args, "min_gene_tree_support", None),
        ),
        report_config=ReportConfig(
            output_dir=args.output,
            output_prefix=file_prefix,
            report_format=args.report_format,
            save_intermediates=args.save_intermediates,
            # 报告风险分桶与实际 HGT 分级同源, 避免图示边界与分级漂移.
            level1_max=(
                args.hgt_threshold_l1l2
                if getattr(args, "hgt_threshold_l1l2", None) is not None
                else hgt_threshold_l1l2
            ),
            level2_max=(
                args.hgt_threshold_l2l3
                if getattr(args, "hgt_threshold_l2l3", None) is not None
                else HGT_LEVEL2_MAX
            ),
            # Coverage warning floor.
            require_evidence_coverage=getattr(args, "require_evidence_coverage", 0.5),
            # Threshold scan switch.
            hgt_scan=bool(getattr(args, "hgt_scan", False)),
            # The scan stability floor must actually reach the
            # Scanner, otherwise --scan-stability-min is a dead flag.
            scan_stability_min=float(
                getattr(args, "scan_stability_min", 0.6) or 0.6
            ),
            # Composition evidence is opt-in via
            # --hgt-steps phylogenetic,composition. Parallel columns only.
            composition_screen=composition_screen_requested,
            # Reportable statement of what actually adjudicated.
            hgt_mode=str(getattr(args, "hgt_mode", "risk") or "risk"),
            # Optional COG category map (None => NA column + note).
            cog_category_map=getattr(args, "cog_category_map", None),
            taxonomy_mustpass=getattr(args, "taxonomy_mustpass", None),
            # Assertion switches wired from the CLI.
            strict_assertions=bool(getattr(args, "strict_assertions", True)),
            allow_assertion_failure=bool(getattr(args, "allow_assertion_failure", False)),
        ),
        taxonomy_config=TaxonomyConfig(
            taxonomy_table=args.taxonomy_table,
            taxonomy_format=args.taxonomy_format,
            taxonomy_source_priority=args.taxonomy_source_priority,
            taxonomy_delimiter_mode=args.taxonomy_delimiter_mode,
            table_sep=args.table_sep,
            ignore_malformed=args.ignore_malformed,
            taxonomy_levels=args.taxonomy_levels,
        ),
    )

    # An explicit --marker-preset overrides the
    # Selection triple (min_occupancy / max_markers / strategy); 'none'
    # (default) keeps the configured behaviour untouched.
    preset_name = getattr(args, "marker_preset", "none")
    if preset_name and preset_name != "none":
        from markerfinder.models.marker import (
            RESOLUTION_PRESETS,
            ResolutionPreset,
        )

        preset = RESOLUTION_PRESETS.get(ResolutionPreset(preset_name))
        if preset is not None:
            pipeline_config.selection_config.min_occupancy = preset.min_occupancy
            pipeline_config.selection_config.max_markers = preset.max_markers
            # An explicit preset is an explicit instruction: it must also win
            # Over the Phase-0 adaptive budget, exactly like --max-markers.
            pipeline_config.selection_config.user_max_markers = preset.max_markers
            pipeline_config.selection_config.selection_strategy = preset.selection_strategy
            logger.info(
                f"  Applied --marker-preset {preset_name}: "
                f"min_occupancy={preset.min_occupancy}, "
                f"max_markers={preset.max_markers}, "
                f"strategy={preset.selection_strategy.value}"
            )

    return pipeline_config
