from __future__ import annotations

import json
import logging
import os
import shutil
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set

from markerfinder.phases import canonical_log_tag
from markerfinder.config import PipelineConfig
from markerfinder.models.genome import Genome
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.tree import Tree
from markerfinder.models.pipeline_types import (
    HGTReport,
    MarkerSelectionResult,
    PhylogeneticResult,
    PreprocessingResult,
)
from markerfinder.models.report import (
    PhaseContext,
    PipelineResult,
    PlainTextReportOutput,
    ReportOutput,
    RuntimeInfo,
)
from markerfinder.modules.hgt_filter import HGTFilterModule
from markerfinder.modules.mag_optimization import MAGOptimizationModule
from markerfinder.modules.marker_selection import AdaptiveMarkerSelectionModule
from markerfinder.modules.phylogenetic_inference import PhylogeneticInferenceModule
from markerfinder.modules.report_generator import ReportGeneratorModule
from markerfinder.utils.db_versioning import DatabaseVersionManager
from markerfinder.utils.state_codec import (
    STATE_FILE_NAME,
    StateSchemaError,
    load_pipeline_state,
    prev_step_product_rel,
    save_pipeline_state,
    state_dir_for,
)
from markerfinder.utils.tree_source import has_usable_tree, prioritize_tree_source
from markerfinder.utils.dependency_check import precheck_all
from markerfinder._version import get_version_string

logger = logging.getLogger(__name__)


class MarkerFinderPipeline:
    def __init__(self, config: PipelineConfig):
        self.config = config
        # Run dependency pre-checks at pipeline startup to detect issues early.
        precheck_all(strict=False)
        self.mag_optimization = MAGOptimizationModule(config.mag_config)
        self.marker_selector = AdaptiveMarkerSelectionModule(config.selection_config)
        self.hgt_filter = HGTFilterModule(
            config.hgt_config,
            tmp_dir=config.tmp_dir,
            gene_trees_dir=str(Path(config.output_dir) / "Phase4_trees" / "gene_trees"),
            cpus=config.cpus,
        )
        self.phylo_inference = PhylogeneticInferenceModule(
            config.phylo_config,
            gene_trees_dir=str(Path(config.output_dir) / "Phase4_trees" / "gene_trees"),
        )
        self.report_generator = ReportGeneratorModule(config.report_config, force=config.force, no_clobber=config.no_clobber)
        if config.ortholog_config.tmp_dir == "/tmp" and config.tmp_dir != "/tmp":
            config.ortholog_config.tmp_dir = config.tmp_dir

    STEP_ORDER = ["scan", "filter", "infer", "report"]

    def run(
        self,
        genomes: List[Genome],
        reference_tree: Optional[object] = None,
        table_taxa: Optional[dict] = None,
        *,
        taxonomy_origin: str = "external table",
    ) -> PipelineResult:
        """Full pipeline run (aggregates all 5 steps in-memory).

        Args:
            genomes: input genomes.
            reference_tree: optional user-provided tree.
            table_taxa: optional taxonomy map (from an external table or
                auto-parsed Format A labels).
            taxonomy_origin: human-readable origin of ``table_taxa`` for the
                Phase 0.2 log line.
        """
        start_time = time.time()
        # The run-level 'what was NOT measured' ledger is per run.
        from markerfinder.utils import unmeasured

        unmeasured.clear()

        if len(genomes) > 10000:
            logger.info(
                f"Large dataset detected ({len(genomes)} genomes). "
                f"Consider using fewer threads to reduce memory usage."
            )

        context = self._prepare_context(reference_tree, table_taxa, taxonomy_origin)

        # Step 1: scan (Phase 0 + 1 + 1.5)
        scan_data = self.run_scan(genomes, context)

        # Step 2: filter (Phase 2)
        filter_data = self.run_filter(genomes, scan_data, context, taxonomy_map=table_taxa)

        # Step 3: infer (Phase 3)
        infer_data = self.run_infer(genomes, scan_data, filter_data, context)

        # Step 4: report (Phase 4; output dirs keep the historical Phase5_* names)
        end_time = time.time()
        runtime_info = RuntimeInfo(start_time=start_time, end_time=end_time, duration=end_time - start_time)
        report_data = self.run_report(scan_data, filter_data, infer_data, context, runtime_info)

        self._write_run_config(start_time, end_time)
        # The log has to carry the same declaration the report does —
        # A reader who only kept the terminal output must still learn which
        # Metrics are missing and for what reason.
        unmeasured.declare_in_log(logger)
        logger.info(f"\nPipeline completed in {runtime_info.duration:.1f}s")

        # 在清理 tmp 前,把有用的中间产物归档到 output/Phase4_intermediate/(若用户开启).
        if self.config.save_intermediates:
            self._save_intermediates()

        # 清理 per-run 临时目录(自动生成且未 --keep-tmp 时). 见 _cleanup_tmp.
        self._cleanup_tmp()

        # ── 可用性信号: 在 __main__.py 汇总检验中判定「管线是否真正成功」 ──
        # 可用标记 = HGT 过滤后剩下(等级 1/2)可用于建树的标记数.
        # 物种树来源优先级判定收敛到 markerfinder.utils.tree_source
        # (prioritize_tree_source), 与报告逻辑共用同一份实现, 避免互相依赖与漂移.
        phylo_result = infer_data.get("phylo_result")
        coalescent = getattr(phylo_result, "coalescent", None)
        coalescent_source = getattr(coalescent, "species_tree_source", None)
        has_coalescent_tree = bool(
            coalescent_source and coalescent_source != "none"
            and has_usable_tree(getattr(coalescent, "species_tree", None))
        )
        has_concat_tree = has_usable_tree(
            getattr(getattr(phylo_result, "supermatrix", None), "tree", None)
        )
        resolved_source = prioritize_tree_source(
            coalescent_source if has_coalescent_tree else None,
            has_concat_tree,
            concat_label="concat",
        )

        return PipelineResult(
            preprocessing=scan_data["preprocessing"],
            marker_selection=scan_data["marker_selection"],
            hgt_report=filter_data["hgt_report"],
            phylogenetic=infer_data["phylo_result"],
            report=report_data["html_output"],
            text_report=report_data["text_output"],
            runtime=runtime_info,
            parameters=self.config.to_dict(),
            usable_marker_count=len(filter_data.get("marker_genes") or {}),
            species_tree_source=resolved_source,
        )

    def _prepare_context(
        self,
        reference_tree: Optional[object],
        table_taxa: Optional[dict],
        taxonomy_origin: str = "external table",
    ) -> PhaseContext:
        """Build the run context and emit the Phase 0.1/0.2 intake log lines.

        Shared by ``run`` and ``run_step``. Both entry points used to carry
        their own copy of this preamble, so any change to reference-tree or
        taxonomy intake had to be made in two places — the shape of the
        drift this file has been checked for.
        """
        Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)
        Path(self.config.tmp_dir).mkdir(parents=True, exist_ok=True)

        context = PhaseContext()
        if reference_tree is not None:
            logger.info(
                f"{canonical_log_tag('0.1')} Using user-provided reference "
                f"species tree ({reference_tree.n_tips} tips)"
            )
            context.species_tree = reference_tree
        if table_taxa:
            logger.info(
                f"{canonical_log_tag('0.2')} Loaded taxonomy for {len(table_taxa)} tips "
                f"from {taxonomy_origin}"
            )
        return context

    # ── Step 1: scan (Phase 0 quality + Phase 1 marker selection + Phase 1.5 extraction) ──

    def run_scan(
        self,
        genomes: List[Genome],
        context: PhaseContext,
    ) -> dict:
        """Phase 0+1+1.5: quality preprocessing, GTDB-TK marker loading, sequence extraction.

        Returns dict with keys: preprocessing, marker_selection, marker_sequences, rank_map.
        """
        logger.info(f"{canonical_log_tag('0')} Quality preprocessing...")
        preprocessing = self.mag_optimization.run(genomes)
        context.quality_data = preprocessing.quality_results
        context.adaptive_params = preprocessing.adaptive_params

        # 显式指定的 --min-occupancy(/--min-marker-coverage)优先于 Phase 0
        # 自适应分层阈值 —— 用户显式输入必须压过自动推断, 否则该参数形同虚设.
        user_min_occ = getattr(self.config.selection_config, "user_min_occupancy", None)
        if user_min_occ is not None and context.adaptive_params is not None:
            context.adaptive_params.min_occupancy = float(user_min_occ)
            logger.info(f"  {canonical_log_tag('0')} User-specified min occupancy override: {user_min_occ}")

        # Same rule for the marker budget (--max-markers / explicit preset).
        # Without this, the recorded selection_config.max_markers was honored
        # Only when Phase 0 happened not to produce adaptive params, i.e. in a
        # Normal run the flag was read, reported, and then discarded.
        user_max_markers = getattr(
            self.config.selection_config, "user_max_markers", None
        )
        if user_max_markers is not None and context.adaptive_params is not None:
            context.adaptive_params.max_markers = int(user_max_markers)
            logger.info(
                f"  {canonical_log_tag('0')} User-specified marker budget "
                f"override: {user_max_markers}"
            )

        if self.config.selection_config.marker_mode == "hmm":
            logger.info(f"{canonical_log_tag('1')} HMM scan against external TIGRFAM/Pfam HMM profiles...")
            if not self.config.selection_config.marker_hmm_dir:
                raise ValueError("--marker-hmm-dir is required when --marker-mode hmm")
            marker_selection = self.marker_selector.run(genomes, context)
            logger.info(f"{canonical_log_tag('1.5')} Extracting marker protein sequences from genome FASTA files...")
            marker_sequences = self.marker_selector.extract_marker_sequences(
                genomes,
                marker_selection.marker_set,
                marker_selection.occupancy_matrix,
                tmp_dir=self.config.tmp_dir,
            )
            gene_trees = None
        else:
            gtdb_dir = self.config.selection_config.gtdb_markers_dir
            if not gtdb_dir:
                raise ValueError("--gtdb-markers-dir is required when --marker-mode gtdb_tk")
            logger.info(f"{canonical_log_tag('1')} GTDB-TK marker loading (ar53/bac120)...")
            # 基因树直接写入 <output>/Phase4_trees/gene_trees 缓存: infer 的
            # 合并法与信息量回填都只认这个目录; 之前写进 tmp/gene_trees 会让 infer
            # 全量重复建树, 且被 HGT 排除标记的树在 tmp 清理后彻底丢失.
            scan_gene_trees_dir = str(
                Path(self.config.output_dir) / "Phase4_trees" / "gene_trees"
            )
            marker_selection, marker_sequences, gene_trees = self.marker_selector.run_gtdb_tk(
                genomes,
                markers_dir=gtdb_dir,
                tmp_dir=self.config.tmp_dir,
                adaptive_params=context.adaptive_params,
                gene_trees_dir=scan_gene_trees_dir,
            )

        # --species-tree 作为 HGT 系统发育步骤的参考物种树.
        # 当用户提供物种树时,Phase 2 将基因树与该物种树做 RF/quartet 比较;
        # 未提供时,Phase 2 退化为 MAD 定根 + 单系比例(需 taxonomy table).
        species_tree_path = self.config.selection_config.species_tree
        if species_tree_path:
            if not os.path.exists(species_tree_path):
                raise FileNotFoundError(f"--species-tree not found: {species_tree_path}")
            context.species_tree = Tree.read(species_tree_path)
            logger.info(
                f"{canonical_log_tag('1')} Reference species tree loaded for HGT screening "
                f"({context.species_tree.n_tips} tips)"
            )

        context.marker_set = marker_selection.marker_set
        logger.info(
            f"  Selected {len(marker_selection.marker_set.markers)} markers "
            f"(mean occupancy: {marker_selection.marker_set.mean_occupancy:.2f})"
        )
        marker_sequences = self._resolve_orthologs(genomes, marker_selection, marker_sequences)

        occupancy_ranking = sorted(
            marker_selection.marker_set.occupancy_scores.items(),
            key=lambda x: x[1], reverse=True,
        )
        rank_map = {cog: rank + 1 for rank, (cog, _) in enumerate(occupancy_ranking)}

        return {
            "preprocessing": preprocessing,
            "marker_selection": marker_selection,
            "marker_sequences": marker_sequences,
            "rank_map": rank_map,
            "gene_trees": gene_trees,
        }

    # ── Step 2: filter (Phase 2) ──

    def run_filter(
        self,
        genomes: List[Genome],
        scan_data: dict,
        context: PhaseContext,
        taxonomy_map: Optional[dict] = None,
    ) -> dict:
        """Phase 2: HGT-aware marker filtering.

        Returns dict with keys: hgt_report, precomputed_levels, marker_genes.
        """
        marker_selection = scan_data["marker_selection"]
        marker_sequences = scan_data["marker_sequences"]
        rank_map = scan_data["rank_map"]

        logger.info(f"{canonical_log_tag('2')} HGT-aware marker filtering...")
        # The structural order gate. consistency/hybrid
        # Refuse to run until every prerequisite artifact exists ( keeps
        # 'risk' the untouchable default). ``hybrid`` needs no separate
        # Decision any more: the specification fixes its semantics at the
        # Phase 2→3 boundary (``hybrid``（两者都要过）— a marker must
        # Pass the risk level screen AND grade CONSISTENT), applied in the
        # Screen block of run_infer alongside the ``excluded`` set.
        if getattr(self.config.hgt_config, "hgt_mode", "risk") != "risk":
            from markerfinder.modules.consistency_screen import guard

            guard(self.config.hgt_config.hgt_mode, self.config.hgt_config)

        hgt_report, precomputed_levels, filter_far_active = self.hgt_filter.run(
            marker_selection.marker_set.markers,
            genomes,
            marker_sequences=marker_sequences,
            occupancy_ranking=rank_map,
            prebuilt_gene_trees=scan_data.get("gene_trees"),
            taxonomy_map=taxonomy_map,
            species_tree=context.species_tree,
            monophyly_rank=self.config.hgt_config.monophyly_rank,
            monophyly_threshold=self.config.hgt_config.monophyly_threshold,
        )
        # Back-fill the real HGT risk into the per-marker quality scores. These were
        # Initialised as placeholders in Phase 1.5 (HGT phase had not run yet); now
        # That the HGT filter has evaluated each marker we can store its true
        # ``overall_risk`` so downstream quality/level reporting reflects reality.
        # UNKNOWN markers were NOT screened — their 0.0 sentinel
        # Must NOT pollute the hgt_cleanliness quality component, so their
        # Hgt_risk_score stays None ("not measured").
        for ev in getattr(hgt_report, "marker_evaluations", []) or []:
            qs = marker_selection.quality_scores.get(ev.marker_id)
            if qs is not None:
                qs.hgt_risk_score = (
                    float(ev.overall_risk) if ev.level is not MarkerLevel.UNKNOWN else None
                )
        context.hgt_evaluations = hgt_report.marker_evaluations
        level_counts = Counter(e.level for e in hgt_report.marker_evaluations)
        logger.info(
            f"  Level 1: {level_counts.get(MarkerLevel.LEVEL_1, 0)}, "
            f"Level 2: {level_counts.get(MarkerLevel.LEVEL_2, 0)}, "
            f"Level 3: {level_counts.get(MarkerLevel.LEVEL_3, 0)}"
        )

        # Taxonomy must-pass gate, on the gene trees this
        # Phase already built. Opt-in via --taxonomy-mustpass; when the flag is
        # Set but the baseline cannot be read, the run REFUSES rather than
        # Reporting a pass (an unloadable gate is not a passed gate).
        mustpass_path = getattr(
            self.report_generator.config, "taxonomy_mustpass", None
        )
        if mustpass_path:
            from markerfinder.exceptions import AssertionFailureError
            from markerfinder.modules.taxonomy_mustpass import (
                enforce_mustpass,
                evaluate_and_write,
                load_mustpass,
                mustpass_marker_scope,
            )

            def _mustpass_newick(obj):
                if obj is None:
                    return None
                if isinstance(obj, str):
                    value = obj.strip()
                    # ``scan_data[\"gene_trees\"]`` carries the *cache path* of each
                    # Per-marker tree, not its text. Handing the path to the
                    # Splitter parser reads it as a one-tip tree, so the gate
                    # Tested nothing and still reported "all required relations
                    # Hold" — a green light produced by a broken input, which is
                    # Exactly what exists to prevent. Resolve a path (or
                    # An inline newick without a.nwk-looking suffix) to content.
                    maybe = Path(value)
                    if value and not value.startswith("(") and maybe.exists():
                        try:
                            return maybe.read_text(encoding="utf-8").strip() or None
                        except OSError as exc:
                            logger.warning(
                                f"  {canonical_log_tag('2')} must-pass could not "
                                f"read gene tree {value}: {exc}"
                            )
                            return None
                    return value or None
                value = getattr(obj, "newick", None)
                if not value and hasattr(obj, "to_newick"):
                    value = obj.to_newick()
                return value or None

            mustpass_trees = {}
            for _mid, _tree in (scan_data.get("gene_trees") or {}).items():
                _nw = _mustpass_newick(_tree)
                if _nw:
                    mustpass_trees[_mid] = _nw
            mustpass_rules = load_mustpass(mustpass_path)
            if not mustpass_rules:
                raise AssertionFailureError(
                    f"--taxonomy-mustpass 指向的对照集无法加载或不含 must_pass "
                    f"条目（{mustpass_path}）。门禁无法执行时不得当作通过。"
                )
            mustpass_scope = mustpass_marker_scope(mustpass_path)
            if mustpass_scope:
                _present = len(mustpass_trees)
                mustpass_trees = {
                    mid: nw for mid, nw in mustpass_trees.items()
                    if mid in set(mustpass_scope)
                }
                logger.info(
                    f"  {canonical_log_tag('2')} must-pass scoped to "
                    f"{len(mustpass_scope)} baseline marker(s); "
                    f"{len(mustpass_trees)} of them present in this run "
                    f"(of {_present} gene trees)."
                )
                if not mustpass_trees:
                    raise AssertionFailureError(
                        "--taxonomy-mustpass 的 markers.ids 与本数据集无交集："
                        "门禁一个标记都没测到，不得当作通过。"
                    )
            else:
                logger.info(
                    f"  {canonical_log_tag('2')} must-pass baseline markers.ids is an unfilled "
                    "skeleton: testing EVERY marker gene tree instead."
                )
            mustpass = evaluate_and_write(
                mustpass_trees,
                taxonomy_map or {},
                mustpass_rules,
                requirements_file=mustpass_path,
                output_dir=self.config.output_dir,
                prefix=self.report_generator.config.output_prefix,
            )
            mustpass_report = mustpass.report
            logger.info(
                f"  {canonical_log_tag('2')} Taxonomy must-pass: {len(mustpass_rules)} rule(s), "
                f"{mustpass_report.markers_checked} tree(s), "
                f"{mustpass_report.groups_tested} group(s) tested, "
                f"{len(mustpass_report.violations)} violation(s), "
                f"{len(mustpass_report.not_checked)} not checked -> "
                f"{mustpass.path.name}"
            )
            enforce_mustpass(mustpass_report)
        else:
            logger.info(
                f"  {canonical_log_tag('2')} Taxonomy must-pass gate NOT ENABLED "
                "(--taxonomy-mustpass unset): broken required relationships "
                "will NOT abort this run."
            )

        excluded = {e.marker_id for e in hgt_report.marker_evaluations if e.level == MarkerLevel.LEVEL_3}
        usable_markers = [m for m in marker_selection.marker_set.markers if m not in excluded]
        marker_genes = {m: marker_sequences.get(m, []) for m in usable_markers}

        return {
            "hgt_report": hgt_report,
            "precomputed_levels": precomputed_levels,
            "far_active": filter_far_active,
            "marker_genes": marker_genes,
        }

    # ── Step 3: infer (Phase 3) ──

    def run_infer(
        self,
        genomes: List[Genome],
        scan_data: dict,
        filter_data: dict,
        context: PhaseContext,
    ) -> dict:
        """Phase 3: phylogenetic inference (concatenation + coalescent).

        Returns dict with key: phylo_result.
        """
        marker_genes = filter_data["marker_genes"]

        logger.info(f"{canonical_log_tag('3')} Phylogenetic inference (concatenation + coalescent)...")
        phylo_result = self.phylo_inference.run(marker_genes, genomes, context)

        # Back-fill phylogenetic informativeness from per-marker gene trees. The
        # Coalescent sub-result holds trees only for markers that passed the HGT
        # Filter, so we also scan the shared gene-trees cache directory (populated
        # By the coalescent builder) to cover every marker that has a gene tree.
        # Coalescent tree-internal support (mean over internal branches) is a direct
        # Proxy for how much phylogenetic signal the marker carries, replacing the
        # Neutral 0.5 placeholder set in Phase 1.5.
        coalescent = getattr(phylo_result, "coalescent", None) or {}
        gene_trees = dict(getattr(coalescent, "gene_trees", {}) or {})

        # The coalescent builder caches each per-marker gene tree on disk. Use that
        # Cache to back-fill informativeness for markers that were excluded by the
        # HGT filter (and thus have no coalescent sub-result tree), so every marker
        # With a gene tree gets a real informativeness value instead of the 0.5 placeholder.
        gene_trees_dir = getattr(getattr(self, "phylo_inference", None), "gene_trees_dir", "") or ""
        if not gene_trees_dir:
            gene_trees_dir = getattr(getattr(getattr(self, "phylo_inference", None), "coalescent", None), "gene_trees_dir", "") or ""
        if gene_trees_dir and Path(gene_trees_dir).is_dir():
            for nwk in Path(gene_trees_dir).glob("*.nwk"):
                mid = nwk.stem
                if mid in gene_trees:
                    continue
                try:
                    from markerfinder.models.tree import Tree as _MFTree
                    gene_trees[mid] = _MFTree.read(str(nwk))
                except Exception:
                    pass

        n_backfilled = 0
        n_pis = 0
        n_support_unmeasured = 0
        composition_rows: dict = {}
        for marker_id, tree in gene_trees.items():
            qs = scan_data.get("marker_selection").quality_scores.get(marker_id)
            if qs is None:
                continue
            _support = _gene_tree_mean_support(tree)
            qs.phylogenetic_informativeness = _support
            if _support is None:
                n_support_unmeasured += 1
            else:
                n_backfilled += 1
            # Back-fill the sequence-informativeness axis from
            # The sidecar written by the gene-tree builder (if present).
            if gene_trees_dir:
                sidecar = Path(gene_trees_dir) / f"{marker_id}.pis"
                if sidecar.exists():
                    try:
                        info = json.loads(sidecar.read_text(encoding="utf-8"))
                        qs.pis = int(info.get("pis", 0))
                        qs.effective_columns = int(info.get("effective_columns", 0))
                        if info.get("rcv") is not None:
                            composition_rows[marker_id] = {
                                "marker_id": marker_id,
                                "rcv": info.get("rcv"),
                                "gc_bias": info.get("gc_bias"),
                                "n_sequences": info.get("n_sequences", "NA"),
                            }
                        n_pis += 1
                    except Exception:
                        pass
        # Apply the provisional PIS floor now that PIS is known.
        # Default floor is 0 = disabled, so shipped behaviour is byte-compatible
        # With the baseline; the flag used to have no consumer at all.
        pis_floor = int(
            getattr(self.config.hgt_config, "min_informative_sites", 0) or 0
        )
        if pis_floor > 0:
            from markerfinder.utils.informative_sites import apply_pis_floor

            # Same trap as the annotation two screens below: ``hgt_report``
            # Is not a local in ``run_infer``, so ``--min-informative-sites`` (a
            # Documented CLI flag) raised NameError before the floor was applied.
            # Nothing in the suite drove Phase 3 with a floor > 0 until the
            # Product-level acceptance tests were written.
            pis_report = filter_data.get("hgt_report")
            if pis_report is None:
                logger.warning(
                    f"  {canonical_log_tag('3')} PIS floor NOT APPLIED: no hgt_report available."
                )
            quality_scores = scan_data.get("marker_selection").quality_scores
            pis_map = {
                marker_id: qs.pis
                for marker_id, qs in quality_scores.items()
                if getattr(qs, "pis", None) is not None
            }
            demoted = apply_pis_floor(
                getattr(pis_report, "marker_evaluations", []) or [],
                pis_map,
                pis_floor,
            )
            if demoted:
                logger.warning(
                    f"  {canonical_log_tag('3')} {len(demoted)} marker(s) below --min-informative-sites"
                    f" {pis_floor} marked inconclusive: "
                    f"{', '.join(sorted(demoted))}"
                )
        # Composition / GC diagnostics go to their own file,
        # As PARALLEL evidence columns. They are never merged into
        # Overall_risk / overall_score (locked by test). Opt-in via
        # --hgt-steps phylogenetic,composition; default off.
        _het_config = (
            getattr(self.config, "heterogeneity_config", None)
            or getattr(self.config, "heterogeneity", None)
        )
        _composition_warn = (
            float(getattr(_het_config, "composition_warn_threshold", 0.15) or 0.15)
            if _het_config is not None else 0.15
        )
        # Two spellings of the same intent used to exist, and only one of
        # Them worked -- see _composition_screen_enabled.
        _composition_requested = _composition_screen_enabled(
            getattr(self.report_generator.config, "composition_screen", False),
            _het_config,
        )
        _warn_about_outlier_metric(_het_config, _composition_requested)
        if _composition_requested:
            from markerfinder.modules.composition import run_stage as run_composition_stage

            # The stage owns its artifact -- the empty case, the file name
            # And the row count are its business, not the orchestrator's.
            _composition = run_composition_stage(
                composition_rows.values(),
                self.config.output_dir,
                self.report_generator.config.output_prefix,
                warn_threshold=_composition_warn,
            )
            if _composition.path is not None:
                logger.info(
                    f"  {canonical_log_tag('3')} Composition diagnostics for "
                    f"{_composition.rows_written} marker(s) written to "
                    f"{_composition.path.name} (parallel columns only)."
                )
        # The cross-framework consistency criterion itself.
        # Runs under --hgt-mode consistency AND hybrid (whose conjunction
        # Semantics -- risk AND consistency, both must pass -- the spec fixes
        # At ); the default risk mode never reaches it, so shipped
        # Products stay identical.
        _hgt_mode = str(
            getattr(self.config.hgt_config, "hgt_mode", "risk") or "risk"
        )
        if _hgt_mode in ("consistency", "hybrid"):
            from markerfinder.modules.consistency_screen import (
                ConsistencyGrade,
                load_cog_category_map,
                run_stage as run_consistency_stage,
                screen,
            )

            def _newick_of(obj):
                if obj is None:
                    return None
                if isinstance(obj, str):
                    return obj or None
                value = getattr(obj, "newick", None)
                if not value and hasattr(obj, "to_newick"):
                    value = obj.to_newick()
                return value or None

            concat_nw = _newick_of(
                getattr(getattr(phylo_result, "supermatrix", None), "tree", None)
            )
            astral_nw = _newick_of(
                getattr(
                    getattr(phylo_result, "coalescent", None),
                    "species_tree", None,
                )
            )
            marker_nw = {}
            for _mid, _tree in gene_trees.items():
                _nw = _newick_of(_tree)
                if _nw:
                    marker_nw[_mid] = _nw
            pis_all = {
                mid: qs.pis
                for mid, qs in scan_data.get(
                    "marker_selection"
                ).quality_scores.items()
                if getattr(qs, "pis", None) is not None
            }

            # ``run_infer`` receives the HGT report through ``filter_data``; it is
            # Not a local here. Reading it from the dict (and saying so when it is
            # Missing) is what the first end-to-end consistency run caught: the
            # Annotation I wired in referenced a name that only exists in
            # ``run_filter``/``run_report``, so ``--hgt-mode consistency`` died with
            # A NameError before writing any product.
            _report_for_notes = filter_data.get("hgt_report")
            if _report_for_notes is None:
                logger.warning(
                    f"  {canonical_log_tag('3')} Consistency per-marker note SKIPPED: no hgt_report in "
                    "filter_data, so the ranking-only demotion could not be "
                    "written into the evaluation notes."
                )

            if not (concat_nw and astral_nw and marker_nw):
                _annotate_risk_role(
                    _report_for_notes,
                    f"fr59: demotion NOT IN EFFECT — {_hgt_mode} criterion could "
                    "not run, markers fall back to risk grading"
                )
                logger.error(
                    f"  {canonical_log_tag('3')} {_hgt_mode} criterion NOT EXECUTED: no usable"
                    f" concatenation tree ({bool(concat_nw)}), ASTRAL tree "
                    f"({bool(astral_nw)}) or per-marker gene trees "
                    f"({len(marker_nw)}). Markers keep their risk grading and "
                    "the consistency column stays NA."
                )
            else:
                consistency_results = screen(
                    marker_nw,
                    concat_nw,
                    astral_nw,
                    stringency=int(
                        getattr(
                            self.config.hgt_config,
                            "consistency_stringency", 1,
                        ) or 1
                    ),
                    min_sites=pis_floor,
                    pis_map=pis_all,
                )
                phylo_result.consistency_grades = {
                    r.marker_id: r.grade.value for r in consistency_results
                }
                if _hgt_mode == "hybrid":
                    # ``hybrid``（两者都要过）— the conjunction joins the
                    # Risk-based ``excluded`` set (run_filter) instead of
                    # Replacing it. The trees of THIS run were built from the
                    # Risk-passing set in a single pass, so the consistency leg
                    # Adjudicates the marker-level verdicts and products below;
                    # Re-running inference on the conjunctive set is a second
                    # Pass this code path does not pretend to have done.
                    hybrid_passed = {
                        r.marker_id for r in consistency_results
                        if r.grade is ConsistencyGrade.CONSISTENT
                    }
                    _annotate_hybrid_verdicts(
                        _report_for_notes,
                        phylo_result.consistency_grades,
                        hybrid_passed,
                    )
                    logger.info(
                        f"  {canonical_log_tag('3')} Hybrid conjunction: risk leg passed "
                        f"{len(marker_nw)} marker(s) into inference; "
                        f"{len(hybrid_passed)} also grade CONSISTENT -> included; "
                        f"{len(marker_nw) - len(hybrid_passed)} excluded by the "
                        "consistency leg (profile below). Single pass: trees of "
                        "this run were built from the risk-passing set."
                    )
                else:
                    _annotate_risk_role(
                        _report_for_notes,
                        "fr59: adjudication = cross-framework consistency "
                        "grade; combined RF+quartet risk is RANKING ONLY"
                    )
                screen_outcome = run_consistency_stage(
                    consistency_results,
                    self.config.output_dir,
                    self.report_generator.config.output_prefix,
                    pis_map=pis_all,
                    category_map=load_cog_category_map(
                        getattr(
                            self.report_generator.config,
                            "cog_category_map", None,
                        )
                    ),
                )
                grade_counts = screen_outcome.grades
                logger.info(
                    f"  {canonical_log_tag('3')} Consistency screen: "
                    f"consistent={grade_counts.get('consistent', 0)}, "
                    f"inconsistent={grade_counts.get('inconsistent', 0)}, "
                    f"inconclusive={grade_counts.get('inconclusive', 0)}; "
                    f"removed-marker profile -> {screen_outcome.profile.path.name}"
                )
        if n_pis:
            logger.info(f"  {canonical_log_tag('3')} Back-filled PIS/effective columns for {n_pis} markers.")
        if n_backfilled or n_support_unmeasured:
            logger.info(
                f"  {canonical_log_tag('3')} Back-filled phylogenetic informativeness for "
                f"{n_backfilled} markers from gene-tree internal support"
                + (
                    f"; {n_support_unmeasured} marker(s) had NO evaluable "
                    f"support and stay unmeasured (None/NA, never 0.5)"
                    if n_support_unmeasured else ""
                )
                + "."
            )

        return {"phylo_result": phylo_result}

    # ── Step 4: report (canonical phase 4) ──

    def run_report(
        self,
        scan_data: dict,
        filter_data: dict,
        infer_data: dict,
        context: PhaseContext,
        runtime_info: RuntimeInfo,
    ) -> dict:
        """Phase 4: report generation.

        Returns dict with keys: html_output, text_output, assertion_report.
        """
        marker_selection = scan_data["marker_selection"]
        hgt_report = filter_data["hgt_report"]
        phylo_result = infer_data["phylo_result"]

        # 's statement is only honest if it reflects the criterion that
        # Actually graded these markers. config_build copies hgt_mode into
        # ReportConfig for CLI runs, but a programmatic or config-file run never
        # Passes through that layer, so the pipeline mirrors it here from the one
        # Authoritative source (the first end-to-end consistency run showed the
        # Summary claiming "risk" while Phase 3 had graded on consistency).
        self.report_generator.config.hgt_mode = str(
            getattr(self.config.hgt_config, "hgt_mode", "risk") or "risk"
        )

        logger.info(f"{canonical_log_tag('4')} Generating reports...")
        html_output, text_output = self.report_generator.run(
            marker_selection, hgt_report, phylo_result, context.quality_data, runtime_info,
            input_dir=self.config.input_dir,
            mode=self.config.mode,
        )

        # Output-numeric assertions run at end of every run,
        # In --check, and in CI. They are observe-only — a failure
        # Aborts the run with EXIT_ASSERTION_FAILED but never re-grades a marker.
        assertion_report = self._run_output_assertions(scan_data, filter_data, infer_data, context)

        # Optional threshold scan — re-grades stored risks,
        # Never re-runs external tools.
        if getattr(self.report_generator.config, "hgt_scan", False):
            from markerfinder.modules.hgt_scan import run_stage as run_scan_stage

            scan_outcome = run_scan_stage(
                hgt_report.marker_evaluations,
                self.config.output_dir,
                self.report_generator.config.output_prefix,
                stability_min=float(
                    getattr(
                        self.report_generator.config, "scan_stability_min", 0.6
                    ) or 0.6
                ),
            )
            logger.info(
                f"  {canonical_log_tag('4')} Threshold scan: "
                f"{len(scan_outcome.results)} band(s) -> {scan_outcome.path.name}"
            )

        return {
            "html_output": html_output,
            "text_output": text_output,
            "assertion_report": assertion_report,
        }

    def _collect_assertion_state(
        self, scan_data: dict, filter_data: dict, infer_data: dict, context: PhaseContext
    ) -> dict:
        """Gather the assertion layer's view of this run."""
        hgt_report = filter_data["hgt_report"]
        evaluations = getattr(hgt_report, "marker_evaluations", []) or []
        steps = [ev.phylo_step for ev in evaluations if ev.phylo_step is not None]

        # A-12: run the built-in zero-discrimination probe (two 4-taxon trees).
        probe = {"conflict_risk": None, "identical_risk": None}
        try:
            from markerfinder.config import HGTConfig
            from markerfinder.models.tree import Tree
            from markerfinder.modules.hgt_filter import (
                HGTDecisionEngine,
                PhylogeneticHGTDetector,
            )

            cfg = HGTConfig()
            det = PhylogeneticHGTDetector(cfg)
            eng = HGTDecisionEngine(cfg)
            ref = Tree(newick="((A,B),(C,D));")
            clash = Tree(newick="((A,C),(B,D));")
            probe["conflict_risk"] = eng.evaluate_marker(
                "probe_clash", det.detect("probe_clash", clash, ref)
            ).overall_risk
            probe["identical_risk"] = eng.evaluate_marker(
                "probe_same", det.detect("probe_same", ref, ref)
            ).overall_risk
        except Exception as e:  # Noqa: BLE001 — probe failure must not crash the run
            logger.debug(f"  zero-discrimination probe skipped: {e}")

        cfg_thresholds = {}
        try:
            cfg_thresholds = dict(getattr(self.config.hgt_config, "level_thresholds", {}) or {})
        except Exception:
            cfg_thresholds = {}

        # A-14: far mode must be visible on every affected marker's notes.
        far_active = bool(filter_data.get("far_active", False))
        far_cards = [
            ev.marker_id for ev in evaluations
            if "far-distance" in (ev.notes or "")
        ] if far_active else []

        return {
            "phylo_steps": steps,
            "hgt_evaluations": [self._assertion_view(ev) for ev in evaluations],
            "evidence_coverage": getattr(hgt_report, "evidence_coverage", None),
            "n_markers": getattr(hgt_report, "total_markers", 0),
            "discrimination_probe": probe,
            "level_thresholds": cfg_thresholds or {"level1_max": 0.25, "level2_max": 0.60},
            "far_active": far_active,
            "far_card_entries": far_cards,
        }

    @staticmethod
    def _assertion_view(ev) -> dict:
        """The decision card is the source the judges read from.

        The assertion state used to read ``ev.overall_risk`` straight off the
        evaluation object — a second, parallel view of the very decision the
        card records. If a back-fill ever updated one and not the other, the
        audit product and the self-check would disagree and nothing would
        notice. Take the value from the card when it carries the key; fall back
        to the attribute only when it does not, and record which was used so a
        divergence stays visible.
        """
        sentinel = object()
        raw_card = getattr(ev, "decision_card", None)
        card = raw_card if isinstance(raw_card, dict) else {}
        risk = sentinel
        for key in ("overall_risk", "risk"):
            if key in card:
                risk = card[key]
                break
        from_card = risk is not sentinel
        if not from_card:
            risk = getattr(ev, "overall_risk", None)
        return {
            "marker_id": card.get("marker_id") or getattr(ev, "marker_id", ""),
            "overall_risk": risk,
            "risk_source": "decision_card" if from_card else "attribute",
        }

    def _run_output_assertions(
        self, scan_data: dict, filter_data: dict, infer_data: dict, context: PhaseContext
    ):
        """Run assertions, persist assertions.tsv, honour strict/bypass flags."""
        from markerfinder.assertions import run_assertions

        state = self._collect_assertion_state(scan_data, filter_data, infer_data, context)
        report = run_assertions(state, mode="run")

        # Persist Phase5_reports/{prefix}.assertions.tsv.
        try:
            from markerfinder.modules.report_generator import write_assertions_tsv

            out_dir = Path(self.config.output_dir) / "Phase5_reports"
            out_dir.mkdir(parents=True, exist_ok=True)
            tsv = out_dir / f"{self.config.output_prefix}.assertions.tsv"
            write_assertions_tsv(report, tsv)
            logger.info(
                f"  Assertions: {report.pass_count} PASS / {report.warn_count} WARN / "
                f"{report.fail_count} FAIL -> {tsv}"
            )
        except Exception as e:  # Noqa: BLE001 — reporting must not mask the run
            logger.warning(
                f"  assertions.tsv NOT WRITTEN — the assertion audit product is "
                f"missing from this run: {e}"
            )

        if report.fail_count:
            failed = [r.assertion_id for r in report.results
                      if not r.passed and r.severity == "fail"]
            if getattr(self.report_generator.config, "allow_assertion_failure", False):
                logger.warning(
                    f"  Assertions FAILED: {failed} — continuing because "
                    "--allow-assertion-failure was set (bypass recorded in report)"
                )
                report.bypassed = True
            else:
                from markerfinder.exceptions import AssertionFailureError

                raise AssertionFailureError(
                    f"output assertions failed: {failed}; "
                    f"re-run with --allow-assertion-failure to bypass "
                    f"(bypass is recorded in the report)"
                )
        return report

    # ── Subcommand support: state persistence ──

    def _cleanup_tmp(self) -> None:
        """清理本趟运行的临时目录.

        仅当 tmp 目录是本进程自动生成的(``tmp_dir_auto=True``, 即用户未显式
        指定 --tmp-dir)且未 ``--keep-tmp`` 时执行: 分步模式 (run_step) 若不清理,
        每个子命令进程都会在系统临时目录泄漏一个 markerfinder-* 目录; 而用户
        显式指定的目录必须保留 —— 聚合 run 也不能清它, 否则与 --tmp-dir 帮助
        文本"显式指定时不自动清理"相悖。两处统一走本助手。
        """
        if not self.config.tmp_dir:
            return
        if not getattr(self.config, "tmp_dir_auto", True):
            return
        if self.config.keep_tmp:
            return
        try:
            shutil.rmtree(self.config.tmp_dir, ignore_errors=True)
            logger.info(f"  Cleaned up tmp dir: {self.config.tmp_dir}")
        except Exception as e:
            logger.debug(f"  tmp cleanup skipped: {e}")

    def _state_path(self) -> Path:
        return state_dir_for(self.config.output_dir) / STATE_FILE_NAME

    def _context_path(self) -> Path:
        return Path(self.config.output_dir) / "Phase5_metadata" / "context.json"

    def _save_state(self, **state) -> None:
        """Persist intermediate state for subcommand chaining (explicit schema).

        Heavy objects are encoded with a typed-JSON codec into individual files
        under ``.markerfinder/objects``; trees are stored as NEWICK files. A
        failure to persist is non-fatal (mirrors the previous pickle behaviour)
        and is logged as a warning rather than crashing the run.
        """
        try:
            save_pipeline_state(self.config.output_dir, dict(state))
        except StateSchemaError as e:
            logger.warning(f"  Failed to save pipeline state: {e}")
            return
        except Exception as e:  # Noqa: BLE001 - persistence must never abort a run
            logger.warning(f"  Failed to save pipeline state: {e}")
            return
        ctx = {
            "completed_steps": state.get("completed_steps", []),
            "timestamp": datetime.now().isoformat(),
            "marker_count": len(state.get("marker_sequences", {}) or {}),
        }
        try:
            self._context_path().parent.mkdir(parents=True, exist_ok=True)
            ctx_target = self._context_path().resolve()
            ctx_target.write_text(
                json.dumps(ctx, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8", newline="\n",
            )
        except OSError:
            pass

    def _load_state(self) -> dict:
        """Load persisted state from previous subcommand (explicit schema)."""
        try:
            state = load_pipeline_state(self.config.output_dir)
        except StateSchemaError as e:
            logger.warning(f"  Failed to load pipeline state: {e}")
            return {}
        if state is None:
            return {}
        return state

    def run_step(
        self,
        step: str,
        genomes: Optional[List[Genome]] = None,
        reference_tree: Optional[object] = None,
        table_taxa: Optional[dict] = None,
        *,
        taxonomy_origin: str = "external table",
    ) -> dict:
        """Run a single pipeline step, loading/saving state for subcommand chaining.

        Args:
            step: one of scan/filter/infer/report.
            genomes: required for scan; loaded from state for subsequent steps.
            reference_tree: optional user-provided reference species tree.
            table_taxa: optional taxonomy table (scan only).
            taxonomy_origin: human-readable origin of ``table_taxa`` (e.g.
                "external table" or "embedded Format A labels"). Used only for
                the Phase 0.2 log line so users can tell how the taxonomy map
                was produced.
        """
        context = self._prepare_context(reference_tree, table_taxa, taxonomy_origin)
        state = self._load_state()

        completed = state.get("completed_steps", [])
        prev_idx = self.STEP_ORDER.index(step) - 1 if step != "scan" else -1
        if prev_idx >= 0:
            prev_step = self.STEP_ORDER[prev_idx]
            if prev_step not in completed:
                raise RuntimeError(
                    f"Step '{prev_step}' must be completed before '{step}'. "
                    f"Completed steps: {completed}"
                )

        # ── Atomicity guard: the previous step's primary product artifact must
        # Exist and be non-empty. This catches corrupt / partially-written state
        # (e.g. a crash mid-step) before we build on top of it, instead of
        # Silently producing dirty downstream results. ──
        prev_product = prev_step_product_rel(step)
        if prev_product is not None:
            artifact = self._state_path().parent / prev_product
            if not artifact.exists() or artifact.stat().st_size == 0:
                raise RuntimeError(
                    f"Step '{step}' cannot start: previous-step product "
                    f"'{prev_product}' is missing or empty. The pipeline state "
                    f"may be corrupt or incomplete; re-run from 'scan' (or the "
                    f"missing step) to rebuild it."
                )

        if step == "scan":
            if genomes is None:
                raise ValueError("genomes is required for the scan step")
            scan_data = self.run_scan(genomes, context)
            state.update(scan_data)
            state["completed_steps"] = ["scan"]
            state["genomes"] = genomes
            state["context"] = context
            state["taxonomy_map"] = table_taxa
            state["start_time"] = time.time()
            self._save_state(**state)
            logger.info("  Step 'scan' completed. Run 'filter' next.")
            self._cleanup_tmp()
            return scan_data

        elif step == "filter":
            genomes = state.get("genomes", genomes)
            scan_data = {k: state[k] for k in
                         ["preprocessing", "marker_selection", "marker_sequences", "rank_map"]
                         if k in state}
            ctx = state.get("context", PhaseContext())
            context.species_tree = ctx.species_tree
            context.quality_data = ctx.quality_data
            context.adaptive_params = ctx.adaptive_params
            context.marker_set = ctx.marker_set
            taxonomy_map = state.get("taxonomy_map")
            filter_data = self.run_filter(genomes, scan_data, context, taxonomy_map=taxonomy_map)
            state.update(filter_data)
            ctx.hgt_evaluations = context.hgt_evaluations
            state["context"] = ctx
            state["completed_steps"] = list(dict.fromkeys(completed + ["filter"]))
            self._save_state(**state)
            logger.info("  Step 'filter' completed. Run 'infer' next.")
            self._cleanup_tmp()
            return filter_data

        elif step == "infer":
            genomes = state.get("genomes", genomes)
            scan_data = {k: state[k] for k in
                         ["preprocessing", "marker_selection", "marker_sequences", "rank_map"]
                         if k in state}
            filter_data = {k: state[k] for k in
                           ["hgt_report", "precomputed_levels", "marker_genes"]
                           if k in state}
            ctx = state.get("context", PhaseContext())
            context.species_tree = ctx.species_tree
            context.hgt_evaluations = ctx.hgt_evaluations
            context.quality_data = ctx.quality_data
            infer_data = self.run_infer(genomes, scan_data, filter_data, context)
            state.update(infer_data)
            state["completed_steps"] = list(dict.fromkeys(completed + ["infer"]))
            self._save_state(**state)
            logger.info("  Step 'infer' completed. Run 'report' next.")
            self._cleanup_tmp()
            return infer_data

        elif step == "report":
            scan_data = {k: state[k] for k in
                         ["preprocessing", "marker_selection", "marker_sequences", "rank_map"]
                         if k in state}
            filter_data = {k: state[k] for k in
                           ["hgt_report", "precomputed_levels", "marker_genes"]
                           if k in state}
            infer_data = {k: state[k] for k in ["phylo_result"] if k in state}
            ctx = state.get("context", PhaseContext())
            context.quality_data = ctx.quality_data
            start_time = state.get("start_time", 0.0)
            end_time = time.time()
            runtime_info = RuntimeInfo(
                start_time=start_time, end_time=end_time,
                duration=end_time - start_time if start_time else 0.0,
            )
            report_data = self.run_report(scan_data, filter_data, infer_data, context, runtime_info)
            if self.config.save_intermediates:
                self._save_intermediates()
            self._write_run_config(start_time, end_time)
            state["completed_steps"] = list(dict.fromkeys(completed + ["report"]))
            self._save_state(**state)
            logger.info("  Step 'report' completed. Pipeline finished.")
            self._cleanup_tmp()
            return report_data

        else:
            raise ValueError(f"Unknown step: {step}. Valid: {self.STEP_ORDER}")

    def _write_run_config(self, start_time: float, end_time: float) -> None:
        """写入 run_config.json，记录完整参数快照和数据库版本哈希，确保可复现性。

        路径参数记录为解析后的绝对路径：快照里一个相对路径无法解释（它是相对于谁
        的？），而重放一份快照恰恰需要知道上一次运行用的是哪里。见
        ``_record_paths_resolved``。
        """
        config_dict = self.config.to_dict()
        config_dict = self._record_paths_resolved(config_dict)

        run_config = {
            "markerfinder_version": get_version_string(),
            "timestamp": datetime.now().isoformat(),
            # start_time 为 0(状态缺失/旧状态)时不能输出天文数字般的
            # "运行时长", 置 0 表示未知.
            "run_duration_seconds": (end_time - start_time) if start_time else 0.0,
            "parameters": config_dict,
            "database_versions": self._collect_database_versions(),
        }

        run_config_path = Path(self.config.output_dir) / "Phase5_metadata" / "run_config.json"
        # 与 report_generator._write_file 一致: 在 no_clobber 模式下保留
        # 已有的 run_config.json(真实运行的可复现性记录), 避免被重复运行
        # 用虚假 duration 污染. 跳过时不补占位时间戳.
        if self.config.no_clobber and run_config_path.exists():
            logger.info(f"  Skipping existing file: {run_config_path}")
            return
        try:
            run_config_path.parent.mkdir(parents=True, exist_ok=True)
            config_target = run_config_path.resolve()
            config_target.write_text(
                json.dumps(run_config, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8", newline="\n",
            )
            logger.info(f"  Written run_config.json to {run_config_path}")
        except OSError as e:
            logger.warning(f"  Failed to write run_config.json: {e}")

    def _record_paths_resolved(self, d: dict) -> dict:
        """Record path fields as the absolute locations this run actually used.

        Storing them relative to ``output_dir`` looked portable but was not
        interpretable: ``--config <previous run>/Phase5_metadata/run_config.json``
        read ``./db`` and ``../../../tmp/...`` against the *current* working
        directory, so a replay reached for a directory nobody asked for (measured
        as ``[Errno 30] Read-only file system`` and as ``db_dir`` pointing inside
        the previous run's output). A snapshot states where the run looked; the
        key set is shared with ``config_loader`` so the reader agrees.
        """
        import os
        from markerfinder.config_loader import SNAPSHOT_PATH_KEYS

        path_keys = set(SNAPSHOT_PATH_KEYS)
        result = {}
        for k, v in d.items():
            if isinstance(v, dict):
                result[k] = self._record_paths_resolved(v)
            elif isinstance(v, str) and k in path_keys and v:
                try:
                    result[k] = os.path.abspath(v)
                except (TypeError, ValueError):
                    result[k] = v
            else:
                result[k] = v
        return result

    def _collect_database_versions(self) -> dict:
        """收集已配置数据库的 SHA256 哈希。"""
        versions: dict = {}
        manager = DatabaseVersionManager(self.config.db_dir)

        # GTDB marker directory (ar53/bac120 HMMs or per-marker FASTAs)
        gtdb_dir = str(Path(self.config.db_dir) / "gtdb_markers")
        if Path(gtdb_dir).exists():
            versions["gtdb_markers"] = manager.verify_database("gtdb_markers", gtdb_dir)

        # HMM library used in hmm mode
        marker_hmm_dir = getattr(self.config.selection_config, "marker_hmm_dir", "") or ""
        if marker_hmm_dir and Path(marker_hmm_dir).exists():
            versions["marker_hmm"] = manager.verify_database("marker_hmm", marker_hmm_dir)

        # A mismatch against db/expected_hashes.json is a WARN
        # On the run path (the run continues, but the fact is recorded).
        for name, info in versions.items():
            if info.get("hash_match") is False:
                logger.warning(
                    f"  Database '{name}' hash mismatch: computed "
                    f"{info['hash'][:12]}... but expected "
                    f"{str(info.get('expected_hash'))[:12]}... — results may not "
                    f"be reproducible against previously published runs."
                )

        return versions

    def _save_intermediates(self) -> None:
        """将 tmp 中有用的中间产物复制到 output/Phase4_intermediate/。

        复制内容:
          - markers/: 每个 marker 的.faa(提取序列)、.aln(MAFFT 比对)、.aln.trim(trimAl 修剪)
          - supermatrix/: 串联后的 concatenated alignment 与 partition 文件
          - quality/: CheckM 质量评估结果(checkm.tsv)与完整 checkm_output/

        仅复制实际存在的文件,并跳过 IQ-TREE/FastTree/ASTRAL 等可随时重建的
        大量日志/检查点文件,以控制 output 体积。
        """
        if not self.config.save_intermediates:
            return
        if not self.config.tmp_dir or not self.config.output_dir:
            return

        tmp_dir = Path(self.config.tmp_dir)
        output_dir = Path(self.config.output_dir)
        intermediate_dir = output_dir / "Phase4_intermediate"

        if not tmp_dir.exists():
            return

        prefix = Path(self.config.output_prefix).name
        copied = 0

        # 1) marker 序列与比对
        markers_dir = intermediate_dir / "markers"
        for ext in (".faa", ".aln", ".aln.trim"):
            for src in tmp_dir.glob(f"*{ext}"):
                # 跳过 supermatrix/concat 前缀的文件
                if src.name.startswith(f"{prefix}_concat") or src.name.startswith(f"{prefix}.concat"):
                    continue
                markers_dir.mkdir(parents=True, exist_ok=True)
                dst = markers_dir / src.name
                try:
                    shutil.copy2(src, dst)
                    copied += 1
                except OSError as e:
                    logger.debug(f"  Failed to copy {src}: {e}")

        # 2) supermatrix 文件
        supermatrix_dir = intermediate_dir / "supermatrix"
        supermatrix_patterns = [
            f"{prefix}.concat.fasta",
            f"{prefix}.partition.nex",
            f"{prefix}_concat.best_model.nex",
            f"{prefix}_concat.best_scheme",
            f"{prefix}_concat.best_scheme.nex",
        ]
        for pattern in supermatrix_patterns:
            src = tmp_dir / pattern
            if src.exists():
                supermatrix_dir.mkdir(parents=True, exist_ok=True)
                try:
                    shutil.copy2(src, supermatrix_dir / src.name)
                    copied += 1
                except OSError as e:
                    logger.debug(f"  Failed to copy {src}: {e}")

        # 3) CheckM 质量结果
        quality_dir = intermediate_dir / "quality"
        checkm_tsv = tmp_dir / "checkm.tsv"
        if checkm_tsv.exists():
            quality_dir.mkdir(parents=True, exist_ok=True)
            try:
                shutil.copy2(checkm_tsv, quality_dir / checkm_tsv.name)
                copied += 1
            except OSError as e:
                logger.debug(f"  Failed to copy {checkm_tsv}: {e}")

        checkm_output_dir = tmp_dir / "checkm_output"
        if checkm_output_dir.exists():
            try:
                dst_checkm = quality_dir / "checkm_output"
                shutil.copytree(checkm_output_dir, dst_checkm, dirs_exist_ok=True)
                copied += 1
            except OSError as e:
                logger.debug(f"  Failed to copy checkm_output: {e}")

        if copied:
            logger.info(f"  Saved {copied} intermediate item(s) to {intermediate_dir}")


    def _resolve_orthologs(
        self,
        genomes: List[Genome],
        marker_selection: MarkerSelectionResult,
        marker_sequences: Dict[str, List],
    ) -> Dict[str, List]:
        # 当前仅 HMM / GTDB-TK 两路径，DIAMOND 直向同源解析预留未启用:
        # 标记序列直接采用, 不做内部直系同源(BBH)解析, 原样返回。保留此薄封装
        # 以便将来接入 DIAMOND 分支时不影响调用方 (pipeline.py L224)。
        logger.info(
            f"{canonical_log_tag('1.5')} Markers taken as-is; internal ortholog resolution "
            "(DIAMOND BBH) not wired in — skipping."
        )
        return marker_sequences


def _composition_screen_enabled(composition_screen: object, heterogeneity: object) -> bool:
    """Either spelling of the intent turns the composition screen on.

    ``hgt_steps: composition`` (which reaches ``ReportConfig.composition_screen``)
    and ``HeterogeneityConfig.enable_composition_screen`` used to be two names for
    one switch where only the first worked -- the shipped MANUAL had to say so out
    loud. Both now count. Defaults (both False) keep the shipped products
    byte-identical.
    """
    return bool(composition_screen) or bool(
        getattr(heterogeneity, "enable_composition_screen", False)
    )


def _warn_about_outlier_metric(heterogeneity: object, enabled: bool) -> None:
    """An unimplemented metric is named, never quietly ignored."""
    if not enabled:
        return
    metric = getattr(heterogeneity, "composition_outlier_metric", "rcv")
    if metric != "rcv":
        logger.warning(
            f"  {canonical_log_tag('3')} Composition outlier metric %r is not implemented; the "
            "only metric is 'rcv'. Using rcv and naming the substitution instead "
            "of ignoring the request.", metric
        )


def _annotate_risk_role(hgt_report: object, note: str) -> None:
    """Mark every evaluation with what actually adjudicated this run.

    ``notes`` reaches both ``hgt_evaluation.tsv`` and the decision card, so the
    demotion of the combined RF+quartet score to ranking-only is per-marker
    evidence rather than one sentence in a summary nobody opens.
    """
    for evaluation in getattr(hgt_report, "marker_evaluations", []) or []:
        existing = getattr(evaluation, "notes", "") or ""
        evaluation.notes = f"{existing}; {note}" if existing else note


def _annotate_hybrid_verdicts(
    hgt_report: object,
    grades: Dict[str, str],
    hybrid_passed: Set[str],
) -> None:
    """Per-marker record of the hybrid conjunction verdict.

    ``hybrid`` means both legs must pass: the risk level screen (Phase 2,
    ``excluded``) AND ``grade == CONSISTENT`` (Phase 3). The note lands in
    ``hgt_evaluation.tsv`` and the decision card, so every marker carries which
    leg decided it -- the same numbers look identical under either criterion
    otherwise (the trap, one mode further).
    """
    for evaluation in getattr(hgt_report, "marker_evaluations", []) or []:
        marker_id = getattr(evaluation, "marker_id", "")
        grade = grades.get(marker_id)
        if getattr(evaluation, "level", None) is MarkerLevel.LEVEL_3:
            verdict = (
                "hybrid: excluded by the risk leg (LEVEL_3); consistency "
                f"grade={grade or 'NA'} cannot rescue — conjunction requires "
                "both."
            )
        elif marker_id in hybrid_passed:
            verdict = (
                "hybrid: passes both legs (risk screen + consistency="
                "consistent)"
            )
        elif grade is None:
            verdict = (
                "hybrid: risk leg passed but no consistency grade (tree "
                "missing or unmeasurable) — NOT included by the conjunction"
            )
        else:
            verdict = (
                f"hybrid: risk leg passed but consistency grade={grade} — "
                "EXCLUDED by the conjunction (both must pass)"
            )
        existing = getattr(evaluation, "notes", "") or ""
        evaluation.notes = f"{existing}; {verdict}" if existing else verdict


def _gene_tree_mean_support(tree: object) -> Optional[float]:
    """Mean internal-branch support of a gene tree as a 0-1 fraction, or None.

    Two independent defects were measured here on ete3 3.1.3 and both are fixed:

    1. **Type contract.** The only caller passes a ``models.tree.Tree`` (built by
       ``_MFTree.read``), which has no ``traverse`` method. Every marker therefore
       fell into the blanket ``except`` and this returned ``None`` for the whole
       run — the PIS / effective-site back-fill never measured anything, and the
       two unit tests asserted exactly that inert contract, so the suite stayed
       green while the feature was dead.
    2. **Scale, plus ete3's default label.** The old code read
       ``node.support`` (which ete3 sets to ``1.0`` when a node carries no label)
       and then clamped with ``min(max(s, 0.0), 1.0)``, mapping IQ-TREE UFBOOT
       percentages onto 1.0. Measured before the fix: a tree with UFBOOT 98/76,
       a tree with UFBOOT 5/3, and a tree with no support labels at all all
       returned exactly ``1.0``. A quantity whose whole purpose is to rank
       markers by informativeness could not tell a 3% branch from a 98% one —
       the zero-discrimination failure mode, re-introduced on the back-fill
       path.

    Labels now come from:meth:`models.tree.Tree.support_labels`, which only
    counts a branch when the Newick really carries a numeric label there. The
    scale is detected rather than assumed: a single label above 1 proves the
    tree is percent-scaled (IQ-TREE UFBOOT/SH), so every label is divided by
    100; a tree whose labels all sit inside 0-1 is read as FastTree SH / ASTRAL
    local support. Labels above 100 are nonsense and yield NOT_MEASURABLE.
    Zero-valued labels are kept: dropping them (as before) biased the mean
    upward.
    """
    def _declare(cause: str) -> None:
        # Name the metric and the cause at run level, not just in a
        # Per-marker NA cell.
        from markerfinder.utils import unmeasured

        unmeasured.record("gene_tree_mean_support (PIS back-fill input)", cause)

    getter = getattr(tree, "support_labels", None)
    if callable(getter):
        labels = [float(v) for v in getter()]
    else:
        newick = getattr(tree, "newick", None)
        if not isinstance(newick, str) or not newick:
            logger.warning(
                "  Gene-tree support NOT_MEASURABLE: unexpected object %s "
                "(no support_labels() and no newick).", type(tree).__name__
            )
            _declare(f"unexpected gene-tree object {type(tree).__name__}")
            return None
        labels = Tree(newick=newick).support_labels()
    if not labels:
        _declare("no readable internal support label on this gene tree")
        return None
    if any(v < 0.0 for v in labels):
        logger.warning(
            "  Gene-tree support NOT_MEASURABLE: negative support label(s) %s.",
            labels,
        )
        _declare("negative support label(s)")
        return None
    if max(labels) > 100.0:
        _declare("support label(s) above 100")
        logger.warning(
            "  Gene-tree support NOT_MEASURABLE: support label(s) above 100 %s.",
            labels,
        )
        return None
    if max(labels) > 1.0:
        # At least one label can only be a percentage, so the whole tree is
        # Percent-scaled (IQ-TREE UFBOOT/SH). A tree whose labels all sit in
        # 0-1 is read as FastTree SH / ASTRAL local support.
        return sum(v / 100.0 for v in labels) / len(labels)
    return sum(labels) / len(labels)
