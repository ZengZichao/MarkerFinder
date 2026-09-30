from __future__ import annotations

import math
import os
import random
import subprocess
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Any

from markerfinder.models.genome import Genome
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.tree import Tree
from markerfinder.models.pipeline_types import (
    PhyloStepResult,
)
from markerfinder.models.pipeline_types import HGTEvaluation, HGTReport
from markerfinder.config import HGTConfig
from markerfinder.taxonomy import AUTO_NO_SCOPE_RANK, AUTO_RANK
from markerfinder.utils.hgt_utils import combine_hgt_scores
from markerfinder.utils.concurrency import parallel_map
import logging

logger = logging.getLogger(__name__)


def _extract_seq(gene_sequence: Any) -> str:
    """从 gene_sequence 中提取蛋白序列字符串, 同时支持 dict 和对象。"""
    if isinstance(gene_sequence, dict):
        seq = gene_sequence.get("seq") or gene_sequence.get("target", "")
    else:
        seq = getattr(gene_sequence, "seq", "")
    if hasattr(seq, "__str__"):
        seq = str(seq)
    return seq if seq else ""


def _extract_id(gene_sequence: Any) -> str:
    """从 gene_sequence 中提取 ID，同时支持 dict 和对象。"""
    if isinstance(gene_sequence, dict):
        return gene_sequence.get("id", gene_sequence.get("target", "unknown"))
    return getattr(gene_sequence, "id", "unknown")


class PhylogeneticHGTDetector:
    """Phylogenetic step: gene-tree vs species-tree (or taxonomy) incongruence detection."""

    def __init__(self, config: HGTConfig):
        self.config = config

    def detect(self, gene_id: str, gene_tree: Tree, species_tree: Tree) -> PhyloStepResult:
        from markerfinder.models.evidence import MeasureState

        # The reference must pass legality checks before it
        # May participate in scoring; an illegal reference routes the marker
        # To UNKNOWN with every component REJECTED — never to a numeric risk.
        from markerfinder.utils.reference import validate_reference_tree

        verdict = validate_reference_tree(species_tree.newick, gene_tree.newick)
        if not verdict.ok:
            reason = f"reference illegal [{verdict.code}]: {verdict.detail}"
            logger.warning(f"  {gene_id}: {reason}")
            return PhyloStepResult(
                gene_id=gene_id,
                rf_distance=None,
                normalized_rf=None,
                quartet_consistency=None,
                rf_state=MeasureState.REJECTED,
                quartet_state=MeasureState.REJECTED,
                rf_reason=reason,
                quartet_reason=reason,
                overall_risk=0.0,
                is_suspicious=False,
                risk_basis="unscreened",
            )

        rf_distance, normalized_rf = self._calculate_rf_distance(gene_tree, species_tree)
        if rf_distance is None:
            # Unmeasurable is UNKNOWN, never the 0.5 placeholder.
            rf_state = MeasureState.NOT_MEASURABLE
            rf_reason = (
                "RF distance uncomputable (ete3 missing and pure-Python "
                "split-set fallback failed or tip cap exceeded)"
            )
            logger.warning(f"  RF distance unmeasurable for {gene_id}: {rf_reason}")
        else:
            rf_state = MeasureState.MEASURED
            rf_reason = ""

        quartet_consistency, quartet_reason = self._calculate_quartet_consistency(
            gene_tree, species_tree
        )
        quartet_state = (
            MeasureState.MEASURED if quartet_consistency is not None
            else MeasureState.NOT_MEASURABLE
        )

        risk: Optional[float]
        risk_basis = ""
        if rf_state is MeasureState.MEASURED and quartet_state is MeasureState.MEASURED:
            risk, _ = combine_hgt_scores(
                {"rf": normalized_rf, "quartet": 1.0 - quartet_consistency},
                getattr(self.config, "hgt_score_weights", {"rf": 0.5, "quartet": 0.5}),
            )
            risk_basis = "two_signal_weighted"
        else:
            # A critical component is unmeasured — no graded risk.
            # Keep the raw measured values on the result; evaluate_marker
            # Routes this marker to UNKNOWN.
            risk = None
            risk_basis = "unscreened"

        return PhyloStepResult(
            gene_id=gene_id,
            rf_distance=rf_distance,
            normalized_rf=normalized_rf,
            quartet_consistency=quartet_consistency,
            rf_state=rf_state,
            quartet_state=quartet_state,
            rf_reason=rf_reason,
            quartet_reason=quartet_reason,
            overall_risk=risk if risk is not None else 0.0,
            is_suspicious=(risk if risk is not None else 0.0) > self.config.phylogenetic_threshold,
            risk_basis=risk_basis,
        )

    def detect_monophyly(
        self,
        marker_id: str,
        gene_tree_newick: str,
        taxonomy_map: Dict[str, Dict],
        level: str = AUTO_RANK,
        threshold: float = 0.5,
    ) -> Optional[PhyloStepResult]:
        """Phylogenetic step via MAD rooting + scope-based monophyly proportion.

        Used when a taxonomy table is supplied instead of comparing the gene
        tree against a species tree. The gene
        tree is first rooted with the MAD criterion (Minimal Ancestor
        Deviation). The taxonomic *scope* of the tree is then inferred from its
        tips: the deepest rank at which every tip shares one label (e.g. all
        tips in the same phylum ⇒ scope = phylum). Under ``level="auto"`` the
        monophyly proportion is computed at the rank **one level below** the
        scope (domain→phylum, phylum→class, class→order, order→family,
        family→genus, genus→species). A low proportion means the marker's gene
        tree is discordant with the taxonomic hierarchy (a hallmark of HGT), so
        the HGT risk is ``1 - proportion``.

        Naming a rank instead of ``auto`` makes that rank authoritative; the
        scope is then only recorded, never applied. Either way the fallback walk
        upward stays, because a rank with no informative split on the tree
        cannot be measured — and under a named rank such a move is flagged as a
        cross-rank comparison rather than silently graded.

        Under ``auto`` the requested rank is therefore only a fallback when the
        tree's scope cannot be determined (e.g. tips span multiple domains); at
        the finest rank (species) no finer rank exists, so the screen is
        degenerate there by construction.

        Args:
            marker_id: marker / gene identifier.
            gene_tree_newick: Newick string of the (unrooted) per-marker tree.
            taxonomy_map: ``{tip: {rank: value,...}}`` mapping (tree tips must
                match the taxonomy keys).
            level: ``"auto"`` (default) to derive the rank from the tree's
                taxonomic scope, or the name of a rank to measure at that rank.
            threshold: monophyly proportion below which the marker is flagged
                as phylogenetically suspicious.

        Returns:
            PhyloStepResult (overall_risk = 1 - proportion) or ``None`` if the
            chosen rank has no taxon with >=2 representatives on the tree.
        """
        from markerfinder.models.evidence import MeasureState

        try:
            from markerfinder.utils.mad_root import mad_root
            from markerfinder.taxonomy import (
                monophyly_proportion,
                detect_scope_rank,
                STANDARD_LEVELS,
            )
            # Ete3 reached only through the sanctioned entry point.
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
        except Exception as e:
            logger.warning(f"  MAD/monophyly unavailable for {marker_id}: {e}")
            # Name the metrics lost and why, at run level too.
            from markerfinder.utils import unmeasured

            unmeasured.record(
                "hgt_monophyly_proportion (MAD-rooted)",
                f"ete3 unusable in this interpreter: {e}",
                marker_id,
            )
            return None

        rooted_nwk = mad_root(gene_tree_newick)

        # Infer the taxonomic scope of this gene tree from its tip labels.
        try:
            t = EteTree(rooted_nwk, format=1)
            tip_ids = [n.name for n in t.get_leaves() if n.name]
        except Exception:
            tip_ids = []
        scope = detect_scope_rank(taxonomy_map, tip_ids)

        # Which rank the proportion is measured at. Two routes, and the user's
        # Word decides between them:
        #
        # * ``level`` names a rank (``--monophyly-rank order``) — that rank is
        # Used. A documented knob that the code replaces with its own automatic
        # Choice is the defect class this project audits (``--max-markers`` was
        # Overridden by the adaptive budget in exactly that way), so an
        # Explicit choice is authoritative.
        # * ``level == "auto"`` (the default) — the rank is derived from the
        # Tree: one below the taxonomic scope of its tips.
        if level != AUTO_RANK:
            test_level = level
            logger.debug(
                f"  {marker_id}: monophyly tested at the requested rank "
                f"{level!r} (scope {scope!r} noted, not applied)"
            )
        elif scope is not None:
            idx = STANDARD_LEVELS.index(scope)
            if idx + 1 < len(STANDARD_LEVELS):
                test_level = STANDARD_LEVELS[idx + 1]
                logger.debug(
                    f"  {marker_id}: tree scope = {scope} -> monophyly tested at "
                    f"'{test_level}' (one rank below scope)"
                )
            else:
                # Scope == species: no finer rank exists; degenerate -> test at
                # Species (trivially monophyletic, no HGT signal at this resolution).
                test_level = "species"
                logger.debug(
                    f"  {marker_id}: tree scope = species (finest rank); "
                    f"testing at 'species' (degenerate)"
                )
        else:
            test_level = AUTO_NO_SCOPE_RANK
            logger.debug(
                f"  {marker_id}: tree scope not determinable (tips span multiple "
                f"domains?) and no rank was requested; measuring at "
                f"{AUTO_NO_SCOPE_RANK!r}"
            )

        # 默认 test_level 无 ≥2 代表时,向更高阶元回退(family→order→...→domain),
        # 直到找到拥有 ≥2 代表的最低阶元. 这避免测试集仅 Haloquadratum 属内有 2 个
        # 代表、其余属全仅 1 的 genomic screen 整体失效.
        fallback_attempts = [test_level] + [lv for lv in reversed(STANDARD_LEVELS)
                                            if lv != test_level]
        test_level_used = test_level
        prop = n_total = n_mono = None
        try:
            from markerfinder.taxonomy import STANDARD_LEVELS as _LV
            fallback_attempts = [test_level] + [lv for lv in reversed(_LV) if lv != test_level]
        except Exception:
            pass

        for attempt_level in fallback_attempts:
            try:
                prop, n_total, n_mono = monophyly_proportion(
                    rooted_nwk, taxonomy_map, attempt_level)
            except Exception as e:
                logger.debug(f"  {marker_id}: monophyly_proportion@{attempt_level} failed: {e}")
                continue
            if prop is not None:
                test_level_used = attempt_level
                if attempt_level != test_level:
                    logger.info(
                        f"  {marker_id}: monophyly rank fallback {test_level!r} -> "
                        f"{attempt_level!r} (no taxon at {test_level} has a "
                        f"informative >=2-vs->=2 split on this tree)"
                    )
                break

        if prop is None:
            logger.info(
                f"  {marker_id}: monophyly screen skipped — no rank from "
                f"{test_level_used!r} up to domain offers an informative split "
                f"(>=2 representatives on both sides of the split)"
            )
            return None

        risk = 1.0 - prop
        logger.info(
            f"  {marker_id}: monophyly@{test_level_used} (scope={scope}) = {prop:.3f} "
            f"({n_mono}/{n_total} taxa monophyletic) -> HGT risk {risk:.3f}"
        )
        # A proportion measured at a rank other than the one the
        # User asked for is NOT comparable with the configured threshold. The
        # Fallback walk stays (robustness), but a mismatch counts as cross-rank
        # Only when the rank was *requested*: under ``auto`` the rank is defined
        # As a function of the tree, so a different rank is not a deviation from
        # A promise. (Flagging the auto case sent every marker of a
        # Scope-detected tree to UNKNOWN and silently disabled the taxonomy
        # Route — measured on the quad4/small8 validation sets.)
        cross_rank = level != AUTO_RANK and test_level_used != level
        return PhyloStepResult(
            gene_id=marker_id,
            rf_distance=None,
            normalized_rf=None,
            quartet_consistency=prop,
            monophyly_proportion=prop,
            monophyly_state=MeasureState.MEASURED,
            overall_risk=risk,
            is_suspicious=risk > threshold,
            risk_basis="monophyly_only",
            test_level_used=test_level_used,
            n_total=n_total,
            n_mono=n_mono,
            cross_rank_comparison=cross_rank,
        )

    def _calculate_rf_distance(self, tree1: Tree, tree2: Tree) -> tuple:
        from markerfinder.utils.tree_utils import calculate_rf_distance
        return calculate_rf_distance(tree1.newick, tree2.newick)

    def _calculate_quartet_consistency(
        self, gene_tree: Tree, species_tree: Tree
    ) -> tuple:
        """计算基因树与物种树之间的四重奏拓扑一致性。

        Returns:
            ``(consistency, reason)`` — ``consistency`` 为 ``None`` 表示不可测
            （不再返回 1.0/0.0 占位），``reason`` 给出可区分的原因。
            ``reason == ""`` 表示测量成功。

        从公共 tip 的全部 quartet 组合空间中无放回随机抽样（固定随机种子，
        结果可复现），最多抽样 1000 个 quartet，以避免固定步长采样的系统性
        偏差并稳定大树的估计。
        """
        tips = set(gene_tree.get_tips()) & set(species_tree.get_tips())
        if len(tips) < 4:
            return None, f"only {len(tips)} shared tips (<4) — no quartet testable"

        try:
            from markerfinder.utils.tree_utils import get_quartet_topology

            tip_list = sorted(tips)
            total_combos = math.comb(len(tip_list), 4)
            max_quartets = min(1000, total_combos)

            rng = random.Random(0)
            sampled_indices = sorted(rng.sample(range(total_combos), max_quartets))

            consistent = 0
            total = 0
            for idx in sampled_indices:
                quartet = self._unrank_quartet(tip_list, idx)
                gene_topo = get_quartet_topology(gene_tree.newick, quartet)
                species_topo = get_quartet_topology(species_tree.newick, quartet)
                if gene_topo is not None and species_topo is not None and gene_topo and species_topo:
                    total += 1
                    if gene_topo == species_topo:
                        consistent += 1

            if total == 0:
                return None, "no quartet topology could be evaluated on either tree"

            return consistent / total, ""

        except ImportError:
            reason = "ete3 unavailable and split-set fallback import failed"
            logger.warning(f"  quartet consistency unmeasurable: {reason}")
            return None, reason
        except Exception as e:
            reason = f"quartet consistency computation failed: {e}"
            logger.warning(f"  {reason}")
            return None, reason

    @staticmethod
    def _unrank_quartet(tip_list: List[str], idx: int) -> tuple:
        """将组合序号映射回 ``tip_list`` 中排序后的 4 元组合。"""
        n = len(tip_list)
        i = 0
        for i in range(n - 3):
            c = math.comb(n - i - 1, 3)
            if idx < c:
                break
            idx -= c
        j = i + 1
        for j in range(i + 1, n - 2):
            c = math.comb(n - j - 1, 2)
            if idx < c:
                break
            idx -= c
        k = j + 1
        for k in range(j + 1, n - 1):
            c = n - k - 1
            if idx < c:
                break
            idx -= c
        l = k + 1 + idx
        return (tip_list[i], tip_list[j], tip_list[k], tip_list[l])


class HGTDecisionEngine:
    """HGT两步筛查综合决策引擎。"""

    def __init__(self, config: HGTConfig):
        self.config = config

    def _unknown(
        self,
        marker_id: str,
        reason_label: str,
        detail: str,
        far_active: bool = False,
    ) -> HGTEvaluation:
        """Single UNKNOWN production path ( principle).

        The phylogenetic HGT screen could not grade this marker. UNKNOWN is
        kept downstream (never excluded, never treated as clean).
        """
        notes = (
            f"Phylogenetic HGT screen could not grade this marker "
            f"({reason_label}: {detail}) — HGT status is UNKNOWN "
            "(not excluded, not clean)."
        )
        return HGTEvaluation(
            marker_id=marker_id,
            overall_risk=0.0,  # Intentional (see trap / test lock)
            level=MarkerLevel.UNKNOWN,
            confidence="unknown",
            step_details={"unscreened": {"reason": reason_label, "detail": detail}},
            notes=notes,
            risk_basis="unscreened",
            decision_card={
                "schema_version": 2,
                "marker_id": marker_id,
                "detector": "none",
                "gene_tree_source": "",
                "trimming_regime": "none",
                "rf": {"value": None, "state": "not_applicable", "method": "", "reason": ""},
                "quartet": {"value": None, "state": "not_applicable", "reason": ""},
                "monophyly": {"value": None, "state": "not_applicable"},
                "rank_used": None,
                "n_total": None,
                "n_mono": None,
                "weights": dict(getattr(self.config, "hgt_score_weights", {}) or {}),
                "n_signals_used": 0,
                "thresholds_in_effect": {},
                "far_active": bool(far_active),
                "risk_basis": "unscreened",
                "overall_risk": 0.0,
                "level": MarkerLevel.UNKNOWN.value,
                "confidence": "unknown",
                "assertion_ids_fired": [],
                "notes": notes,
                "unknown_reason": reason_label,
            },
        )

    def evaluate_marker(
        self,
        marker_id: str,
        phylogenetic_result: Optional[PhyloStepResult] = None,
        silent: bool = False,
        far_active: bool = False,
        skip_reason: str = "",
    ) -> HGTEvaluation:
        """评估单个标记的 HGT 风险并分级 (仅基于系统发育步骤).

        Args:
            marker_id: 标记 id。
            phylogenetic_result: 系统发育步骤结果。
            silent: 为 True 时抑制 far-distance 放宽 INFO 日志,以避免每个
                标记重复打印 far-distance 模式消息(far-distance 检测本身已在
                数据集级别记录一次).最终分档调用应传 True.
            far_active: 显式传入的"远缘数据集"判定结果。由 ``HGTFilterModule.run``
                在数据集级别计算后传入,避免通过 config 私有属性做跨函数隐式通道
                (因此也不会泄漏到后续运行)。保留对 ``config._far_distance_active``
                的读取仅为兼容直接单测该引擎的用例;生产路径一律走此参数。
            skip_reason: 当 ``phylogenetic_result`` 为 None 时的具体成因说明
                (如"未提供参照树"、"基因树构建失败"、"无可检阶元"),用于区分
                UNKNOWN 的不同来源,避免报告中一概归因于"未提供参照"。
        """
        from markerfinder.models.evidence import MeasureState, UnknownReason

        if phylogenetic_result is None:
            reason = skip_reason.strip() or UnknownReason.NO_REFERENCE.value
            detail = skip_reason.strip() or (
                "no reference species tree or taxonomy table was provided"
            )
            return self._unknown(marker_id, reason, detail, far_active=far_active)

        # Partial unmeasurability must not be graded. The
        # Monophyly path is complete by itself; the species-tree path needs
        # Both RF and quartet.
        unmet = phylogenetic_result.unmeasured_components()
        if unmet:
            return self._unknown(
                marker_id,
                unmet[0].value,
                phylogenetic_result.reason_for(unmet[0]),
                far_active=far_active,
            )

        # Cross-rank monophyly proportions must not share the
        # Configured threshold — degrade to UNKNOWN, keeping the measured
        # Values in the notes for auditability.
        if getattr(phylogenetic_result, "cross_rank_comparison", False):
            used = getattr(phylogenetic_result, "test_level_used", None)
            detail = (
                f"monophyly proportion was measured at rank {used!r}, not the "
                f"configured rank — cross-rank proportions are not comparable "
                f"(measured value "
                f"{phylogenetic_result.monophyly_proportion:.3f} kept for audit)"
            )
            ev = self._unknown(
                marker_id, "cross_rank_incomparable", detail, far_active=far_active
            )
            return ev

        risk = phylogenetic_result.overall_risk

        # Thresholds come solely from configuration (HGTConfig.level_thresholds).
        # The fallback reads the canonical defaults defined on HGTConfig, so the
        # Engine no longer carries its own hard-coded 0.25/0.60 magic numbers and
        # There is a single source of truth for HGT level thresholds.
        _defaults = HGTConfig().level_thresholds
        thresholds = self.config.level_thresholds or {}
        level1_max = thresholds.get("level1_max", _defaults.get("level1_max", 0.25))
        level2_max = thresholds.get("level2_max", _defaults.get("level2_max", 0.60))

        # 远缘数据集放宽: 当输入跨多个高阶分类单元时,系统发育风险饱和,
        # 会把所有标记归为 等级 3. 此时放宽 level2_max 以保留足够标记建树.
        # 远缘判定由调用方显式传入 (far_active), 不依赖 config 私有状态,
        # 因此不会泄漏到后续运行. 保留对 _far_distance_active 的读取仅为兼容
        # 直接单测本引擎的用例.
        far_mode = getattr(self.config, "adaptive_far_thresholds", False) and (
            far_active or getattr(self.config, "_far_distance_active", False)
        )
        if far_mode:
            level2_max = max(level2_max, getattr(self.config, "level2_max_far", 0.95))
            if not silent:
                logger.info(f"  [HGT] far-distance mode active for {marker_id}: "
                            f"level2_max relaxed to {level2_max:.2f}")

        step_notes: List[str] = []
        monophyly_only = (
            phylogenetic_result.monophyly_state is MeasureState.MEASURED
        )
        n_signals_used = 2 if not monophyly_only else 1
        if monophyly_only and n_signals_used < 2:
            # A single complete signal (monophyly proportion) is the
            # Legitimate evidence of this path; it grades with capped
            # Confidence (below) rather than via weighted synthesis.
            step_notes.append(
                "graded from monophyly proportion only (single-signal path)"
            )

        if risk < level1_max:
            level = MarkerLevel.LEVEL_1
        elif risk < level2_max:
            level = MarkerLevel.LEVEL_2
        else:
            level = MarkerLevel.LEVEL_3

        # Confidence is NOT a synonym of level. It is derived
        # From the marker's relative distance to the nearest grading boundary
        # (how close the risk is to flipping bands) and is capped for
        # Single-signal evidence paths.
        nearest_boundary = min(abs(risk - level1_max), abs(risk - level2_max))
        relative_distance = nearest_boundary / max(level2_max, 1e-9)
        if relative_distance >= 0.4:
            confidence = "high"
        elif relative_distance >= 0.15:
            confidence = "medium"
        else:
            confidence = "low"
        if monophyly_only and confidence == "high":
            confidence = "medium"
            step_notes.append("confidence capped to medium: single-signal evidence")

        step_notes.append(
            f"risk={risk:.4f} graded against thresholds "
            f"[L1<{level1_max:.2f}, L2<{level2_max:.2f}] "
            f"(n_signals_used={n_signals_used})"
        )
        if far_mode:
            step_notes.append(
                f"far-distance dataset: level2_max relaxed to {level2_max:.2f}"
            )

        risk_basis = phylogenetic_result.risk_basis or (
            "monophyly_only" if monophyly_only else "two_signal_weighted"
        )
        step_details = {
            "phylogenetic": {
                "risk": risk,
                "threshold_exceeded": risk > self.config.phylogenetic_threshold,
                "n_signals_used": n_signals_used,
                "thresholds_in_effect": {
                    "level1_max": level1_max,
                    "level2_max": level2_max,
                },
                "far_active": bool(far_mode),
            }
        }

        notes = " ".join(step_notes)
        evaluation = HGTEvaluation(
            marker_id=marker_id, overall_risk=risk,
            level=level, confidence=confidence, step_details=step_details,
            notes=notes, risk_basis=risk_basis,
            phylo_step=phylogenetic_result,
        )
        # Fill the per-marker decision card.
        res = phylogenetic_result
        rf_state = res.rf_state.value if res is not None else "unknown"
        quartet_state = res.quartet_state.value if res is not None else "unknown"
        monophyly_state = res.monophyly_state.value if res is not None else "unknown"
        evaluation.decision_card = {
            "schema_version": 2,
            "marker_id": marker_id,
            "detector": "taxonomy_monophyly" if monophyly_only else "species_tree_rf",
            "gene_tree_source": "",  # Filled by HGTFilterModule.run (cache info)
            "trimming_regime": "none",
            "rf": {
                "value": res.normalized_rf if res is not None else None,
                "state": rf_state,
                "method": "robinson_foulds",
                "reason": res.rf_reason if res is not None else "",
            },
            "quartet": {
                "value": res.quartet_consistency if res is not None else None,
                "state": quartet_state,
                "reason": res.quartet_reason if res is not None else "",
            },
            "monophyly": {
                "value": res.monophyly_proportion if res is not None else None,
                "state": monophyly_state,
            },
            "rank_used": getattr(res, "test_level_used", None),
            "n_total": getattr(res, "n_total", None),
            "n_mono": getattr(res, "n_mono", None),
            "weights": dict(getattr(self.config, "hgt_score_weights", {}) or {}),
            "n_signals_used": n_signals_used,
            "thresholds_in_effect": {
                "level1_max": level1_max,
                "level2_max": level2_max,
            },
            "far_active": bool(far_mode),
            "risk_basis": risk_basis,
            "overall_risk": risk,
            "level": level.value,
            "confidence": confidence,
            "assertion_ids_fired": [],
            "notes": notes,
        }
        return evaluation

    def generate_hgt_report(self, evaluations: List[HGTEvaluation], phylogenetic_enabled: bool = True,
                            far_active: bool = False) -> HGTReport:
        from markerfinder.models.evidence import UnknownReason

        level_counts = Counter(e.level for e in evaluations)
        n_measured = sum(1 for e in evaluations if e.level is not MarkerLevel.UNKNOWN)
        total = len(evaluations)
        reason_counts: Dict[str, int] = Counter()
        for e in evaluations:
            if e.level is MarkerLevel.UNKNOWN and e.step_details.get("unscreened"):
                label = e.step_details["unscreened"].get("reason", "unknown")
                try:
                    label = UnknownReason(label).value
                except ValueError:
                    pass
                reason_counts[label] += 1
        coverage = (n_measured / total) if total else 0.0
        return HGTReport(
            total_markers=total,
            level1_count=level_counts.get(MarkerLevel.LEVEL_1, 0),
            level2_count=level_counts.get(MarkerLevel.LEVEL_2, 0),
            level3_count=level_counts.get(MarkerLevel.LEVEL_3, 0),
            unknown_count=level_counts.get(MarkerLevel.UNKNOWN, 0),
            mean_risk=(sum(e.overall_risk for e in evaluations) / len(evaluations) if evaluations else 0.0),
            high_confidence_count=sum(1 for e in evaluations if e.confidence == "high"),
            marker_evaluations=evaluations,
            phylogenetic_enabled=phylogenetic_enabled,
            n_actually_measured=n_measured,
            evidence_coverage=coverage,
            unknown_reason_counts=dict(reason_counts),
            far_active=far_active,
        )


class HGTFilterModule:
    """路径二: HGT感知标记过滤主模块。"""

    def __init__(self, config: HGTConfig, tmp_dir: str = "/tmp", gene_trees_dir: str = "", cpus: int = 1):
        self.config = config
        # Hold run-time cpus on the instance, NOT on the shared config object,
        # So we never mutate caller-owned configuration state.
        self.cpus = cpus
        self.tmp_dir = tmp_dir
        self.gene_trees_dir = gene_trees_dir or os.path.join(tmp_dir, "gene_trees")
        self.decision_engine = HGTDecisionEngine(config)

    def run(
        self,
        marker_ids: List[str],
        genomes: List[Genome] = None,
        marker_sequences: Optional[Dict[str, List]] = None,
        occupancy_ranking: Optional[Dict[str, int]] = None,
        precomputed_levels: Optional[Dict[str, MarkerLevel]] = None,
        prebuilt_gene_trees: Optional[Dict[str, str]] = None,
        taxonomy_map: Optional[Dict] = None,
        species_tree: Optional[object] = None,
        monophyly_rank: str = AUTO_RANK,
        monophyly_threshold: float = 0.5,
    ) -> tuple:
        """执行 HGT 筛查 (仅系统发育步骤).

        对每个标记基因构建基因树, 基于单系比例或物种树不一致性评估
        HGT 风险并分级.

        Returns:
            (HGTReport, Dict[str, MarkerLevel]) — 报告和预计算 等级 映射
        """
        evaluations: List[HGTEvaluation] = []
        computed_levels: Dict[str, MarkerLevel] = {}

        # 远缘判定: 基于 taxonomy_map 的高阶分类单元多样性. 当输入跨多个属/目时,
        # 系统发育风险饱和,放宽 level2_max 以保留标记.
        far_active = False
        if taxonomy_map and getattr(self.config, "adaptive_far_thresholds", False):
            try:
                from markerfinder.taxonomy import parse_taxonomy
                parsed = []
                for tip, tax in taxonomy_map.items():
                    if isinstance(tax, dict):
                        parsed.append(tax)
                    else:
                        parsed.append(parse_taxonomy(tax))
                genera = {t.get("genus") for t in parsed if t.get("genus")}
                orders = {t.get("order") for t in parsed if t.get("order")}
                n_genomes = len(taxonomy_map)
                ratio = len(genera) / n_genomes if n_genomes else 0.0
                if ratio >= getattr(self.config, "far_distance_genera_ratio", 0.70) or \
                   len(orders) >= getattr(self.config, "far_distance_min_orders", 2):
                    far_active = True
                    logger.info(f"  [HGT] far-distance dataset detected: {len(genera)} genera, "
                                f"{len(orders)} orders across {n_genomes} genomes; relaxing level2_max "
                                f"to {getattr(self.config, 'level2_max_far', 0.95):.2f}")
            except Exception as e:
                logger.debug(f"  [HGT] far-distance detection skipped: {e}")

        # 先并发构建所有需要的基因树 (每个标记独立的 MAFFT+FastTree, 写入各自
        # {marker_id}.nwk, 文件名互斥, 可安全并发), 再串行做检测/分级. 结果与
        # 原逐标记串行构建完全一致.
        need_tree: Dict[str, List] = {}
        if self.config.enable_phylogenetic and (species_tree is not None or taxonomy_map):
            for marker_id in marker_ids:
                seqs = (marker_sequences or {}).get(marker_id, [])
                if seqs:
                    need_tree[marker_id] = seqs

        gene_tree_source: Dict[str, str] = {}

        def _build_task(task):
            mid, seqs = task
            if prebuilt_gene_trees and prebuilt_gene_trees.get(mid):
                try:
                    tree = Tree.read(prebuilt_gene_trees[mid])
                    gene_tree_source[mid] = "prebuilt"
                    return mid, tree
                except Exception:
                    pass
            tree, source = _hgt_build_gene_tree_with_source(
                mid, seqs, self.tmp_dir, self.cpus, self.gene_trees_dir,
            )
            gene_tree_source[mid] = source
            return mid, tree

        built = parallel_map(
            _build_task, list(need_tree.items()),
            max_workers=max(1, min(self.cpus, len(need_tree))),
        )
        gene_tree_cache: Dict[str, Optional[Tree]] = {mid: t for mid, t in built}

        for marker_id in marker_ids:
            phylogenetic_res = None
            skip_reason = ""
            if not self.config.enable_phylogenetic:
                skip_reason = "the phylogenetic HGT step is disabled (--hgt-steps)"
            elif species_tree is None and taxonomy_map is None:
                # 保持默认成因文案: 未提供任何参照.
                skip_reason = ""
            else:
                # HGT 系统发育步骤二选一:
                # 1) 用户提供了参考物种树 -> 基因树 vs 物种树 (RF + quartet)
                # 2) 未提供物种树 -> MAD 定根 + 单系比例 (需 taxonomy_map)
                marker_seqs = (marker_sequences or {}).get(marker_id, [])
                gene_tree = gene_tree_cache.get(marker_id)
                if gene_tree is None and (species_tree is not None or taxonomy_map):
                    gene_tree = self._get_gene_tree(marker_id, marker_seqs, prebuilt_gene_trees)
                if gene_tree is None:
                    skip_reason = (
                        "the marker gene tree could not be built "
                        "(<4 sequences, or MAFFT/FastTree unavailable or failed)"
                    )
                else:
                    detector = PhylogeneticHGTDetector(self.config)
                    if species_tree is not None:
                        phylogenetic_res = detector.detect(
                            marker_id, gene_tree, species_tree
                        )
                    else:
                        phylogenetic_res = detector.detect_monophyly(
                            marker_id,
                            gene_tree.newick,
                            taxonomy_map,
                            monophyly_rank,
                            monophyly_threshold,
                        )
                        if phylogenetic_res is None:
                            skip_reason = (
                                "no taxonomic rank on this tree offers an "
                                "informative monophyly split (>=2 representatives "
                                "on both sides); monophyly screen not applicable"
                            )

            evaluation = self.decision_engine.evaluate_marker(
                marker_id, phylogenetic_res, silent=True, far_active=far_active,
                skip_reason=skip_reason,
            )
            evaluations.append(evaluation)
            computed_levels[marker_id] = evaluation.level

        # Fill tree provenance into each decision card.
        # Occupancy joins the card as an INDEPENDENT
        # Evidence column (never merged into overall_risk).
        for ev in evaluations:
            if ev.decision_card:
                ev.decision_card["gene_tree_source"] = gene_tree_source.get(ev.marker_id, "")
                ev.decision_card["trimming_regime"] = "none"  # Phase 2 HGT trees are untrimmed
                rank = (occupancy_ranking or {}).get(ev.marker_id)
                ev.decision_card["occupancy"] = rank

        report = self.decision_engine.generate_hgt_report(
            evaluations,
            phylogenetic_enabled=self.config.enable_phylogenetic,
            far_active=far_active,
        )
        # Expose the far-distance decision for the assertion layer (A-14).
        return report, computed_levels, far_active

    def _get_gene_tree(
        self,
        marker_id: str,
        marker_seqs: List,
        prebuilt_gene_trees: Optional[Dict[str, str]],
    ) -> Optional[Tree]:
        """Return a per-marker gene tree, preferring a pre-built one.

        A pre-built tree (e.g. from GTDB-TK alignments, already MAFFT-aligned,
        trimmed and FastTree-built) is read when available; otherwise one is
        built on the fly from the marker sequences.

        Returns:
            Tree object, or None if no tree could be obtained.
        """
        if prebuilt_gene_trees and marker_id in prebuilt_gene_trees and prebuilt_gene_trees[marker_id]:
            try:
                return Tree.read(prebuilt_gene_trees[marker_id])
            except Exception:
                pass
        return self._build_gene_tree(marker_id, marker_seqs)

    def _build_gene_tree(self, marker_id: str, marker_seqs: List) -> Optional[Tree]:
        """为指定标记基因构建快速系统发育树 (FastTree)，支持缓存复用。

        构建的树写入 {gene_trees_dir}/{marker_id}.nwk 持久化缓存，
        供合并法 (CoalescentInference, canonical phase 3) 复用，避免重复建树。
        (历史注释误写 "Phase 4"；目录名 Phase4_trees 属历史别名, 保留)
        """
        return _hgt_build_gene_tree(
            marker_id, marker_seqs, self.tmp_dir, self.cpus, self.gene_trees_dir,
        )


def _hgt_build_gene_tree(
    marker_id: str,
    marker_seqs: List,
    tmp_dir: str,
    cpus: int,
    gene_trees_dir: str,
) -> Optional[Tree]:
    """为指定标记基因构建快速系统发育树 (FastTree)，支持缓存复用。"""
    tree, _source = _hgt_build_gene_tree_with_source(
        marker_id, marker_seqs, tmp_dir, cpus, gene_trees_dir,
    )
    return tree


def _marker_input_id(marker_seqs: List) -> str:
    """Sha256 over the (id, sequence) pairs a gene tree would be built from.

    Ordering-independent, so a marker file that is re-written with its records
    in another order still identifies the same input.
    """
    import hashlib

    digest = hashlib.sha256()
    for seq in sorted(marker_seqs, key=lambda s: str(_extract_id(s))):
        digest.update(f"{_extract_id(seq)}\0{_extract_seq(seq)}\n".encode("utf-8"))
    return digest.hexdigest()


def _hgt_build_gene_tree_with_source(
    marker_id: str,
    marker_seqs: List,
    tmp_dir: str,
    cpus: int,
    gene_trees_dir: str,
) -> tuple:
    """Source-aware worker: returns ``(tree, source)`` where
    ``source`` ∈ {"cached", "rebuilt", "prebuilt", ""} ("" = no tree)."""
    if not marker_seqs or len(marker_seqs) < 4:
        return None, ""

    cache_dir = Path(gene_trees_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{marker_id}.nwk"
    stamp_path = cache_dir / f"{marker_id}.nwk.input_sha256"
    input_id = _marker_input_id(marker_seqs)

    # 缓存命中：仅当缓存树确实是从这份输入构建出来的才算命中。
    # A cache entry keyed by marker id alone is handed back to whatever run next
    # Uses this output directory (``--force`` over an existing directory is a
    # Documented workflow), and a tree built from different sequences would then
    # Silently become this run's input. The id of the sequences is stored beside
    # The tree and a mismatch rebuilds it.
    if cache_path.exists() and cache_path.stat().st_size > 0:
        stamped = ""
        if stamp_path.exists():
            stamped = stamp_path.read_text(encoding="utf-8").strip()
        if stamped != input_id:
            logger.info(
                f"  {marker_id}: cached gene tree was built from different "
                "sequences; rebuilding"
            )
        else:
            try:
                return Tree.read(str(cache_path)), "cached"
            except Exception:
                pass  # 缓存损坏，重建

    fasta_path = os.path.join(tmp_dir, f"hgt_l3_{marker_id}.faa")
    aln_path = os.path.join(tmp_dir, f"hgt_l3_{marker_id}.aln")

    try:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)

        fasta_lines: List[str] = []
        for s in marker_seqs:
            seq_str = _extract_seq(s)
            sid = _extract_id(s)
            fasta_lines.append(f">{sid}\n{seq_str}\n")
        fasta_target = Path(fasta_path).resolve()
        fasta_target.write_text("".join(fasta_lines), encoding="utf-8", newline="\n")

        aln_target = Path(aln_path).resolve()
        with aln_target.open("w", encoding="utf-8", newline="\n") as out:
            subprocess.run(
                ["mafft", "--auto", "--quiet", "--thread", "1", fasta_path],
                stdout=out, stderr=subprocess.DEVNULL,
                check=True, timeout=300,
            )

        # FastTree 二进制名在各平台/安装中大小写不统一(FastTree / fasttree),
        # 依次尝试,首个命中即使用;蛋白序列使用 WAG 模型.
        ft_ok = False
        for ft_exe in ("FastTree", "fasttree"):
            cache_target = Path(cache_path).resolve()
            try:
                if ft_exe == "FastTree":
                    subprocess.run(
                        ["FastTree", "-wag", "-quiet", "-seed", "1", "-out", str(cache_path), aln_path],
                        capture_output=True, check=True, timeout=300,
                    )
                else:
                    subprocess.run(
                        ["fasttree", "-wag", "-quiet", "-seed", "1", "-out", str(cache_path), aln_path],
                        capture_output=True, check=True, timeout=300,
                    )
                ft_ok = True
                break
            except FileNotFoundError as e:
                logger.debug(f"  FastTree binary {ft_exe} not found: {e}")
                continue
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
                break
        if not ft_ok:
            logger.debug(f"  FastTree failed for {marker_id}")

        if cache_path.exists() and cache_path.stat().st_size > 0:
            stamp_path.write_text(input_id + "\n", encoding="utf-8")
            return Tree.read(str(cache_path)), "rebuilt"
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass
    finally:
        for p in [fasta_path, aln_path]:
            Path(p).unlink(missing_ok=True)

    return None, ""
