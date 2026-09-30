from __future__ import annotations

import json
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from markerfinder.phases import canonical_log_tag
from markerfinder.config import ReportConfig
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.pipeline_types import (
    HGTReport,
    MarkerSelectionResult,
    PhylogeneticResult,
    QualityData,
)
from markerfinder.models.report import RuntimeInfo, ReportOutput, PlainTextReportOutput, DownloadLink
from markerfinder.utils.tree_source import (
    coalescent_is_degenerate,
    has_usable_tree,
    prioritize_tree_source,
)

import logging

logger = logging.getLogger(__name__)


# NOTE: 物种树来源优先级判定已收敛到 markerfinder.utils.tree_source
# (prioritize_tree_source / coalescent_is_degenerate)。本模块仅引用公共实现,
# 不再各自维护一份相同的优先级逻辑, 避免与 pipeline 的可用性信号判定漂移.


def _fmt(value) -> str:
    """Render an optional metric for TSV output.

    ``None`` means "not measured" and MUST render as ``NA`` — never as
    ``0.0000``, which a reader cannot distinguish from a real zero.
    """
    if value is None:
        return "NA"
    return f"{value:.4f}"


def _fmt_html(value) -> str:
    """ HTML counterpart of:func:`_fmt` (missing value renders N/A)."""
    if value is None:
        return "N/A"
    return f"{value:.4f}"


def _plain(value) -> str:
    """Render a scalar card value for TSV (None -> NA, no float formatting)."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_assertions_tsv(report, tsv_path) -> None:
    """Persist the assertion report as Phase5_reports/{prefix}.assertions.tsv."""
    lines = ["\t".join(row) + "\n" for row in report.to_tsv_rows()]
    target = Path(tsv_path).resolve()
    target.write_text("".join(lines), encoding="utf-8", newline="\n")


class InteractiveReportGenerator:
    def __init__(self, config: ReportConfig, force: bool = False, no_clobber: bool = False):
        self.config = config
        self.force = force
        self.no_clobber = no_clobber

    def generate(
        self,
        marker_result: MarkerSelectionResult,
        hgt_report: HGTReport,
        phylo_result: PhylogeneticResult,
        quality_data: Optional[QualityData],
        runtime_info: RuntimeInfo,
    ) -> ReportOutput:
        logger.info(f"{canonical_log_tag('4')} Generating interactive HTML report")
        prefix_name = Path(self.config.output_prefix).name
        output_path = f"Phase5_reports/{prefix_name}.report.html"
        Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)
        (Path(self.config.output_dir) / "Phase5_reports").mkdir(parents=True, exist_ok=True)

        html = self._build_html(marker_result, hgt_report, phylo_result, quality_data, runtime_info)
        full_path = str(Path(self.config.output_dir) / output_path)

        if self.no_clobber and Path(full_path).exists():
            logger.info(f"  Skipping existing file: {full_path}")
            return ReportOutput(html_path=full_path, data_files=[])

        html_target = Path(full_path).resolve()
        html_target.write_text(html, encoding="utf-8", newline="\n")

        return ReportOutput(html_path=full_path, data_files=[])

    def _build_html(self, marker_result, hgt_report, phylo_result, quality_data, runtime_info) -> str:
        prefix_name = Path(self.config.output_prefix).name
        n_markers = len(marker_result.marker_set.markers) if marker_result.marker_set else 0
        n_genomes = marker_result.occupancy_matrix.n_genomes if marker_result.occupancy_matrix else 0

        # Bootstrap 饱和提示: 小数据集(标记数 < 12) + 内部分支 UFBOOT 全为 100 时,
        # 不要把 100 过度解读(小数据集 IQ-TREE3 UFBOOT 天然趋向饱和).
        bootstrap_note = ""
        concat_tree = getattr(getattr(phylo_result, "supermatrix", None), "tree", None)
        avg_support = None
        if concat_tree is not None and hasattr(concat_tree, "get_average_support"):
            try:
                avg_support = concat_tree.get_average_support()
            except Exception:  # Noqa: BLE001 - a note, never a crash
                avg_support = None
        # ``None`` means the supports were not readable (e.g. no ete3 in this
        # Interpreter). Unmeasurable is not 0.0 and is definitely not "saturated"
        #, so the heuristic stays silent rather than guessing.
        if (avg_support is not None
                and n_markers > 0 and n_markers < 12 and avg_support >= 99.99):
            bootstrap_note = (
                "<p><strong>Note:</strong> concatenation-tree UFBOOT support values are "
                "saturated (all ~100). This is common for small datasets and should not be "
                "over-interpreted as exceptionally strong phylogenetic signal.</p>"
            )

        # 真实可用物种树的优先级: astral/consensus (真合并物种树) >
        # Concat (supermatrix) > first_gene_tree(退化替代) > none.
        # 关键: 'first_gene_tree' 是 ASTRAL-III 缺失时的退化替代,可靠性低于
        # IQ-TREE3 串联法 concat,因此仅在无 concat 时才作为兜底展示.
        coalescent_source = getattr(getattr(phylo_result, "coalescent", None),
                                   "species_tree_source", None)
        has_coalescent_tree = bool(
            coalescent_source and coalescent_source != "none"
            and has_usable_tree(getattr(getattr(phylo_result, "coalescent", None), "species_tree", None))
        )
        has_concat_tree = has_usable_tree(
            getattr(getattr(phylo_result, "supermatrix", None), "tree", None)
        )
        tree_source = prioritize_tree_source(
            coalescent_source if has_coalescent_tree else None,
            has_concat_tree,
        )
        # 退化合并树 → 在 HTML 中的展示名明确标识
        if has_coalescent_tree and coalescent_is_degenerate(coalescent_source):
            tree_source_display = f"{tree_source} [fallback; install ASTRAL-III for a proper coalescent tree]"
        else:
            tree_source_display = tree_source

        # HGT 汇总表(按 等级 分组). 必须包含 "unknown": 筛查未运行时标记为
        # MarkerLevel.UNKNOWN(被保留而非剔除), 缺桶会直接 KeyError 崩溃.
        level_rows: Dict[str, List[str]] = {"level_1": [], "level_2": [], "level_3": [], "unknown": []}
        for ev in hgt_report.marker_evaluations:
            level_rows.setdefault(ev.level.value, []).append(ev.marker_id)
        level_counts = {t: len(rows) for t, rows in level_rows.items()}

        def _table(title: str, rows: List[List[str]], headers: List[str]) -> str:
            if not rows:
                return f"<h2>{title}</h2><p><em>None.</em></p>"
            out = [f"<h2>{title}</h2>",
                   "<table><thead><tr>"]
            out += [f"<th>{h}</th>" for h in headers]
            out.append("</tr></thead><tbody>")
            for row in rows:
                out.append("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>")
            out.append("</tbody></table>")
            return "".join(out)

        # Quality summary from CheckM / default estimates.
        quality_source = getattr(quality_data, "quality_source", "unknown") if quality_data else "unknown"
        quality_rows = []
        if quality_data and quality_data.genome_qualities:
            for gq in quality_data.genome_qualities:
                quality_rows.append([
                    gq.genome_id if hasattr(gq, "genome_id") else "",
                    f"{gq.completeness:.2f}",
                    f"{gq.contamination:.2f}",
                    f"{gq.quality_score:.2f}",
                ])
        quality_table = _table(
            "Quality Summary", quality_rows,
            ["Genome", "Completeness (%)", "Contamination (%)", "Quality score"],
        )
        quality_note = f"<p>Quality source: <code>{quality_source}</code></p>"

        # Species tree files: list both concat and coalescent outputs.
        species_tree_rows = []
        if has_concat_tree:
            species_tree_rows.append([
                "Concatenation (supermatrix)",
                f"Phase4_trees/{prefix_name}.species_tree_concat.newick",
                "yes" if tree_source == "concat (supermatrix)" else "",
            ])
        if has_coalescent_tree:
            species_tree_rows.append([
                f"Coalescent ({coalescent_source})",
                f"Phase4_trees/{prefix_name}.species_tree_astral.newick",
                "yes" if tree_source == coalescent_source else "",
            ])
        species_tree_table = _table(
            "Species Trees", species_tree_rows,
            ["Method", "File", "Primary"],
        )

        # Conflict detection summary.
        conflict_report = getattr(phylo_result, "conflict_report", None)
        rf_str = "N/A"
        quartet_agreement_str = "N/A"
        topology_note_str = ""
        if conflict_report is not None and getattr(conflict_report, "normalized_rf", None) is not None:
            rf_str = f"{conflict_report.normalized_rf:.4f}"
        if conflict_report is not None and getattr(conflict_report, "quartet_agreement", None) is not None:
            quartet_agreement_str = f"{conflict_report.quartet_agreement:.4f}"
        if conflict_report is not None:
            topology_note_str = getattr(conflict_report, "topology_note", "") or ""
        conflict_rows = [
            ["Normalized RF distance (concat vs coalescent)", rf_str],
            ["Quartet agreement (concat vs coalescent)", quartet_agreement_str],
        ]
        if topology_note_str:
            conflict_rows.append(["Topology interpretation", topology_note_str])
        conflict_table = _table(
            "Conflict Detection", conflict_rows,
            ["Metric", "Value"],
        )

        # HGT step availability note.
        hgt_step_notes: List[str] = []
        if not hgt_report.phylogenetic_enabled:
            hgt_step_notes.append("<strong>Phylogenetic step disabled</strong>.")
        else:
            hgt_step_notes.append("Phylogenetic HGT step was enabled.")
        hgt_availability_html = "<p>" + "<br>".join(hgt_step_notes) + "</p>"

        # IQ-TREE composition test note.
        iqtree_note = (
            "<p><strong>Note on IQ-TREE composition test:</strong> IQ-TREE may report "
            "'failed composition chi2 test' for GC/codon-specialised clades (e.g. halophilic archaea). "
            "This warning is typically benign and does not invalidate the inferred tree.</p>"
        )

        # HGT risk distribution bar chart (inline CSS).
        # 分桶边界取自配置(与实际分级同源, 见 ReportConfig.level1_max/level2_max),
        # 不再硬编码全局常量; UNKNOWN 标记未经过风险度量(overall_risk 恒 0),
        # 单独以灰色条呈现, 不计入 Level 1 桶避免误导.
        l1_max = float(getattr(self.config, "level1_max", 0.25))
        l2_max = float(getattr(self.config, "level2_max", 0.60))
        bucket_l1 = bucket_l2 = bucket_l3 = bucket_unknown = 0
        for ev in hgt_report.marker_evaluations:
            if ev.level == MarkerLevel.UNKNOWN:
                bucket_unknown += 1
            elif ev.overall_risk < l1_max:
                bucket_l1 += 1
            elif ev.overall_risk < l2_max:
                bucket_l2 += 1
            else:
                bucket_l3 += 1
        total_eval = len(hgt_report.marker_evaluations) or 1
        bar_html = '<div style="margin:10px 0;">'
        risk_bars = [
            (f"0.00 - {l1_max:.2f} (Level 1)", bucket_l1, "#2ecc71"),
            (f"{l1_max:.2f} - {l2_max:.2f} (Level 2)", bucket_l2, "#f1c40f"),
            (f"≥ {l2_max:.2f} (Level 3)", bucket_l3, "#e74c3c"),
        ]
        if bucket_unknown:
            risk_bars.append(("unknown (screen skipped)", bucket_unknown, "#95a5a6"))
        for label, count, color in risk_bars:
            pct = count / total_eval * 100
            bar_html += (
                f'<div style="margin:4px 0;">'
                f'<span style="display:inline-block;width:180px;font-size:12px;">{label}</span>'
                f'<span style="display:inline-block;width:{pct:.1f}%;background:{color};" title="{count} markers">&nbsp;</span>'
                f'<span style="margin-left:8px;font-size:12px;">{count}</span>'
                f'</div>'
            )
        bar_html += "</div>"
        risk_distribution_html = f"<h2>HGT Risk Distribution</h2>{bar_html}"

        hgt_rows = [[ev.marker_id, f"{ev.overall_risk:.4f}", ev.level.value, ev.confidence, ev.notes]
                    for ev in hgt_report.marker_evaluations]
        hgt_table = _table(
            "HGT Evaluation", hgt_rows,
            ["Marker", "Overall risk", "HGT risk level", "HGT evidence confidence", "Notes"],
        )

        marker_rows = []
        if marker_result.marker_set:
            for mid in marker_result.marker_set.markers:
                occ = marker_result.marker_set.occupancy_scores.get(mid, 0)
                qs = marker_result.quality_scores.get(mid)
                if qs:
                    q_score = _fmt_html(qs.overall_score)
                    q_level = qs.level.value
                else:
                    q_score, q_level = "NA", "NA"
                marker_rows.append([mid, f"{occ:.4f}", q_score, q_level])
        marker_table = _table(
            "Marker Summary", marker_rows,
            ["Marker", "Occupancy", "Marker quality score", "Marker quality level"],
        )

        return f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>MarkerFinder Report</title>
<style>body{{font-family:Arial;margin:20px}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ddd;padding:8px}}th{{background:#3498db;color:#fff}}h2{{margin-top:30px;color:#2c3e50}}.note{{background:#f9f9f9;border-left:4px solid #3498db;padding:10px;margin:10px 0}}</style>
</head><body>
<h1>MarkerFinder Report</h1>
<p>Generated: {datetime.now().isoformat()}</p>
<p>Runtime: {runtime_info.duration:.1f}s | Genomes: {n_genomes} | Markers: {n_markers}</p>
<p>HGT: L1={level_counts["level_1"]} L2={level_counts["level_2"]} L3={level_counts["level_3"]} UNK={level_counts["unknown"]} | Species tree source: {tree_source_display}</p>
{quality_note}
{quality_table}
{species_tree_table}
{conflict_table}
{bootstrap_note}
<div class="note">{hgt_availability_html}</div>
{risk_distribution_html}
{hgt_table}
<div class="note">{iqtree_note}</div>
{marker_table}
</body></html>"""


class PlainTextReportGenerator:
    def __init__(self, config: ReportConfig, force: bool = False, no_clobber: bool = False):
        self.config = config
        self.force = force
        self.no_clobber = no_clobber

    def generate(
        self,
        marker_result: MarkerSelectionResult,
        hgt_report: HGTReport,
        phylo_result: PhylogeneticResult,
        quality_data: Optional[QualityData],
        runtime_info: RuntimeInfo,
        input_dir: str = "",
        mode: str = "",
    ) -> PlainTextReportOutput:
        logger.info(f"{canonical_log_tag('4')} Generating plain text outputs")
        output_dir = Path(self.config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        reports_dir = output_dir / "Phase5_reports"
        trees_dir = output_dir / "Phase4_trees"
        alignments_dir = output_dir / "Phase4_alignments"
        reports_dir.mkdir(parents=True, exist_ok=True)
        trees_dir.mkdir(parents=True, exist_ok=True)
        alignments_dir.mkdir(parents=True, exist_ok=True)
        prefix = Path(self.config.output_prefix).name
        output_files: List[str] = []

        output_files.append(self._write_marker_summary(
            marker_result, reports_dir, prefix,
            consistency_grades=getattr(phylo_result, "consistency_grades", None)
            or {},
        ))
        output_files.append(self._write_hgt_evaluation(hgt_report, reports_dir, prefix))
        # Machine-readable per-marker decision cards.
        output_files.append(
            self._write_decision_cards(hgt_report, Path(self.config.output_dir), prefix)
        )

        if phylo_result.supermatrix and phylo_result.supermatrix.tree:
            path = str(trees_dir / f"{prefix}.species_tree_concat.newick")
            self._write_file(path, phylo_result.supermatrix.tree.to_newick())
            output_files.append(path)

        if phylo_result.supermatrix and phylo_result.supermatrix.partition:
            path = str(alignments_dir / f"{prefix}.partition.nex")
            if not (self.no_clobber and Path(path).exists()):
                phylo_result.supermatrix.partition.write_nexus(path)
            output_files.append(path)

        if phylo_result.coalescent and phylo_result.coalescent.species_tree:
            coalescent_source = getattr(phylo_result.coalescent, "species_tree_source", None) or "none"
            # 仅当 ASTRAL-III 真正产出合并物种树时才写入文件;
            # Source == 'none' 表示 ASTRAL-III 失败或样本不足,此时不输出空树文件.
            if coalescent_source == "astral":
                path = str(trees_dir / f"{prefix}.species_tree_{coalescent_source}.newick")
                self._write_file(path, phylo_result.coalescent.species_tree.to_newick())
                output_files.append(path)

        if phylo_result.coalescent and phylo_result.coalescent.gene_trees:
            path = str(trees_dir / f"{prefix}.gene_trees.newick")
            # 标准 multi-newick: 每行一棵基因树、无 FASTA 风格 `>` 头, 以便下游
            # 树工具(含 ASTRAL)直接解析. 基因顺序与 coalescent.gene_trees 一致.
            lines = []
            for gene_id, tree in phylo_result.coalescent.gene_trees.items():
                newick = tree.to_newick() if hasattr(tree, "to_newick") else str(tree)
                lines.append(newick + "\n")
            self._write_file(path, "".join(lines))
            output_files.append(path)

        output_files.append(self._write_pipeline_summary(
            marker_result, hgt_report, phylo_result, runtime_info, reports_dir, prefix,
            quality_source=getattr(quality_data, "quality_source", "checkm") if quality_data else "checkm",
            input_dir=input_dir,
            mode=mode,
        ))

        return PlainTextReportOutput(file_paths=output_files)

    def _write_file(self, path: str, content: str) -> str:
        if self.no_clobber and Path(path).exists():
            logger.info(f"  Skipping existing file: {path}")
            return path
        target = Path(path).resolve()
        target.write_text(content, encoding="utf-8", newline="\n")
        return path

    def _write_marker_summary(self, marker_result, output_dir, prefix,
                              consistency_grades=None) -> str:
        path = str(output_dir / f"{prefix}.marker_summary.tsv")
        # Pis + effective_columns appended on the right
        # (existing column order untouched).
        # Consistency_grade appended last; NA whenever the
        # Criterion did not run (the default risk mode), never a borrowed value.
        consistency_grades = consistency_grades or {}
        lines = [
            "marker_id\toccupancy_score\tmarker_quality_score\tmarker_quality_level"
            "\tpis\teffective_columns\tconsistency_grade\n"
        ]
        if marker_result.marker_set:
            for mid in marker_result.marker_set.markers:
                score = marker_result.quality_scores.get(mid)
                occ = marker_result.marker_set.occupancy_scores.get(mid, 0)
                grade = consistency_grades.get(mid, "NA")
                if score:
                    lines.append(
                        f"{mid}\t{occ:.4f}\t{_fmt(score.overall_score)}\t"
                        f"{score.level.value}\t{_plain(score.pis)}\t"
                        f"{_plain(score.effective_columns)}\t{grade}\n"
                    )
                else:
                    lines.append(f"{mid}\t{occ:.4f}\tNA\tNA\tNA\tNA\tNA\n")
        return self._write_file(path, "".join(lines))

    def _write_hgt_evaluation(self, hgt_report, output_dir, prefix) -> str:
        path = str(output_dir / f"{prefix}.hgt_evaluation.tsv")
        # Existing 5 columns stay first; appends the
        # Decision-card columns on the right.
        header = (
            "marker_id\toverall_risk\thgt_risk_level\thgt_evidence_confidence\tnotes\t"
            "detector\trank_used\tn_total\tn_mono\t"
            "rf\trf_state\tquartet\tquartet_state\tmonophyly\tmonophyly_state\t"
            "n_signals_used\tweights\trisk_basis\tthresholds_in_effect\tfar_active\t"
            "trimming_regime\tgene_tree_source\tupstream_shared\tassertion_ids_fired\n"
        )
        lines = [header]
        for ev in hgt_report.marker_evaluations:
            notes = ev.notes.replace("\t", " ").replace("\n", " ")
            card = getattr(ev, "decision_card", {}) or {}
            rf = card.get("rf", {}) or {}
            quartet = card.get("quartet", {}) or {}
            mono = card.get("monophyly", {}) or {}
            weights = card.get("weights", {}) or {}
            thresholds = card.get("thresholds_in_effect", {}) or {}
            row = [
                ev.marker_id,
                _fmt(ev.overall_risk) if ev.overall_risk is not None else "NA",
                ev.level.value,
                ev.confidence,
                notes,
                str(card.get("detector", "")),
                _plain(card.get("rank_used")),
                _plain(card.get("n_total")),
                _plain(card.get("n_mono")),
                _plain(rf.get("value")),
                str(rf.get("state", "")),
                _plain(quartet.get("value")),
                str(quartet.get("state", "")),
                _plain(mono.get("value")),
                str(mono.get("state", "")),
                _plain(card.get("n_signals_used")),
                json.dumps(weights, sort_keys=True) if weights else "NA",
                str(card.get("risk_basis", "")),
                json.dumps(thresholds, sort_keys=True) if thresholds else "NA",
                _plain(card.get("far_active")),
                str(card.get("trimming_regime", "")),
                str(card.get("gene_tree_source", "")),
                _plain(card.get("upstream_shared")),
                ",".join(card.get("assertion_ids_fired", []) or []) or "NA",
            ]
            lines.append("\t".join(str(c) for c in row) + "\n")
        return self._write_file(path, "".join(lines))

    def _write_decision_cards(self, hgt_report, output_dir, prefix) -> str:
        """Phase5_evidence/decision_<n>.json per marker."""
        evidence_dir = Path(output_dir) / "Phase5_evidence"
        evidence_dir.mkdir(parents=True, exist_ok=True)
        n = 0
        for ev in getattr(hgt_report, "marker_evaluations", []) or []:
            card = getattr(ev, "decision_card", None)
            if not card:
                continue
            target = (evidence_dir / f"decision_{n + 1:04d}.json").resolve()
            import json as _json

            target.write_text(
                _json.dumps(card, ensure_ascii=False, indent=2, default=str),
                encoding="utf-8", newline="\n",
            )
            n += 1
        logger.info(f"  Decision cards written: {n} -> {evidence_dir}")
        return str(evidence_dir)

    def _write_pipeline_summary(self, marker_result, hgt_report, phylo_result, runtime_info, output_dir, prefix, quality_source: str = "checkm", input_dir: str = "", mode: str = "") -> str:
        path = str(output_dir / f"{prefix}.pipeline_summary.txt")
        n_genomes = marker_result.occupancy_matrix.n_genomes if marker_result.occupancy_matrix else 0
        n_markers = len(marker_result.marker_set.markers) if marker_result.marker_set else 0
        strategy = marker_result.marker_set.strategy.value if marker_result.marker_set else "unknown"

        # 真实可用物种树的优先级: astral/consensus (真合并物种树) >
        # Concat (supermatrix) > first_gene_tree(退化替代) > none.
        coalescent_source = getattr(getattr(phylo_result, "coalescent", None),
                                   "species_tree_source", None)
        has_coalescent_tree = bool(
            coalescent_source and coalescent_source != "none"
            and has_usable_tree(getattr(getattr(phylo_result, "coalescent", None), "species_tree", None))
        )
        has_concat_tree = has_usable_tree(
            getattr(getattr(phylo_result, "supermatrix", None), "tree", None)
        )
        tree_source = prioritize_tree_source(
            coalescent_source if has_coalescent_tree else None,
            has_concat_tree,
        )

        # Level marker lists. "unknown" 必须存在: 筛查未运行时标记为
        # MarkerLevel.UNKNOWN, 缺桶会 KeyError 使报告阶段崩溃.
        level_markers: Dict[str, List[str]] = {"level_1": [], "level_2": [], "level_3": [], "unknown": []}
        for ev in hgt_report.marker_evaluations:
            level_markers.setdefault(ev.level.value, []).append(ev.marker_id)

        # Concat tree average support
        concat_support_str = "N/A"
        concat_tree = getattr(getattr(phylo_result, "supermatrix", None), "tree", None)
        if concat_tree is not None and hasattr(concat_tree, "get_average_support"):
            # ``None`` is NOT_MEASURABLE and used to be formatted with ":.2f",
            # I.e. the honest answer only appeared because the resulting
            # TypeError was swallowed by the bare except below. Judge it
            # Explicitly; "N/A" now means unmeasurable by construction.
            try:
                _concat_support = concat_tree.get_average_support()
            except Exception:  # Noqa: BLE001 - a report cell, never a crash
                _concat_support = None
            if _concat_support is not None:
                concat_support_str = f"{_concat_support:.2f}"

        # Coalescent gene tree count
        n_gene_trees = 0
        if phylo_result.coalescent and phylo_result.coalescent.gene_trees:
            n_gene_trees = len(phylo_result.coalescent.gene_trees)

        # RF distance if available
        rf_str = "N/A"
        quartet_agreement_str = "N/A"
        topology_note = ""
        if hasattr(phylo_result, "conflict_report") and phylo_result.conflict_report:
            rf = getattr(phylo_result.conflict_report, "normalized_rf", None)
            if rf is not None:
                rf_str = f"{rf:.4f}"
            qa = getattr(phylo_result.conflict_report, "quartet_agreement", None)
            if qa is not None:
                quartet_agreement_str = f"{qa:.4f}"
            topology_note = getattr(phylo_result.conflict_report, "topology_note", "") or ""

        # HGT step availability note
        hgt_step_notes: List[str] = []
        if not hgt_report.phylogenetic_enabled:
            hgt_step_notes.append("Phylogenetic step disabled.")
        else:
            hgt_step_notes.append("Phylogenetic HGT step was enabled.")
        hgt_step_availability = " ".join(hgt_step_notes)

        # Phylogenetic inference notes
        phylo_notes: List[str] = []
        if topology_note:
            phylo_notes.append(topology_note)
        # IQ-TREE composition chi2 test commonly fails on GC/codon-specialised clades.
        # This is a known, usually benign warning from the IQ-TREE log and does not
        # Prevent tree construction.
        phylo_notes.append(
            "IQ-TREE may report 'failed composition chi2 test' for GC/codon-specialised "
            "clades (e.g. halophilic archaea). This warning is typically benign and does "
            "not invalidate the inferred tree."
        )
        phylo_note_text = " ".join(phylo_notes)

        lines = [
            "=" * 60,
            "MarkerFinder Pipeline Summary",
            "=" * 60,
            "",
            f"Generated: {datetime.now().isoformat()}",
            f"Runtime: {runtime_info.duration:.1f}s",
            "",
        ]
        # A low evidence-coverage run must be flagged at the
        # Very top — grading conclusions on a mostly-unmeasured marker set is
        # Exactly the retracted-paper failure mode.
        _min_cov = float(getattr(self.config, "require_evidence_coverage", 0.5) or 0.5)
        _cov = float(getattr(hgt_report, "evidence_coverage", 0.0) or 0.0)
        if hgt_report.total_markers and _cov < _min_cov:
            lines.append(
                f"WARNING: EVIDENCE COVERAGE LOW — only {_cov:.0%} of markers were "
                f"truly measured (floor {_min_cov:.0%}); HGT-related conclusions "
                f"rest on a small measured subset."
            )
            lines.append("")
        # Coverage alone does not tell a reader WHAT is missing or WHY.
        # The run-level ledger is written by the call sites that could not
        # Measure, so this declaration cannot be added later by guesswork.
        from markerfinder.utils import unmeasured as _unmeasured

        _declaration = _unmeasured.render_lines()
        if _declaration:
            lines.extend(_declaration)
        # In consistency/hybrid mode the combined RF+quartet score must be
        # Demoted to "ranking only", and the reader has to be able to tell which
        # Criterion decided this run -- the same numbers look identical either way.
        _mode = str(getattr(self.config, "hgt_mode", "risk") or "risk")
        if _mode in ("consistency", "hybrid"):
            _grades = getattr(phylo_result, "consistency_grades", None) or {}
            if _grades and _mode == "hybrid":
                lines.append(
                    "Adjudication: HYBRID conjunction — a marker passes only if "
                    "it clears the risk level screen AND grades CONSISTENT "
                    "both must pass; per-marker verdicts are in "
                    f"the evaluation notes, {len(_grades)} marker(s) graded; "
                    "see consistency_grade column and excluded_profile.tsv). "
                    "The weighted RF+quartet risk score is used for RANKING "
                    "ONLY and does not decide a marker's status."
                )
            elif _grades:
                lines.append(
                    "Adjudication: cross-framework CONSISTENCY grade "
                    f"({len(_grades)} marker(s) graded; see consistency_grade "
                    "column and excluded_profile.tsv). The weighted RF+quartet "
                    "risk score is used for RANKING ONLY and does not decide a "
                    "marker's status."
                )
            else:
                lines.append(
                    f"WARNING: consistency demotion NOT IN EFFECT — {_mode} mode was "
                    "requested but no marker received a consistency grade "
                    "(missing trees or an illegal reference), so the run fell back "
                    "to risk grading. Read every level below as risk-based."
                )
            lines.append("")
        # How independent the two legs actually are
        # Must be READABLE in the product. It was computed and stored on
        # ConflictReport.independence but never surfaced anywhere — the same
        # "computed, stored, then dropped" channel defect as baseline.
        _indep = getattr(
            getattr(phylo_result, "conflict_report", None), "independence", None
        )
        if _indep is not None:
            lines.append(
                f"Leg independence: shared gene trees "
                f"{_indep.shared_gene_tree_ratio:.0%} (cache reuse across the two "
                f"legs); marker-set Jaccard {_indep.marker_set_jaccard:.2f}; "
                f"same trimming regime: "
                f"{'yes' if _indep.same_trimming_regime else 'no'}; reference "
                f"built from the markers under test: "
                f"{'yes' if _indep.ref_built_from_tested_markers else 'no'}"
            )
            if (_indep.shared_gene_tree_ratio >= 0.5
                    or _indep.ref_built_from_tested_markers):
                lines.append(
                    "NOTE: the two legs are NOT independent upstream — agreement "
                    "between them is weaker evidence than independent replication "
                    "(risk R2), so recommendation confidence is capped at medium "
                    "rather than reported as high."
                )
            lines.append("")
        # The dataset-level verdict itself, computed inside the
        # Run. The cap note above used to claim a property of a recommendation
        # Nothing ever showed -- the first product-level pass found the
        # Verdict object had zero consumers ( family), so the summary
        # Asserted a fact about thin air. It is now on the page, and the cap
        # Note and the verdict can be checked against each other.
        _reco = getattr(phylo_result, "tree_recommendation", None)
        if _reco is not None:
            _reco_conf = str(getattr(_reco, "confidence", "unknown") or "unknown")
            _reco_reason = str(getattr(_reco, "reason", "") or "")
            _reco_tree = getattr(_reco, "recommended_tree", None)
            lines.append(
                f"Tree recommendation: confidence={_reco_conf}; "
                f"recommended_tree={'yes' if _reco_tree is not None else 'none'}; "
                f"reason: {_reco_reason}"
            )
            if (_indep is not None
                    and (_indep.shared_gene_tree_ratio >= 0.5
                         or _indep.ref_built_from_tested_markers)
                    and _reco_conf == "high"):
                lines.append(
                    "WARNING: the independence note above promises a medium cap, "
                    "but the recommendation object says 'high' — the two "
                    "disagree, which is a wiring defect, not a judgement call."
                )
            lines.append("")
        # Far-distance relaxation is declared at report top,
        # Not only as a per-marker note.
        if getattr(hgt_report, "far_active", False):
            lines.append(
                "WARNING: FAR-DISTANCE DATASET — this run triggered adaptive far-distance "
                f"mode: the Level 2/3 boundary was relaxed from "
                f"{self.config.level2_max:.2f} to "
                f"{float(getattr(self.config, 'level2_max_far', 0.95)):.2f}; every affected "
                f"marker's notes record 'far-distance'."
            )
            lines.append("")
        lines += [
            "Input/Output:",
            f"  Input directory:  {input_dir or 'N/A'}",
            f"  Output directory: {self.config.output_dir}",
            f"  Analysis mode:    {mode or 'N/A'}",
            f"  Selection strategy: {strategy}",
            f"  Intermediate files saved: {'yes (see Phase4_intermediate/)' if self.config.save_intermediates else 'no'}",
            "",
            "Dataset:",
            f"  Genomes: {n_genomes}",
            f"  Markers selected: {n_markers}",
            f"  Quality source: {quality_source}",
            "",
            "HGT screening:",
            # Evidence coverage leads the HGT section; below
            # The configured floor it is a prominent warning, not a footnote.
            f"  Evidence coverage: {getattr(hgt_report, 'evidence_coverage', 0.0):.2f} "
            f"({getattr(hgt_report, 'n_actually_measured', 0)}/{hgt_report.total_markers} markers truly measured)",
            f"  Level 1 (clean):     {hgt_report.level1_count}",
            f"  Level 2 (suspicious): {hgt_report.level2_count}",
            f"  Level 3 (excluded):  {hgt_report.level3_count}",
            f"  Unknown (screen skipped, kept): {getattr(hgt_report, 'unknown_count', sum(1 for e in hgt_report.marker_evaluations if e.level.value == 'unknown'))}",
            f"  Mean HGT risk: {hgt_report.mean_risk:.4f}",
            "",
            "HGT step availability:",
            f"  {hgt_step_availability}",
            "",
            "Level 1 markers:",
            f"  {', '.join(level_markers['level_1']) if level_markers['level_1'] else 'None'}",
            "",
            "Level 2 markers:",
            f"  {', '.join(level_markers['level_2']) if level_markers['level_2'] else 'None'}",
            "",
            "Level 3 markers (excluded from tree inference):",
            f"  {', '.join(level_markers['level_3']) if level_markers['level_3'] else 'None'}",
            "",
            "Unknown markers (HGT screen skipped; kept for tree inference):",
            f"  {', '.join(level_markers['unknown']) if level_markers['unknown'] else 'None'}",
            "",
            "Phylogenetic inference:",
            f"  Primary species tree source: {tree_source}",
            f"  Concatenation tree avg. UFBOOT support: {concat_support_str}",
            f"  Coalescent gene trees used: {n_gene_trees}",
            f"  Concat vs coalescent normalized RF distance: {rf_str}",
            f"  Concat vs coalescent quartet agreement: {quartet_agreement_str}",
            "",
            "Phylogenetic inference notes:",
            f"  {phylo_note_text}",
            "",
            "HGT evidence confidence key:",
            "  high = phylogenetic step produced a clean signal",
            "  low  = phylogenetic step could not be run or produced no signal",
            "",
        ]
        content = "\n".join(lines)
        return self._write_file(path, content)


class ReportGeneratorModule:
    def __init__(self, config: ReportConfig, force: bool = False, no_clobber: bool = False):
        self.config = config
        self.force = force
        self.no_clobber = no_clobber
        self.html_generator = InteractiveReportGenerator(config, force=force, no_clobber=no_clobber)
        self.text_generator = PlainTextReportGenerator(config, force=force, no_clobber=no_clobber)

    def run(
        self,
        marker_result: MarkerSelectionResult,
        hgt_report: HGTReport,
        phylo_result: PhylogeneticResult,
        quality_data: Optional[QualityData],
        runtime_info: RuntimeInfo,
        input_dir: str = "",
        mode: str = "",
    ) -> tuple:
        html_output = self.html_generator.generate(
            marker_result, hgt_report, phylo_result, quality_data, runtime_info
        )
        text_output = self.text_generator.generate(
            marker_result, hgt_report, phylo_result, quality_data, runtime_info,
            input_dir=input_dir, mode=mode,
        )
        return html_output, text_output
