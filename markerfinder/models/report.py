from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class RuntimeInfo:
    start_time: float = 0.0
    end_time: float = 0.0
    duration: float = 0.0


@dataclass
class DownloadLink:
    path: str = ""
    label: str = ""
    format: str = ""


@dataclass
class ReportData:
    summary: str = ""
    coverage_heatmap: str = ""
    quality_dashboard: str = ""
    phylo_tree: Optional[object] = None
    hgt_panel: str = ""
    conflict_panel: str = ""
    download_links: List[DownloadLink] = field(default_factory=list)


@dataclass
class ReportOutput:
    html_path: str = ""
    data_files: List[str] = field(default_factory=list)


@dataclass
class PlainTextReportOutput:
    file_paths: List[str] = field(default_factory=list)


@dataclass
class PhaseContext:
    quality_data: Optional[object] = None
    adaptive_params: Optional[object] = None
    marker_set: Optional[object] = None
    hgt_evaluations: Optional[List] = None
    coalescent_result: Optional[object] = None
    species_tree: Optional[object] = None


@dataclass
class PipelineResult:
    preprocessing: Optional[object] = None
    marker_selection: Optional[object] = None
    hgt_report: Optional[object] = None
    phylogenetic: Optional[object] = None
    report: Optional[ReportOutput] = None
    text_report: Optional[PlainTextReportOutput] = None
    runtime: RuntimeInfo = field(default_factory=RuntimeInfo)
    parameters: Optional[Dict] = None
    # 可用性信号: 在 __main__.py 汇总检验中用于判定管线是否"真正成功".
    usable_marker_count: int = 0  # 通过 HGT 过滤(等级 1/2)可用于建树的标记数
    species_tree_source: str = "none"  # 物种树实际来源 (astral/concat/consensus/first_gene_tree/none)
