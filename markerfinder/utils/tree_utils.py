"""Tree manipulation utilities."""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from markerfinder.exceptions import PhyloToolUnavailable
from markerfinder.models.tree import Tree
from markerfinder.utils.etree import TreeMeasureError

logger = logging.getLogger(__name__)


def strip_nhx_annotations(newick: str) -> str:
    """剥离 Newick 字符串中的 NHX / extended-newick 注释（防御性）。

    ASTRAL-III (>=5.7) 在 ``-t 2`` 模式下会在内部节点附加形如
        '[q1=0.97;q2=0.02;pp1=0.999;QC=3;EN=30]'
    的 quartet 支持度 / posterior 注释. MarkerFinder 已改用 ``-t 1`` 运行
    ASTRAL (LPP 作为常规 Newick 节点标签, 直接保留为分支支持度), 正常产出
    不再含 NHX 注释. 本函数作为防御性后路: 若上游被改为 ``-t 2`` 或读入外部
    含 NHX 注释的树, 仍能把这些方括号注释剥去, 避免 ete3 抛
    'Unexpected newick format' 解析失败.

    处理两步:
      1. 剥除单引号包裹的注释块: '...'(含内部方括号内容)整体移除.
      2. 再剥除残余的裸 '[...]' 方括号片段.
    两步剥除后,可能留下紧接的空单引号 (形如 Name'':length) 一并清理.
    """
    if not newick:
        return newick
    # 第 1 步: 剥除 '...' 块 (ASTRAL 注释以单引号包裹;非贪婪,注释内无单引号).
    stripped = re.sub(r"'[^']*'", "", newick)
    # 第 2 步: 剥除残余的裸 '[...]' 方括号片段.
    stripped = re.sub(r"\[.*?\]", "", stripped)
    return stripped.strip()


def parse_newick_tips(newick: str) -> List[str]:
    """从 Newick 字符串中提取所有叶节点名称"""
    inner = newick.strip().rstrip(";")
    tokens = re.findall(r"[(),]([^(),:]+)", inner)
    tips: List[str] = []
    seen: Set[str] = set()
    for t in tokens:
        name = t.strip()
        if name and name not in seen:
            tips.append(name)
            seen.add(name)
    return tips


def calculate_rf_distance(tree1_newick: str, tree2_newick: str) -> Tuple[Optional[int], Optional[float]]:
    """计算两棵树之间的 Robinson-Foulds 距离

    测量顺序：先走 ete3；ete3 不可导入或计算失败时，
    回退到：mod:`markerfinder.utils.etree` 的纯 Python split-set 真实测量
    （≤64 共享尖端）；两路都失败才返回 ``(None, None)`` 并在日志写明原因
    （不得静默降级）。

    Returns:
        ``(raw_rf, normalized_rf)``. On failure (empty input after stripping
        annotations, both measurement paths unavailable, or shared-tip cap
        exceeded) returns ``(None, None)``
        — distinct from a real result of ``(0, 0.0)`` (identical trees), so
        callers must treat ``None`` as "uncomputable / not applicable", never as
        a zero-distance match.
    """
    t1_newick = strip_nhx_annotations(tree1_newick)
    t2_newick = strip_nhx_annotations(tree2_newick)
    if not t1_newick or not t2_newick:
        logger.warning("RF distance: empty newick after stripping annotations")
        return None, None

    try:
        # The sanctioned ete3 entry point, so "dependency missing" has
        # Exactly one observable failure type instead of per-site guessing.
        from markerfinder.utils.etree import require_ete3

        ete3 = require_ete3()
        t1 = ete3.Tree(t1_newick)
        t2 = ete3.Tree(t2_newick)
        result = t1.robinson_foulds(t2, unrooted_trees=True)
        rf = result[0]
        max_rf = result[1]
        norm_rf = rf / max_rf if max_rf > 0 else 0.0
        return rf, norm_rf
    except PhyloToolUnavailable:
        logger.info(
            "ete3 unavailable — RF distance falls back to pure-Python "
            "split-set measurement (method='splits-python')"
        )
    except Exception as e:
        logger.warning(
            f"ete3 RF computation failed ({e}) — falling back to "
            "pure-Python split-set measurement (method='splits-python')"
        )

    try:
        from markerfinder.utils import etree as _etree

        return _etree.rf_distance(t1_newick, t2_newick)
    except TreeMeasureError as e:
        logger.warning(f"RF distance not measurable (split-set fallback): {e}")
        _declare_rf_unmeasured(str(e))
        return None, None
    except Exception as e:
        logger.warning(f"RF distance computation failed: {e}")
        _declare_rf_unmeasured(str(e))
        return None, None


def _declare_rf_unmeasured(cause: str) -> None:
    """RF distance unmeasurable must be declared, not just logged.

    Called on the path where BOTH ete3 and the pure-Python split-set measurement
    failed, i.e. the number really is absent for this pair of trees.
    """
    from markerfinder.utils import unmeasured

    unmeasured.record("rf_distance (Robinson-Foulds)", cause)


def get_quartet_topology(newick: str, quartet: Tuple[str, str, str, str]) -> Optional[str]:
    """获取 4-物种组合在树上的拓扑（规范化形式，对书写顺序不变）

    测量顺序与：func:`calculate_rf_distance` 相同：先 ete3，失败回退纯
    Python split-set。两路输出统一为规范化串（对内配对、配对间均排序），
    保证两种环境的度量可比（ 契约）。
    """
    try:
        # Single sanctioned ete3 entry point.
        from markerfinder.utils.etree import require_ete3

        ete3 = require_ete3()
        t = ete3.Tree(newick)
        t.prune(list(quartet), preserve_branch_length=True)
        # 规范化 ete3 输出（write(format=9) 依赖存储顺序，非旋转不变）。
        from markerfinder.utils import etree as _etree

        canonical = _etree.canonical_small_newick(t.write(format=9))
        return canonical if canonical is not None else t.write(format=9)
    except PhyloToolUnavailable:
        logger.info(
            "ete3 unavailable — quartet topology falls back to pure-Python "
            "split-set measurement (method='splits-python')"
        )
    except Exception as e:
        logger.warning(
            f"ete3 quartet topology failed ({e}) — falling back to "
            "pure-Python split-set measurement (method='splits-python')"
        )

    from markerfinder.utils import etree as _etree

    topology = _etree.quartet_topology(newick, quartet)
    if topology is None:
        from markerfinder.utils import unmeasured

        unmeasured.record(
            "quartet_topology",
            "neither ete3 nor the pure-Python split-set path could measure it "
            "(tip cap exceeded or unparseable Newick)",
        )
    return topology


def merge_newick_files(input_files: List[str], output_file: str) -> str:
    """合并多个 Newick 文件为单文件（每行一棵树）"""
    lines: List[str] = []
    for fpath in input_files:
        with open(fpath, encoding="utf-8", newline="") as f:
            for line in f:
                line = line.strip()
                if line:
                    lines.append(line + "\n")
    target = Path(output_file).resolve()
    target.write_text("".join(lines), encoding="utf-8", newline="\n")
    return output_file
