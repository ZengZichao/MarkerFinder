"""Every configuration knob must be consumed, or be declared and documented.

The single most common defect this shape of suite finds is a switch that exists but is
never read: ``--scan-stability-min``, ``--min-informative-sites``,
``enable_composition_screen``, ``ConflictDetector.quartet_threshold``, and the
baseline's own ("存为属性后从未被读取"). A dataclass field with no reader is
worse than a missing field, because a reader of ``config.example.yaml`` or the
MANUAL reasonably believes setting it changes something.

This is a ratchet, not a report: a new unread config field fails the suite, and an
"intentionally reserved" field only passes if it is listed here WITH a reason AND
named in the MANUAL's not-wired section -- so the excuse has to be written down in
the place users actually read.
"""

from __future__ import annotations

import ast
from pathlib import Path

from markerfinder import config as config_module

PACKAGE = Path(__file__).resolve().parents[2] / "markerfinder"
MANUAL = PACKAGE.parent / "MANUAL.CN.md"

# Fields that are deliberately NOT consumed by the pipeline, and why. Anything
# Here must also appear by name in MANUAL.CN.md's "尚未接线" section.
RESERVED: dict[str, str] = {
    "verify_db": "数据库哈希校验目前只在 `--check` 自检里做，run 路径不消费此开关",
    "databases": "同上：`--check` 自行发现 db/ 下的库，run 路径不读该映射",
    "min_completeness": "完整度门槛走 Phase 0 自适应占居率路径，MAGConfig.min_completeness 未被查询",
    "use_diamond": "DIAMOND 直系同源解析未接线（Phase 1.5 明确记日志跳过）",
    "graph_clustering": "同上：图聚类直系同源解析未接线",
    "min_allele_freq": "仅保留给未接入主流水线的 MAGHeterogeneityHandler（其 docstring 自陈）",
    "consensus_threshold": "仅保留给未接入主流水线的 MAGHeterogeneityHandler；"
                           "主流程每标记只取最优命中，不做共识序列重建",
    "minor_allele_threshold": "仅保留给未接入主流水线的 MAGHeterogeneityHandler；"
                              " minority-allele 判定在流水线中无调用点",
    "identity_threshold": "仅保留给未接入主流水线的 MAGHeterogeneityHandler；"
                          "等位身份阈值在流水线中无调用点",
}

CONFIG_CLASSES = sorted(
    name
    for name, value in vars(config_module).items()
    if isinstance(value, type) and name.endswith("Config")
)


def _fields_by_class() -> dict:
    tree = ast.parse((PACKAGE / "config.py").read_text(encoding="utf-8"))
    out = {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name.endswith("Config"):
            out[node.name] = [
                stmt.target.id
                for stmt in node.body
                if isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name)
            ]
    return out


def _consumed_names() -> set:
    """Attribute reads, keyword arguments and string tokens across the package.

    Deliberately generous: a field reached through ``getattr(cfg, "x")``, a dict
    key in a serialised product, or a template placeholder all count as consumed.
    The point is to catch fields that are named NOWHERE outside their own
    declaration -- those cannot possibly do anything.
    """
    names = set()
    for path in sorted(PACKAGE.rglob("*.py")):
        if path.name == "config.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute):
                names.add(node.attr)
            elif isinstance(node, ast.keyword) and node.arg:
                names.add(node.arg)
            elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                for token in node.value.replace('"', " ").replace("'", " ").split():
                    names.add(token.strip(".,:;()[]{}"))
                names.add(node.value)
    return names


def test_control_the_scanner_fires_on_a_planted_unread_knob():
    """Without this, a scanner that always returns "consumed" makes the ratchet
    vacuously green -- the exact failure mode this file exists to prevent."""
    fields = {"SomeConfig": ["definitely_not_used_anywhere_zzz", "output_dir"]}
    consumed = _consumed_names()
    unread = [f for f in fields["SomeConfig"] if f not in consumed]
    assert unread == ["definitely_not_used_anywhere_zzz"], unread


def test_config_classes_are_discovered():
    assert CONFIG_CLASSES, "no *Config classes found -- the scan is broken"
    for expected in ("PipelineConfig", "HGTConfig", "ReportConfig", "HeterogeneityConfig"):
        assert expected in CONFIG_CLASSES, expected


def test_every_config_field_is_read_or_declared_reserved_and_documented():
    consumed = _consumed_names()
    undocumented = []
    unreserved = []
    manual_text = MANUAL.read_text(encoding="utf-8")
    for cls, fields in sorted(_fields_by_class().items()):
        for field in fields:
            if field in consumed:
                continue
            if field not in RESERVED:
                unreserved.append(f"{cls}.{field}")
                continue
            if field not in manual_text:
                undocumented.append(f"{cls}.{field}")
    assert not unreserved, (
        f"unread configuration knobs (no attribute read, no keyword use, no "
        f"string reference anywhere): {unreserved}. Wire them into the decision "
        f"they claim to control, or add them to RESERVED with a reason AND name "
        f"them in MANUAL.CN.md's 尚未接线 section."
    )
    assert not undocumented, (
        f"reserved knobs missing from MANUAL.CN.md: {undocumented}"
    )


def test_reserved_reasons_are_not_placeholders():
    for field, reason in RESERVED.items():
        assert len(reason) >= 12, (field, reason)


def test_control_reserved_list_matches_reality():
    """A reserved entry whose field has since been wired must be removed, or the
    allow-list quietly becomes a hiding place."""
    consumed = _consumed_names()
    stale = [field for field in RESERVED if field in consumed]
    assert not stale, (
        f"these knobs are now consumed and should leave RESERVED: {stale}"
    )

# ── the same rule for the data model ─────────────────────────────────────
#
# Only classes the production code actually constructs: a class nobody builds is
# Dead code (a different clean-up), while a built-and-returned object whose field
# Nobody reads is the pattern -- computed, stored, then dropped. The scan
# Below distinguishes them, and the verdict on each reserved field records which
# Of the two it is.

MODEL_RESERVED = {
    # Branch-level conflict classification (HGT vs ILS vs method bias) is
    # Explicitly NOT implemented; ConflictDetector's own docstring says the
    # Fields are reserved placeholders kept for API compatibility (
    # Removed the enum members that had no story left to tell).
    "n_quartet_conflicts": "分支级冲突分型未实现：ConflictDetector 类文档自陈这些是保留占位",
    "quartet_conflicts": "同上：分支级冲突分型的保留列表，无生产者也无消费者",
    "n_conflicting_branches": "同上：分支级冲突分型的保留计数",
    "conflicting_branches": "同上：分支级冲突分型的保留列表",
    "gene_tree_agreement": "保留字段：逐支基因树一致率属未实现的分支级分析",
    "conflict_summary": "保留字段：分支级冲突摘要未实现，报告层从不读取",
    # Declared, never filled by any producer in the package (verified by a
    # Write-scan): advertising them as measured numbers would be worse, so they
    # Stay listed here as known-unused surface rather than silently drifting.
    "taxonomic_group": "无生产者：分类信息走 taxonomy_map 路径，不回填到 Genome 上",
    "gc_content": "无生产者：GC 由 composition.py 按标记计算（同名函数，不是本字段），"
                  "Genome/GenomeQuality 上的这两处从未被写入或读取",
    "gc_std": "无生产者：GC 诊断由 composition.py 按标记计算，不在 Genome 上",
    "strain_heterogeneity": "无生产者：CheckM 预计算表只提供完整度/污染度",
    "n_contigs": "无生产者：CheckM 预计算表与 .faa 输入都不提供 contig 计数",
    "genome_size": "无生产者：未统计序列总长，只有完整度/污染度参与判定",
    "coding_density": "无生产者：编码密度需要 CDS 注释，本流程只消费蛋白 .faa",
    "allele_frequencies": "无生产者：直系同源解析（DIAMOND/图聚类）未接线，Phase 1.5 明确跳过",
    "effective_information": "无生产者：矩阵优化仅记录剪除计数",
    "weight": "无生产者：分区权重由 IQ-TREE 自行决定，MarkerFinder 不写权重",
    "quartet_support": "无生产者也无消费者：本轮已如实注明为保留字段（models/tree.py）",
    "flagged_branches": "无生产者：分支级标注属未实现的分支级分析",
    "coalescent_result": "无生产者：PhaseContext 只承载跨阶段的树与状态对象",
    "is_ingroup": "无生产者：群内/群外判定由分类学层级路径完成，不回填 BlastHit",
    "taxonomic_distance": "无生产者：远缘度由分类层级路径与 hgt_config.far 判定算，不回填命中对象",
}


def _model_classes_and_fields() -> dict:
    out = {}
    for path in sorted((PACKAGE / "models").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                fields = [
                    stmt.target.id
                    for stmt in node.body
                    if isinstance(stmt, ast.AnnAssign)
                    and isinstance(stmt.target, ast.Name)
                ]
                if fields:
                    out[node.name] = fields
    return out


def _production_constructed() -> set:
    built = set()
    for path in sorted(PACKAGE.rglob("*.py")):
        if path.parent.name == "models":
            continue
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                built.add(node.func.id)
    return built


def _names_used_everywhere() -> set:
    """Names the PRODUCTION code refers to: attribute reads, keywords, strings.

    Two exclusions, both learned the hard way:
    * this file -- its reserved dicts name every reserved entry, so scanning it
      would make each key "used", the primary test vacuous and the staleness test
      permanently red. Self-reference is the kind of green that means nothing.
    * tests/ -- a test that mentions a field is not a pipeline that consumes it,
      which is the exact distinction this ratchet exists to enforce.
    """
    here = Path(__file__).resolve()
    names = set()
    for path in sorted(PACKAGE.rglob("*.py")):
        if path.resolve() == here:
            continue
        if path.name in {"config.py", "phases.py"}:
            # Config.py holds the declarations themselves; phases.py is a
            # Registry whose entire job is to name other things.
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute):
                names.add(node.attr)
            elif isinstance(node, ast.keyword) and node.arg:
                names.add(node.arg)
            elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                names.add(node.value)
                for token in node.value.replace('"', " ").split():
                    names.add(token.strip(".,:;()[]{}"))
    return names


def test_control_the_model_scanner_finds_live_classes_and_skips_dead_ones():
    fields = _model_classes_and_fields()
    built = _production_constructed()
    assert "HGTEvaluation" in fields and "HGTEvaluation" in built
    # A model class nothing in the pipeline builds must not be counted as live.
    dead = [name for name in fields if name not in built]
    assert dead, "no dead model classes at all -- the scan is probably broken"


def test_live_model_fields_are_read_or_declared_reserved():
    fields = _model_classes_and_fields()
    built = _production_constructed()
    used = _names_used_everywhere()
    unaccounted = []
    for cls, names in sorted(fields.items()):
        if cls not in built:
            continue
        for field in names:
            if field in used or field in MODEL_RESERVED:
                continue
            unaccounted.append(f"{cls}.{field}")
    assert not unaccounted, (
        f"fields on production-built model objects that nothing reads and "
        f"nothing declares: {unaccounted}. Either consume them, or add them to "
        f"MODEL_RESERVED with the reason (which producer is missing, or which "
        f"unimplemented feature they were reserved for)."
    )


def test_reserved_model_fields_must_be_stale_checked_and_explained():
    used = _names_used_everywhere()
    stale = sorted(name for name in MODEL_RESERVED if name in used)
    assert not stale, (
        f"these model fields are now consumed and must leave MODEL_RESERVED: "
        f"{stale}"
    )
    for name, reason in MODEL_RESERVED.items():
        assert len(reason) >= 12, (name, reason)
