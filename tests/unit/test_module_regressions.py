"""Regression tests for the fixes documented in ``MarkerFinder-代码审查.md``.

Each test maps to a finding id (F1-F21) from the review report. The module is
ete3-import-safe: ete3 is only imported lazily inside the library (and even
then every use has a non-ete3 fallback), so these tests also run in
environments where ete3 is unavailable.
"""

from __future__ import annotations

import json
import types
from pathlib import Path
from unittest.mock import patch

import pytest

from markerfinder.cli.constants import EXIT_ARG_ERROR, EXIT_DATA_ERROR
from markerfinder.config import (
    HGTConfig,
    PipelineConfig,
    PhylogeneticConfig,
    ReportConfig,
    SelectionConfig,
)
from markerfinder.models.marker import MarkerLevel, SelectedMarkerSet
from markerfinder.models.pipeline_types import (
    HGTEvaluation,
    HGTReport,
    MarkerSelectionResult,
    PhylogeneticResult,
)
from markerfinder.models.report import RuntimeInfo
from markerfinder.models.tree import Tree


# ---------------------------------------------------------------------------
# F1 — report generators must accept MarkerLevel.UNKNOWN
# ---------------------------------------------------------------------------

def _unknown_hgt_report() -> HGTReport:
    return HGTReport(
        total_markers=2,
        marker_evaluations=[
            HGTEvaluation(marker_id="mk1", overall_risk=0.0,
                          level=MarkerLevel.UNKNOWN, confidence="unknown",
                          notes="screen skipped"),
            HGTEvaluation(marker_id="mk2", overall_risk=0.1,
                          level=MarkerLevel.LEVEL_1, confidence="high"),
        ],
        unknown_count=1,
    )


def _marker_result() -> MarkerSelectionResult:
    return MarkerSelectionResult(
        marker_set=SelectedMarkerSet(markers=["mk1", "mk2"],
                                     occupancy_scores={"mk1": 0.9, "mk2": 0.8}),
        quality_scores={},
    )


def test_f1_plain_text_report_handles_unknown(tmp_path):
    from markerfinder.modules.report_generator import PlainTextReportGenerator

    out = PlainTextReportGenerator(ReportConfig(output_dir=str(tmp_path))).generate(
        _marker_result(), _unknown_hgt_report(), PhylogeneticResult(), None,
        RuntimeInfo(duration=1.0),
    )
    txt = Path(out.file_paths[-1]).read_text(encoding="utf-8")
    assert "Unknown" in txt
    assert "mk1" in txt


def test_f1_html_report_handles_unknown(tmp_path):
    from markerfinder.modules.report_generator import InteractiveReportGenerator

    out = InteractiveReportGenerator(ReportConfig(output_dir=str(tmp_path))).generate(
        _marker_result(), _unknown_hgt_report(), PhylogeneticResult(), None,
        RuntimeInfo(duration=1.0),
    )
    html = Path(out.html_path).read_text(encoding="utf-8")
    assert "UNK=1" in html
    assert "unknown" in html


def test_f1_unknown_markers_kept_for_inference(tmp_path):
    """审查修复 F1 前置链: 无参照时 HGT 过滤给出 UNKNOWN 且标记全部保留."""
    from markerfinder.modules.hgt_filter import HGTFilterModule

    mod = HGTFilterModule(HGTConfig(), tmp_dir=str(tmp_path))
    report, levels, _far = mod.run(
        ["mk1", "mk2"],
        marker_sequences={"mk1": [{"id": "a", "seq": "M"}],
                          "mk2": [{"id": "a", "seq": "M"}]},
    )
    assert all(lv == MarkerLevel.UNKNOWN for lv in levels.values())
    assert report.unknown_count == 2
    excluded = {e.marker_id for e in report.marker_evaluations
                if e.level == MarkerLevel.LEVEL_3}
    assert not excluded


# ---------------------------------------------------------------------------
# F2 — StateSchemaError must degrade gracefully, not raise NameError
# ---------------------------------------------------------------------------

def test_f2_save_state_degrades_on_schema_error(tmp_path, monkeypatch):
    import markerfinder.pipeline as pipeline_mod
    from markerfinder.utils.state_codec import StateSchemaError

    pipe = pipeline_mod.MarkerFinderPipeline.__new__(pipeline_mod.MarkerFinderPipeline)
    pipe.config = PipelineConfig(output_dir=str(tmp_path))

    def _boom(output_dir, state):
        raise StateSchemaError("simulated encode failure")

    monkeypatch.setattr(pipeline_mod, "save_pipeline_state", _boom)
    pipe._save_state(completed_steps=["scan"])  # Must not raise NameError


def test_f2_load_state_degrades_on_schema_error(tmp_path, monkeypatch):
    import markerfinder.pipeline as pipeline_mod
    from markerfinder.utils.state_codec import StateSchemaError

    pipe = pipeline_mod.MarkerFinderPipeline.__new__(pipeline_mod.MarkerFinderPipeline)
    pipe.config = PipelineConfig(output_dir=str(tmp_path))

    def _boom(output_dir):
        raise StateSchemaError("simulated corrupt state")

    monkeypatch.setattr(pipeline_mod, "load_pipeline_state", _boom)
    assert pipe._load_state() == {}


# ---------------------------------------------------------------------------
# F3 — run_config.json replay flattens the nested "parameters" mapping
# ---------------------------------------------------------------------------

def _write_run_config(tmp_path: Path, params: dict) -> Path:
    payload = {
        "markerfinder_version": "0.1.0",
        "timestamp": "2026-01-01T00:00:00",
        "run_duration_seconds": 12.5,
        "parameters": params,
        "database_versions": {},
    }
    p = tmp_path / "run_config.json"
    p.write_text(json.dumps(payload), encoding="utf-8")
    return p


def test_f3_load_config_file_flattens_parameters(tmp_path):
    from markerfinder.config_loader import load_config_file

    p = _write_run_config(tmp_path, {"mode": "mag_adaptive", "threads": 4})
    data = load_config_file(str(p))
    assert data["mode"] == "mag_adaptive"
    assert data["threads"] == 4
    assert "parameters" not in data


def test_f3_apply_config_file_replays_run_config(tmp_path):
    from markerfinder.cli.parser import _apply_config_file, _build_parser

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "g.faa").write_text(">a\nMKV\n", encoding="utf-8")
    cfg = _write_run_config(tmp_path, {"mode": "mag_adaptive", "threads": 4})

    parser = _build_parser()
    args = parser.parse_args(["--config", str(cfg), "-i", str(genome_dir),
                              "-o", str(tmp_path / "out")])
    args = _apply_config_file(args, parser)
    assert args.mode == "mag_adaptive"
    assert args.threads == 4


def test_f3_apply_config_file_warns_on_unknown_keys(tmp_path, caplog):
    from markerfinder.cli.parser import _apply_config_file, _build_parser

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "g.faa").write_text(">a\nMKV\n", encoding="utf-8")
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("not_a_real_key: 1\nmode: expanded\n", encoding="utf-8")

    parser = _build_parser()
    args = parser.parse_args(["--config", str(cfg), "-i", str(genome_dir),
                              "-o", str(tmp_path / "out")])
    args = _apply_config_file(args, parser)
    assert any("not_a_real_key" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# F4 — --multi-tree-mode first/last/random actually select a tree
# ---------------------------------------------------------------------------

def test_f4_split_newick_trees_basic_and_quoted():
    from markerfinder.validation import split_newick_trees

    trees = split_newick_trees("((A,B),(C,D));\n((E,F),(G,H));\n")
    assert len(trees) == 2
    quoted = "((A,B),'weird;name');((C,D),(E,F));"
    assert len(split_newick_trees(quoted)) == 2


def test_f4_load_tree_tree_index_selectsRequested_tree(tmp_path):
    from markerfinder.validation import load_tree

    f = tmp_path / "multi.nwk"
    f.write_text("((A,B),(C,D));\n((E,F),(G,H));\n", encoding="utf-8")
    first = load_tree(str(f), validate=True, tree_index=0)
    last = load_tree(str(f), validate=True, tree_index=1)
    assert set(first.get_tips()) == {"A", "B", "C", "D"}
    assert set(last.get_tips()) == {"E", "F", "G", "H"}


def test_f4_load_tree_tree_index_out_of_range(tmp_path):
    from markerfinder.exceptions import TreeValidationError
    from markerfinder.validation import load_tree

    f = tmp_path / "multi.nwk"
    f.write_text("((A,B),(C,D));\n", encoding="utf-8")
    with pytest.raises(TreeValidationError, match="out of range"):
        load_tree(str(f), validate=True, tree_index=5)


def test_f4_handle_tree_input_last_mode(tmp_path):
    import markerfinder.cli.validation as cli_val

    f = tmp_path / "multi.nwk"
    f.write_text("((A,B),(C,D));\n((E,F),(G,H));\n", encoding="utf-8")
    args = types.SimpleNamespace(
        species_tree=str(f), tree=None, multi_tree_mode="last",
        strip_annotations=False,
    )
    tree = cli_val._handle_tree_input(args)
    assert set(tree.get_tips()) == {"E", "F", "G", "H"}


def test_f4_handle_tree_input_first_mode(tmp_path):
    import markerfinder.cli.validation as cli_val

    f = tmp_path / "multi.nwk"
    f.write_text("((A,B),(C,D));\n((E,F),(G,H));\n", encoding="utf-8")
    args = types.SimpleNamespace(
        species_tree=str(f), tree=None, multi_tree_mode="first",
        strip_annotations=False,
    )
    tree = cli_val._handle_tree_input(args)
    assert set(tree.get_tips()) == {"A", "B", "C", "D"}


def test_f4_handle_tree_input_split_raises(tmp_path):
    from markerfinder.cli.validation import _handle_tree_input
    from markerfinder.exceptions import MultiTreeError

    f = tmp_path / "multi.nwk"
    f.write_text("((A,B),(C,D));\n((E,F),(G,H));\n", encoding="utf-8")
    args = types.SimpleNamespace(
        species_tree=str(f), tree=None, multi_tree_mode="split",
        strip_annotations=False,
    )
    with pytest.raises(MultiTreeError, match="split is not supported"):
        _handle_tree_input(args)


# ---------------------------------------------------------------------------
# F5 — placeholder trees must not count as usable species trees
# ---------------------------------------------------------------------------

def test_f5_supermatrix_no_alignments_returns_tree_none(tmp_path):
    from markerfinder.modules.phylogenetic_inference import SupermatrixInference

    sm = SupermatrixInference(PhylogeneticConfig(tmp_dir=str(tmp_path)))
    result = sm.run({"C1": [], "C2": []}, [])
    assert result.tree is None
    assert result.partition is not None


def test_f5_coalescent_no_gene_trees_returns_species_tree_none(tmp_path):
    from markerfinder.modules.phylogenetic_inference import CoalescentInference

    ci = CoalescentInference(PhylogeneticConfig(tmp_dir=str(tmp_path)))
    result = ci.run({"C1": []}, [])
    assert result.species_tree is None
    assert result.species_tree_source == "none"


def test_f5_has_usable_tree_rejects_placeholder():
    from markerfinder.utils.tree_source import has_usable_tree, prioritize_tree_source

    assert not has_usable_tree(None)
    assert not has_usable_tree(Tree(newick="();"))
    assert has_usable_tree(Tree(newick="((A,B),(C,D));"))

    # 占位树不再能把来源抬成 concat → 汇总检验能如实报"无物种树".
    assert prioritize_tree_source(None, has_usable_tree(None), concat_label="concat") == "none"


def test_f5_small_dataset_cannot_report_success(tmp_path):
    """3 个基因组 × 每标记 <4 序列 → 无任何真树 → species_tree_source 必须是 none."""
    from markerfinder.modules.phylogenetic_inference import (
        CoalescentInference,
        PhylogeneticInferenceModule,
        SupermatrixInference,
    )

    cfg = PhylogeneticConfig(tmp_dir=str(tmp_path), coalescent_mode="post-filter")
    marker_genes = {
        "MK1": [{"id": f"g{i}", "seq": "MKV"} for i in range(3)],
    }
    sup = SupermatrixInference(cfg).run(marker_genes, [])
    coa = CoalescentInference(cfg, gene_trees_dir=str(tmp_path / "gt")).run(marker_genes, [])

    from markerfinder.utils.tree_source import has_usable_tree
    has_concat = has_usable_tree(sup.tree)
    coa_source = coa.species_tree_source
    has_coal = bool(coa_source and coa_source != "none"
                    and has_usable_tree(coa.species_tree))
    resolved = None
    from markerfinder.utils.tree_source import prioritize_tree_source
    resolved = prioritize_tree_source(coa_source if has_coal else None, has_concat,
                                      concat_label="concat")
    assert resolved == "none"


# ---------------------------------------------------------------------------
# F6 — TaxonomyConflictError is handled by the CLI (exit 3, no traceback)
# ---------------------------------------------------------------------------

def test_f6_load_taxonomy_conflict_exits_with_data_error(tmp_path):
    import markerfinder.cli.validation as cli_val

    f = tmp_path / "taxa.tsv"
    f.write_text("A\td__B\x00acteria;p__X\n", encoding="utf-8")
    args = types.SimpleNamespace(
        taxonomy_table=str(f), table_sep=None, taxonomy_format="table",
        taxonomy_delimiter_mode="reverse", ignore_malformed=False,
    )
    with pytest.raises(SystemExit) as exc:
        cli_val._load_taxonomy(args)
    assert exc.value.code == EXIT_DATA_ERROR


# ---------------------------------------------------------------------------
# F7 — stepwise runs clean their tmp dir; explicit tmp dirs are kept
# ---------------------------------------------------------------------------

def _bare_pipeline(config: PipelineConfig):
    import markerfinder.pipeline as pipeline_mod

    pipe = pipeline_mod.MarkerFinderPipeline.__new__(pipeline_mod.MarkerFinderPipeline)
    pipe.config = config
    return pipe


def test_f7_cleanup_removes_auto_tmp(tmp_path):
    cfg = PipelineConfig(output_dir=str(tmp_path / "out"))
    auto_tmp = tmp_path / "auto-tmp"
    auto_tmp.mkdir()
    (auto_tmp / "junk.txt").write_text("x", encoding="utf-8")
    cfg.tmp_dir = str(auto_tmp)
    cfg.tmp_dir_auto = True
    cfg.keep_tmp = False
    _bare_pipeline(cfg)._cleanup_tmp()
    assert not auto_tmp.exists()


def test_f7_cleanup_keeps_explicit_and_keep_tmp(tmp_path):
    explicit = tmp_path / "explicit-tmp"
    explicit.mkdir()
    cfg = PipelineConfig(output_dir=str(tmp_path / "out"), tmp_dir=str(explicit),
                         tmp_dir_auto=False, keep_tmp=False)
    _bare_pipeline(cfg)._cleanup_tmp()
    assert explicit.exists()

    auto = tmp_path / "auto2"
    auto.mkdir()
    cfg2 = PipelineConfig(output_dir=str(tmp_path / "out2"), tmp_dir=str(auto),
                          tmp_dir_auto=True, keep_tmp=True)
    _bare_pipeline(cfg2)._cleanup_tmp()
    assert auto.exists()


def _seed_state_objects(tmp_path: Path, keys) -> None:
    """Run_step 的原子性守卫要求上一步产物文件存在且非空, 先在磁盘上种好."""
    obj_dir = Path(tmp_path) / ".markerfinder" / "objects"
    obj_dir.mkdir(parents=True, exist_ok=True)
    for key in keys:
        (obj_dir / f"{key}.json").write_text("{}", encoding="utf-8")


def test_f7_run_step_cleanup_invoked(tmp_path, monkeypatch):
    """Run_step 的 report 分支结束时必须调用 _cleanup_tmp."""
    import markerfinder.pipeline as pipeline_mod

    calls = []
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_cleanup_tmp",
                        lambda self: calls.append(True))

    from markerfinder.models.pipeline_types import (
        HGTReport,
        MarkerSelectionResult as MSR,
        PhylogeneticResult as PR,
    )
    from markerfinder.models.report import PhaseContext, PlainTextReportOutput, ReportOutput

    _seed_state_objects(tmp_path, ("preprocessing", "marker_selection",
                                   "marker_sequences", "rank_map", "hgt_report",
                                   "precomputed_levels", "marker_genes", "phylo_result"))

    pipe = pipeline_mod.MarkerFinderPipeline.__new__(pipeline_mod.MarkerFinderPipeline)
    pipe.config = PipelineConfig(output_dir=str(tmp_path))

    scan_data = {"preprocessing": None, "marker_selection": MSR(), "marker_sequences": {},
                 "rank_map": {}}
    filter_data = {"hgt_report": HGTReport(), "precomputed_levels": {}, "marker_genes": {}}
    infer_data = {"phylo_result": PR()}
    state = {**scan_data, **filter_data, **infer_data,
             "completed_steps": ["scan", "filter", "infer"], "start_time": 1.0,
             "context": PhaseContext(quality_data=None)}
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_load_state",
                        lambda self: dict(state))
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_save_state",
                        lambda self, **kw: None)

    from markerfinder.models.report import RuntimeInfo

    def _fake_report(self, sd, fd, idata, ctx, runtime):
        return {"html_output": ReportOutput(html_path=""),
                "text_output": PlainTextReportOutput(file_paths=[])}
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "run_report", _fake_report)
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_write_run_config",
                        lambda self, s, e: None)
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_save_intermediates",
                        lambda self: None)

    pipe.run_step("report")
    assert calls, "report step must clean up its tmp dir"


# ---------------------------------------------------------------------------
# F8 — --min-occupancy reaches SelectionConfig and overrides adaptive params
# ---------------------------------------------------------------------------

def _parse_with_genome(tmp_path, *extra):
    from markerfinder.cli.parser import _build_parser

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir(exist_ok=True)
    (genome_dir / "g.faa").write_text(">a\nMKV\n", encoding="utf-8")
    parser = _build_parser()
    return parser, parser.parse_args(["-i", str(genome_dir), "-o", str(tmp_path / "out"),
                                      *extra])


def test_f8_min_occupancy_new_flag(tmp_path):
    from markerfinder.cli.config_build import _build_pipeline_config

    parser, args = _parse_with_genome(tmp_path, "--min-occupancy", "0.55")
    cfg = _build_pipeline_config(args, parser)
    assert cfg.selection_config.user_min_occupancy == 0.55
    assert cfg.selection_config.min_occupancy == 0.55


def test_f8_legacy_min_marker_coverage_alias(tmp_path):
    from markerfinder.cli.config_build import _build_pipeline_config

    parser, args = _parse_with_genome(tmp_path, "--min-marker-coverage", "0.4")
    cfg = _build_pipeline_config(args, parser)
    assert cfg.selection_config.user_min_occupancy == 0.4


def test_f8_min_occupancy_defaults_unchanged(tmp_path):
    from markerfinder.cli.config_build import _build_pipeline_config

    parser, args = _parse_with_genome(tmp_path)
    cfg = _build_pipeline_config(args, parser)
    assert cfg.selection_config.user_min_occupancy is None
    assert cfg.selection_config.min_occupancy == 0.75


def test_f8_user_occupancy_overrides_adaptive_params(tmp_path):
    import markerfinder.pipeline as pipeline_mod
    from markerfinder.models.marker import SelectedMarkerSet
    from markerfinder.models.pipeline_types import (
        AdaptiveParams,
        MarkerSelectionResult,
        PreprocessingResult,
    )

    cfg = PipelineConfig(
        output_dir=str(tmp_path),
        selection_config=SelectionConfig(user_min_occupancy=0.55),
    )
    genomes = []
    pipe = pipeline_mod.MarkerFinderPipeline.__new__(pipeline_mod.MarkerFinderPipeline)
    pipe.config = cfg

    adaptive = AdaptiveParams(min_occupancy=0.35)
    pre = PreprocessingResult(adaptive_params=adaptive)
    pipe.mag_optimization = types.SimpleNamespace(run=lambda g: pre)

    # Marker_mode hmm without hmm dir would raise; use gtdb with a stubbed selector.
    cfg.selection_config.marker_mode = "gtdb_tk"
    cfg.selection_config.gtdb_markers_dir = "unused"
    ms_result = MarkerSelectionResult(marker_set=SelectedMarkerSet())
    stub = types.SimpleNamespace()
    stub.run_gtdb_tk = lambda *a, **k: (ms_result, {}, None)
    pipe.marker_selector = stub
    pipe._resolve_orthologs = lambda genomes, ms, seqs: seqs

    context = pipeline_mod.PhaseContext()
    pipe.run_scan(genomes, context)
    assert context.adaptive_params.min_occupancy == 0.55


# ---------------------------------------------------------------------------
# F9 — scan-phase gene trees land in the output cache used by infer
# ---------------------------------------------------------------------------

def test_f9_run_gtdb_tk_uses_output_gene_trees_dir(tmp_path, monkeypatch):
    import markerfinder.modules.marker_selection as ms_mod
    import markerfinder.utils.gtdb_tk_markers as gtm
    from markerfinder.models.genome import GeneState, Genome, OccupancyMatrix

    seen = {}

    matrix = OccupancyMatrix(genomes=["g1", "g2", "g3", "g4"], cogs=["MK1"])
    for g in matrix.genomes:
        matrix.set(g, "MK1", GeneState.SINGLE_COPY)

    mseq = {"MK1": [{"id": g, "genome_id": g, "seq": "MKVL"} for g in matrix.genomes]}
    occ = {"MK1": 1.0}

    monkeypatch.setattr(gtm, "load_gtdb_markers",
                        lambda markers_dir, genome_ids=None: (matrix, mseq, occ))

    def fake_build(seqs, out_dir, method="fasttree", min_tips=4, trim_threshold=0.2, cpus=1):
        seen["out_dir"] = out_dir
        return {"MK1": None}

    monkeypatch.setattr(gtm, "build_gene_trees", fake_build)

    cfg = SelectionConfig(marker_mode="gtdb_tk", gtdb_markers_dir=str(tmp_path))
    module = ms_mod.AdaptiveMarkerSelectionModule(cfg)
    genomes = [Genome(id=f"g{i}", fasta_path=str(tmp_path / f"g{i}.faa")) for i in range(4)]

    expected = str(Path(tmp_path / "out") / "Phase4_trees" / "gene_trees")
    module.run_gtdb_tk(genomes, markers_dir=str(tmp_path), tmp_dir=str(tmp_path),
                       adaptive_params=None, gene_trees_dir=expected)
    assert seen["out_dir"] == expected


# ---------------------------------------------------------------------------
# F11 — `filter --resume` without -i fails with an actionable error
# ---------------------------------------------------------------------------

def test_f11_resume_without_input_returns_arg_error(tmp_path, caplog):
    from markerfinder.cli.commands import _execute_pipeline

    args = types.SimpleNamespace(resume=True, redo=False, output=str(tmp_path),
                                 input=None, command="filter")
    rc = _execute_pipeline("filter", args, object(), None, None, "external table",
                           __import__("logging").getLogger(__name__))
    assert rc == EXIT_ARG_ERROR


# ---------------------------------------------------------------------------
# F13 — unbalanced parentheses / empty tips are critical validation errors
# ---------------------------------------------------------------------------

def test_f13_unbalanced_newick_is_critical(tmp_path):
    from markerfinder.exceptions import TreeValidationError
    from markerfinder.validation import validate_tree_file

    f = tmp_path / "bad.nwk"
    f.write_text("((A,B);", encoding="utf-8")
    with pytest.raises(TreeValidationError, match="Critical"):
        validate_tree_file(str(f))


def test_f13_empty_tip_is_critical(tmp_path):
    from markerfinder.exceptions import TreeValidationError
    from markerfinder.validation import validate_tree_file

    f = tmp_path / "empty_tip.nwk"
    f.write_text("((:0.1,B:0.2),(C:0.3,D:0.4));", encoding="utf-8")
    with pytest.raises(TreeValidationError, match="Critical"):
        validate_tree_file(str(f))


# ---------------------------------------------------------------------------
# F15 — self-test version comparison is numeric, not lexicographic
# ---------------------------------------------------------------------------

def test_f15_version_tuple_numeric_compare():
    from markerfinder.cli.self_test import _version_tuple

    assert _version_tuple("1.100") > _version_tuple("1.99")
    assert _version_tuple("2.0") > _version_tuple("1.81")
    assert _version_tuple("1.9") < _version_tuple("1.10")
    assert _version_tuple("abc") == (0,)


# ---------------------------------------------------------------------------
# F17 — detect_scope_rank type hints resolve (Iterable imported)
# ---------------------------------------------------------------------------

def test_f17_scope_rank_type_hints_resolve():
    import typing

    import markerfinder.taxonomy as taxonomy

    hints = typing.get_type_hints(taxonomy.detect_scope_rank)
    assert "tip_labels" in hints


# ---------------------------------------------------------------------------
# F19 — completed_steps dedup + sane run duration for missing start_time
# ---------------------------------------------------------------------------

def test_f19_run_config_duration_guard(tmp_path):
    import time

    pipe = _bare_pipeline(PipelineConfig(output_dir=str(tmp_path)))
    pipe._write_run_config(0.0, time.time())
    data = json.loads((Path(tmp_path) / "Phase5_metadata" / "run_config.json").read_text())
    assert data["run_duration_seconds"] == 0.0


def test_f19_completed_steps_dedup(tmp_path, monkeypatch):
    import markerfinder.pipeline as pipeline_mod
    from markerfinder.models.pipeline_types import (
        HGTReport,
        MarkerSelectionResult as MSR,
    )
    from markerfinder.models.report import PhaseContext

    pipe = pipeline_mod.MarkerFinderPipeline.__new__(pipeline_mod.MarkerFinderPipeline)
    pipe.config = PipelineConfig(output_dir=str(tmp_path))

    _seed_state_objects(tmp_path, ("marker_sequences",))

    saved = {}
    state = {
        "preprocessing": None, "marker_selection": MSR(), "marker_sequences": {},
        "rank_map": {}, "hgt_report": HGTReport(), "precomputed_levels": {},
        "marker_genes": {}, "genomes": [], "context": PhaseContext(),
        "taxonomy_map": None, "start_time": 1.0,
        "completed_steps": ["scan", "filter", "infer"],
    }
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_load_state",
                        lambda self: dict(state))
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_save_state",
                        lambda self, **kw: saved.update(kw))
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "_cleanup_tmp",
                        lambda self: None)
    monkeypatch.setattr(pipeline_mod.MarkerFinderPipeline, "run_filter",
                        lambda self, g, sd, ctx, taxonomy_map=None: {
                            "hgt_report": HGTReport(), "precomputed_levels": {},
                            "marker_genes": {}})

    # 重跑 filter: completed 里已有 infer, 追加 filter 后不得产生重复项.
    pipe.run_step("filter", genomes=[])
    assert saved["completed_steps"].count("filter") == 1


# ---------------------------------------------------------------------------
# F20 — parse_fasta warns on duplicate ids (keeps last occurrence)
# ---------------------------------------------------------------------------

def test_f20_parse_fasta_duplicate_id_warns_and_keeps_last(tmp_path, caplog):
    from markerfinder.utils.io import parse_fasta

    f = tmp_path / "dup.faa"
    f.write_text(">a\nMKV\n>a\nLLLL\n", encoding="utf-8")
    with caplog.at_level("WARNING", logger="markerfinder.utils.io"):
        seqs = parse_fasta(str(f))
    assert seqs == {"a": "LLLL"}
    assert "Duplicate FASTA id 'a'" in caplog.text


# ---------------------------------------------------------------------------
# F21 — HTML risk buckets follow configured thresholds, UNKNOWN handled
# ---------------------------------------------------------------------------

def test_f21_risk_buckets_use_configured_thresholds(tmp_path):
    from markerfinder.modules.report_generator import InteractiveReportGenerator

    report = HGTReport(marker_evaluations=[
        HGTEvaluation(marker_id="m1", overall_risk=0.30, level=MarkerLevel.LEVEL_1,
                      confidence="high"),
        HGTEvaluation(marker_id="m2", overall_risk=0.50, level=MarkerLevel.LEVEL_2,
                      confidence="high"),
    ])
    gen = InteractiveReportGenerator(ReportConfig(output_dir=str(tmp_path),
                                                  level1_max=0.40, level2_max=0.70))
    out = gen.generate(_marker_result(), report, PhylogeneticResult(), None,
                       RuntimeInfo(duration=1.0))
    html = Path(out.html_path).read_text(encoding="utf-8")
    assert "0.00 - 0.40 (Level 1)" in html
    assert "0.40 - 0.70 (Level 2)" in html
    assert "≥ 0.70 (Level 3)" in html
    # Risk=0.30 在 0.40 阈值下属于 Level 1 桶(旧常量 0.25 会把它算进 Level 2).
    level1_bar = [ln for ln in html.splitlines() if "0.00 - 0.40 (Level 1)" in ln][0]
    assert "title=\"1 markers\"" in level1_bar


# ---------------------------------------------------------------------------
# F14 — UNKNOWN carries a specific skip reason (not a one-size-fits-all note)
# ---------------------------------------------------------------------------

def test_f14_skip_reason_distinguishes_causes(tmp_path):
    from markerfinder.modules.hgt_filter import HGTFilterModule

    # 提供了物种树但基因树不可得(<4 序列) → 成因必须写明基因树问题,
    # 而不是一概归咎"未提供参照".
    mod = HGTFilterModule(HGTConfig(), tmp_dir=str(tmp_path))
    report, _, _far = mod.run(
        ["mk1"],
        marker_sequences={"mk1": [{"id": "a", "seq": "M"}]},
        species_tree=Tree(newick="((a,b),(c,d));"),
    )
    notes = report.marker_evaluations[0].notes
    assert report.marker_evaluations[0].level == MarkerLevel.UNKNOWN
    assert "gene tree" in notes
