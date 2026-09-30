"""Additional unit tests for phylogenetic inference helpers and modules."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from markerfinder.config import PhylogeneticConfig
from markerfinder.models.tree import Tree
from markerfinder.models.pipeline_types import (
    ConflictReport,
    PhylogeneticInferenceError,
    PhylogeneticResult,
    SupermatrixResult,
    CoalescentResult,
)
from markerfinder.modules.phylogenetic_inference import (
    _concatenate_fasta_alignments,
    _run_astral,
    _run_fasttree,
    _run_iqtree3,
    _run_mafft,
    _run_trimal,
    _write_fasta_for_seqs,
    _write_partition_nexus,
    CoalescentInference,
    ConflictDetector,
    PhylogeneticInferenceModule,
    SupermatrixInference,
)


def _simple_seqs():
    return [{"id": f"s{i}", "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in range(6)]


def _aln_fasta(ids):
    lines = []
    for sid in ids:
        lines.append(f">{sid}")
        lines.append("ACDEFGHIKLMNPQRSTVWY" * 3)
    return "\n".join(lines) + "\n"


def _tree_newick(ids):
    if len(ids) == 2:
        return f"({ids[0]},{ids[1]});"
    return f"({ids[0]},{ids[1]},({ids[2]},{ids[3]}));"


class TestWriteFastaForSeqs:
    def test_writes_file(self, tmp_path):
        out = tmp_path / "out.faa"
        _write_fasta_for_seqs([{"id": "a", "seq": "ACGT"}], str(out))
        text = out.read_text()
        assert ">a" in text
        assert "ACGT" in text


class TestMafftTrimalHelpers:
    def test_run_mafft_writes_output(self, tmp_path):
        inp = tmp_path / "in.faa"
        out = tmp_path / "out.aln"
        inp.write_text(">s1\nACGT\n")

        def runner(cmd, **kwargs):
            if cmd[0] == "mafft":
                out.resolve().write_text(">s1\nACGT\n", encoding="utf-8")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            assert _run_mafft(str(inp), str(out), cpus=1) is True

    def test_run_mafft_not_found(self, tmp_path, caplog):
        # 设计契约: MAFFT 缺失时 _run_mafft 优雅返回 False(上层据此跳过该标记 /
        # Coalescent 得到 0 基因树), 而非抛异常中断整条流程. 这与
        # Test_phylogenetic.TestCoalescentInference.test_run_with_seqs 一致.
        # 失败会通过 warning 日志可见, 不属"静默降级为错误结论"(对照
        # 真正修复的物种树占位/ HGT 冒充 Level3 / ABSENT / 质量空解析).
        import logging

        with caplog.at_level(
            logging.WARNING,
            logger="markerfinder.modules.phylogenetic_inference",
        ):
            with patch("subprocess.run", side_effect=FileNotFoundError):
                assert (
                    _run_mafft(str(tmp_path / "in"), str(tmp_path / "out")) is False
                )
            assert any("MAFFT not found" in r.message for r in caplog.records)

    def test_run_trimal_writes_output(self, tmp_path):
        inp = tmp_path / "in.aln"
        out = tmp_path / "out.trim"
        inp.write_text(">s1\nACGT\n")

        def runner(cmd, **kwargs):
            if cmd[0] == "trimal":
                out.write_text(inp.read_text())
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            assert _run_trimal(str(inp), str(out)) is True

    def test_run_trimal_not_found(self, tmp_path):
        with patch("subprocess.run", side_effect=FileNotFoundError):
            assert _run_trimal(str(tmp_path / "in"), str(tmp_path / "out")) is False


class TestTreeBuilders:
    def test_run_fasttree_writes_tree(self, tmp_path):
        aln = tmp_path / "aln"
        tree = tmp_path / "tree.nwk"
        aln.write_text(">s1\nACGT\n")

        def runner(cmd, **kwargs):
            tree.write_text("(s1);\n")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = _run_fasttree(str(aln), str(tree))
            assert result is not None
            assert result.n_tips >= 1

    def test_run_fasttree_not_found(self, tmp_path):
        with patch("subprocess.run", side_effect=FileNotFoundError):
            assert _run_fasttree(str(tmp_path / "aln"), str(tmp_path / "out")) is None

    def test_run_iqtree3_writes_tree(self, tmp_path):
        aln = tmp_path / "concat.fasta"
        prefix = str(tmp_path / "test")
        aln.write_text(">s1\nACGT\n")

        def runner(cmd, **kwargs):
            Path(f"{prefix}.treefile").write_text("(s1);\n")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = _run_iqtree3(str(aln), prefix, fast_mode=True)
            assert result is not None

    def test_run_iqtree3_failed(self, tmp_path):
        with patch("subprocess.run", side_effect=FileNotFoundError):
            assert _run_iqtree3(str(tmp_path / "aln"), str(tmp_path / "pre")) is None

    def test_run_astral_writes_tree(self, tmp_path):
        gene_trees = tmp_path / "gt.nwk"
        out_tree = tmp_path / "astral.nwk"
        gene_trees.write_text("(s1,s2);\n(s1,s2);\n")

        def runner(cmd, **kwargs):
            out_tree.write_text("(s1,s2);\n")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = _run_astral(str(gene_trees), str(out_tree))
            assert result is not None

    def test_run_astral_not_found(self, tmp_path):
        with patch("subprocess.run", side_effect=FileNotFoundError):
            assert _run_astral(str(tmp_path / "in"), str(tmp_path / "out")) is None


class TestConcatenateAndPartition:
    def test_concatenate_alignments(self, tmp_path):
        aln1 = tmp_path / "a.aln"
        aln2 = tmp_path / "b.aln"
        aln1.write_text(">s1\nAAAA\n>s2\nAAAA\n")
        aln2.write_text(">s1\nCCCC\n>s2\nCCCC\n")
        out = tmp_path / "concat.fasta"
        species, sites = _concatenate_fasta_alignments(
            {"M1": str(aln1), "M2": str(aln2)}, str(out)
        )
        assert sorted(species) == ["s1", "s2"]
        assert sites == 8
        text = out.read_text()
        assert "AAAACCCC" in text

    def test_write_partition_nexus(self, tmp_path):
        out = tmp_path / "part.nex"
        _write_partition_nexus({"M1": 100, "M2": 200}, str(out))
        text = out.read_text()
        assert "charset M1 = 1-100;" in text
        assert "charset M2 = 101-300;" in text


class TestSupermatrixInference:
    def _mock_run(self, tmp_path, ids):
        def runner(cmd, **kwargs):
            exe = cmd[0].lower()
            if "mafft" in exe:
                out_path = kwargs.get("stdout") or (cmd[-1].replace(".faa", ".aln"))
                if hasattr(out_path, "write"):
                    out_path.write(_aln_fasta(ids))
                else:
                    Path(out_path).write_text(_aln_fasta(ids))
            elif "trimal" in exe:
                in_path = cmd[cmd.index("-in") + 1]
                out_path = cmd[cmd.index("-out") + 1]
                Path(out_path).write_text(Path(in_path).read_text())
            elif "iqtree3" in exe or "iqtree" in exe:
                prefix = cmd[cmd.index("-pre") + 1]
                Path(f"{prefix}.treefile").write_text(_tree_newick(ids))
            return MagicMock(returncode=0)
        return runner

    def test_supermatrix_runs(self, tmp_path):
        ids = ["s0", "s1", "s2", "s3"]
        seqs = [{"id": i, "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in ids]
        cfg = PhylogeneticConfig(tmp_dir=str(tmp_path), fast_mode=True)
        sm = SupermatrixInference(cfg)
        with patch("subprocess.run", side_effect=self._mock_run(tmp_path, ids)):
            result = sm.run({"M1": seqs}, [])
        assert result.tree is not None
        assert result.partition is not None

    def test_supermatrix_falls_back_to_fasttree(self, tmp_path):
        ids = ["s0", "s1", "s2", "s3"]
        seqs = [{"id": i, "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in ids]
        cfg = PhylogeneticConfig(tmp_dir=str(tmp_path))
        sm = SupermatrixInference(cfg)

        def runner(cmd, **kwargs):
            exe = cmd[0].lower()
            if "iqtree" in exe:
                # No treefile -> failure
                return MagicMock(returncode=1)
            if "mafft" in exe:
                out = kwargs.get("stdout")
                out.write(_aln_fasta(ids))
            elif "trimal" in exe:
                Path(cmd[cmd.index("-out") + 1]).write_text(
                    Path(cmd[cmd.index("-in") + 1]).read_text())
            elif "fasttree" in exe or cmd[0] == "FastTree":
                Path(cmd[cmd.index("-out") + 1]).write_text(_tree_newick(ids))
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = sm.run({"M1": seqs}, [])
        assert result.tree is not None

    def test_supermatrix_raises_when_all_tree_builders_fail(self, tmp_path):
        # 当 IQ-TREE3 与 FastTree 均缺失/失败时, 不应再静默返回空树;
        # 必须显式抛出 PhylogeneticInferenceError 让调用方感知工具失败.
        ids = ["s0", "s1", "s2", "s3"]
        seqs = [{"id": i, "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in ids]
        cfg = PhylogeneticConfig(tmp_dir=str(tmp_path))
        sm = SupermatrixInference(cfg)

        def runner(cmd, **kwargs):
            exe = cmd[0].lower()
            if "mafft" in exe:
                kwargs["stdout"].write(_aln_fasta(ids))
            elif "trimal" in exe:
                Path(cmd[cmd.index("-out") + 1]).write_text(
                    Path(cmd[cmd.index("-in") + 1]).read_text())
            elif "iqtree" in exe or "fasttree" in exe or cmd[0] == "FastTree":
                # Both tree builders unavailable -> propagate FileNotFoundError
                # (the helpers swallow it internally and return None).
                raise FileNotFoundError("tree builder missing")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            with pytest.raises(PhylogeneticInferenceError):
                sm.run({"M1": seqs}, [])


class TestCoalescentInference:
    def _mock_run(self, tmp_path, ids):
        def runner(cmd, **kwargs):
            exe = cmd[0].lower()
            if "mafft" in exe:
                kwargs["stdout"].write(_aln_fasta(ids))
            elif "trimal" in exe:
                Path(cmd[cmd.index("-out") + 1]).write_text(
                    Path(cmd[cmd.index("-in") + 1]).read_text())
            elif "fasttree" in exe or cmd[0] == "FastTree":
                out = cmd[cmd.index("-out") + 1]
                Path(out).write_text(_tree_newick(ids))
            elif "astral" in exe:
                out = cmd[cmd.index("-o") + 1]
                Path(out).write_text(_tree_newick(ids))
            return MagicMock(returncode=0)
        return runner

    def test_coalescent_runs_with_fasttree(self, tmp_path):
        ids = ["s0", "s1", "s2", "s3"]
        seqs = [{"id": i, "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in ids]
        cfg = PhylogeneticConfig(tmp_dir=str(tmp_path), use_fasttree=True)
        ci = CoalescentInference(cfg, gene_trees_dir=str(tmp_path / "gts"))
        with patch("subprocess.run", side_effect=self._mock_run(tmp_path, ids)):
            result = ci.run({"M1": seqs}, [])
        assert result.n_total_genes == 1
        assert result.species_tree is not None

    def test_coalescent_reuses_cached_gene_tree(self, tmp_path):
        ids = ["s0", "s1", "s2", "s3"]
        seqs = [{"id": i, "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in ids]
        cache_dir = tmp_path / "gts"
        cache_dir.mkdir()
        (cache_dir / "M1.nwk").write_text(_tree_newick(ids))
        cfg = PhylogeneticConfig(tmp_dir=str(tmp_path), use_fasttree=True)
        ci = CoalescentInference(cfg, gene_trees_dir=str(cache_dir))

        def runner(cmd, **kwargs):
            if cmd[0].lower() == "astral":
                out = cmd[cmd.index("-o") + 1]
                Path(out).write_text(_tree_newick(ids))
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = ci.run({"M1": seqs}, [])
        assert result.n_total_genes == 1

    def test_coalescent_astral_failure(self, tmp_path):
        ids = ["s0", "s1", "s2", "s3"]
        seqs = [{"id": i, "seq": "ACDEFGHIKLMNPQRSTVWY" * 3} for i in ids]
        cfg = PhylogeneticConfig(tmp_dir=str(tmp_path), use_fasttree=True)
        ci = CoalescentInference(cfg)

        def runner(cmd, **kwargs):
            exe = cmd[0].lower()
            if "mafft" in exe:
                kwargs["stdout"].write(_aln_fasta(ids))
            elif "trimal" in exe:
                Path(cmd[cmd.index("-out") + 1]).write_text(
                    Path(cmd[cmd.index("-in") + 1]).read_text())
            elif "fasttree" in exe or cmd[0] == "FastTree":
                out = cmd[cmd.index("-out") + 1]
                Path(out).write_text(_tree_newick(ids))
            elif "astral" in exe:
                raise FileNotFoundError("astral")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            result = ci.run({"M1": seqs}, [])
        assert result.species_tree_source == "none"


class TestConflictDetector:
    def test_detect_conflicts(self):
        cd = ConflictDetector()
        report = cd.detect_conflicts(
            Tree(newick="((A,B),(C,D));"),
            Tree(newick="((A,B),(C,D));"),
            {},
        )
        assert isinstance(report, ConflictReport)
        assert report.normalized_rf == 0.0
        assert report.quartet_agreement == 1.0
        assert "identical" in report.topology_note

    def test_topology_note_moderate(self):
        cd = ConflictDetector()
        note = cd._topology_note("t1", "t2", 0.2, 0.8)
        assert "Moderate" in note

    def test_topology_note_strong(self):
        cd = ConflictDetector()
        note = cd._topology_note("t1", "t2", 0.5, 0.5)
        assert "Strong" in note

    def test_topology_note_negative_rf(self):
        cd = ConflictDetector()
        note = cd._topology_note("", "", -1, None)
        assert "could not be computed" in note


class TestPhylogeneticInferenceModule:
    def test_recommend_with_conflict(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig())
        pr = PhylogeneticResult(
            supermatrix=SupermatrixResult(tree=Tree(newick="(A,B);")),
            coalescent=CoalescentResult(species_tree=Tree(newick="(A,B);"), species_tree_source="astral"),
            conflict_report=ConflictReport(normalized_rf=0.05),
        )
        rec = mod.recommend_tree(pr)
        assert rec.confidence == "high"

    def test_recommend_moderate_conflict(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig())
        pr = PhylogeneticResult(
            supermatrix=SupermatrixResult(tree=Tree(newick="(A,B);")),
            coalescent=CoalescentResult(species_tree=Tree(newick="(A,B);"), species_tree_source="astral"),
            conflict_report=ConflictReport(normalized_rf=0.2),
        )
        rec = mod.recommend_tree(pr)
        assert rec.confidence == "medium"

    def test_recommend_strong_conflict(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig())
        pr = PhylogeneticResult(
            supermatrix=SupermatrixResult(tree=Tree(newick="(A,B);")),
            coalescent=CoalescentResult(species_tree=Tree(newick="(A,B);"), species_tree_source="astral"),
            conflict_report=ConflictReport(normalized_rf=0.5),
        )
        rec = mod.recommend_tree(pr)
        assert rec.confidence == "low"
