import pytest
from unittest.mock import patch, MagicMock

from markerfinder.config import PhylogeneticConfig
from markerfinder.models.tree import Tree, TreeRecommendation
from markerfinder.models.pipeline_types import PhylogeneticResult, SupermatrixResult, CoalescentResult
from markerfinder.modules.phylogenetic_inference import (
    SupermatrixInference, CoalescentInference, ConflictDetector, PhylogeneticInferenceModule,
)


class TestSupermatrixInference:
    def test_run(self):
        # 审查修复 F5: 无可用比对时不再返回占位树 ";"(会被真值判定当成有效
        # 物种树造成假成功), 而是 tree=None; 分区信息保留.
        sm = SupermatrixInference(PhylogeneticConfig())
        r = sm.run({"C1": [], "C2": []}, [])
        assert r is not None
        assert r.partition is not None
        assert r.tree is None


class TestCoalescentInference:
    def test_run_empty(self):
        ci = CoalescentInference(PhylogeneticConfig())
        r = ci.run({"C1": [], "C2": []}, [])
        assert r.n_total_genes == 0  # Empty seqs → no gene trees

    def test_run_with_seqs(self):
        ci = CoalescentInference(PhylogeneticConfig(use_fasttree=False))
        seqs = [{"id": f"s{i}", "seq": "ACGT" * 30} for i in range(6)]
        r = ci.run({"C1": seqs}, [])
        assert r.n_total_genes == 0  # Mafft not available → 0 gene trees

    def test_gene_tree_respects_config_ufboot(self, tmp_path):
        config = PhylogeneticConfig(
            ufboot_replicates=2000,
            fast_mode=False,
            use_fasttree=False,
            iqtree_timeout=1200,
        )
        ci = CoalescentInference(config, gene_trees_dir=str(tmp_path))
        seqs = [{"id": f"s{i}", "seq": "ACGT" * 30} for i in range(6)]
        with patch("markerfinder.modules.phylogenetic_inference._run_iqtree3") as mock_iqtree:
            mock_iqtree.return_value = Tree(newick="(s0,s1,(s2,s3,s4,s5));")
            with patch("markerfinder.modules.phylogenetic_inference._run_mafft", return_value=True):
                with patch("markerfinder.modules.phylogenetic_inference._run_trimal", return_value=False):
                    r = ci.run({"C1": seqs}, [])
                    mock_iqtree.assert_called_once()
                    _, kwargs = mock_iqtree.call_args
                    assert kwargs["ufboot"] == 2000
                    assert kwargs["fast_mode"] is False
                    assert kwargs["timeout"] == 1200


class TestConflictDetector:
    def test_detect(self):
        cd = ConflictDetector()
        # A 2-tip tree is not a legal reference — conflict
        # Detection is refused (all metrics None) instead of scored.
        r = cd.detect_conflicts(Tree(newick="(A,B);"), Tree(newick="(A,B);"), {})
        assert r.normalized_rf is None
        assert "reference illegal" in r.topology_note

    def test_detect_legal_trees(self):
        cd = ConflictDetector()
        r = cd.detect_conflicts(
            Tree(newick="((A,B),(C,D));"), Tree(newick="((A,B),(C,D));"), {}
        )
        assert r.normalized_rf == 0.0
        assert r.quartet_agreement == 1.0


class TestPhylogeneticInferenceModule:
    def test_run_both(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig(coalescent_mode="post-filter"))
        r = mod.run({"C1": []}, [])
        assert r.supermatrix is not None

    def test_run_coalescent_off(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig(coalescent_mode="off"))
        r = mod.run({"C1": []}, [])
        assert r.supermatrix is not None
        assert r.coalescent is None

    def test_run_empty_markers(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig())
        r = mod.run({}, [])
        assert r.supermatrix is None

    def test_recommend(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig(coalescent_mode="off"))
        r = PhylogeneticResult(supermatrix=SupermatrixResult(tree=Tree(newick="(A,B);")))
        rec = mod.recommend_tree(r)
        assert rec.recommended_tree is not None
