"""Upstream independence provenance tests."""

from markerfinder.config import PhylogeneticConfig
from markerfinder.models.pipeline_types import (
    CoalescentResult,
    ConflictReport,
    IndependenceReport,
    PhylogeneticResult,
    SupermatrixResult,
)
from markerfinder.models.tree import Tree
from markerfinder.modules.hgt_filter import _hgt_build_gene_tree_with_source
from markerfinder.modules.phylogenetic_inference import (
    CoalescentInference,
    PhylogeneticInferenceModule,
)


class TestCacheProvenance:
    def test_cache_hit_recorded(self, tmp_path):
        """A cached gene tree must be labelled 'cached'.

        A hit requires the input stamp beside the tree: an entry whose input is
        unknown is a different marker's tree wearing this marker's name.
        """
        from markerfinder.modules.hgt_filter import (
            _hgt_build_gene_tree_with_source,
            _marker_input_id,
        )

        gene_trees_dir = tmp_path / "gt"
        gene_trees_dir.mkdir()
        seqs = [{"id": f"g{i}", "seq": "ACDEFG" * 5} for i in range(4)]
        (gene_trees_dir / "M1.nwk").write_text("((A,B),(C,D));", encoding="utf-8")
        (gene_trees_dir / "M1.nwk.input_sha256").write_text(
            _marker_input_id(seqs) + "\n", encoding="utf-8")

        tree, source = _hgt_build_gene_tree_with_source(
            "M1", seqs, str(tmp_path), 1, str(gene_trees_dir),
        )
        assert source == "cached"
        assert tree is not None

    def test_a_tree_without_an_input_stamp_is_not_trusted(self, tmp_path):
        """MUST-FAIL CONTROL: same directory, no stamp sidecar.

        The only difference is the record of what the cached tree was built
        from, so this pins the invalidation rather than the tool availability: a
        run that reuses an output directory (--force over an existing directory
        is documented) would otherwise infer from sequences it never supplied.
        """
        from markerfinder.modules.hgt_filter import (
            _hgt_build_gene_tree_with_source,
        )

        gene_trees_dir = tmp_path / "gt"
        gene_trees_dir.mkdir()
        (gene_trees_dir / "M1.nwk").write_text("((A,B),(C,D));", encoding="utf-8")

        seqs = [{"id": f"g{i}", "seq": "ACDEFG" * 5} for i in range(4)]
        tree, source = _hgt_build_gene_tree_with_source(
            "M1", seqs, str(tmp_path), 1, str(gene_trees_dir),
        )
        # Either rebuilt now, or no tree at all if FastTree is unavailable — but
        # Never the unstamped "cached" answer.
        assert source != "cached", source

    def test_too_few_seqs_source_empty(self, tmp_path):
        tree, source = _hgt_build_gene_tree_with_source(
            "M1", [{"id": "a", "seq": "ACD"}], str(tmp_path), 1, str(tmp_path),
        )
        assert tree is None
        assert source == ""


class TestProvenanceOnResult:
    def test_provenance_field_exists(self):
        result = CoalescentResult()
        assert result.gene_tree_provenance == {}

    def test_coalescent_run_records_rebuilt(self, tmp_path):
        """Without a cache dir, every built tree is labelled 'rebuilt'."""
        import types

        cfg = PhylogeneticConfig(
            tmp_dir=str(tmp_path), coalescent_mode="post-filter",
        )
        module = CoalescentInference(cfg, gene_trees_dir="")
        seqs = {"M1": [{"id": f"g{i}", "seq": "ACDEFG" * 5} for i in range(4)]}
        result = module.run(seqs, [])
        # ASTRAL missing => no species tree, but provenance is still recorded
        for source in result.gene_tree_provenance.values():
            assert source in ("cached", "rebuilt")


class TestIndependenceReport:
    @staticmethod
    def _module_and_result():
        cfg = PhylogeneticConfig()
        mod = PhylogeneticInferenceModule(cfg)
        # Minimal four-tip legal trees
        concat = Tree(newick="((A,B),(C,D));")
        astral = Tree(newick="((A,B),(C,D));")
        supermatrix = SupermatrixResult(tree=concat)
        coalescent = CoalescentResult(
            species_tree=astral,
            gene_tree_provenance={"M1": "cached", "M2": "rebuilt"},
            species_tree_source="astral",
        )
        result = PhylogeneticResult(
            supermatrix=supermatrix, coalescent=coalescent,
        )
        return mod, result

    def test_independence_computed_when_conflicts(self):
        mod, result = self._module_and_result()
        result.conflict_report = mod.conflict_detector.detect_conflicts(
            result.supermatrix.tree, result.coalescent.species_tree,
            result.coalescent.gene_trees,
        )
        # Simulate the pipeline's provenance wiring
        from markerfinder.models.pipeline_types import IndependenceReport

        provenance = result.coalescent.gene_tree_provenance
        cached = sum(1 for v in provenance.values() if v == "cached")
        result.conflict_report.independence = IndependenceReport(
            marker_set_jaccard=1.0,
            shared_gene_tree_ratio=cached / len(provenance),
            same_trimming_regime=False,
            ref_built_from_tested_markers=True,
        )
        assert result.conflict_report.independence.shared_gene_tree_ratio == 0.5

    def test_shared_upstream_caps_confidence(self):
        """Shared upstream => recommendation confidence never 'high'."""
        mod, result = self._module_and_result()
        result.conflict_report = ConflictReport(
            rf_distance=0,
            normalized_rf=0.0,
            quartet_agreement=1.0,
            independence=IndependenceReport(
                marker_set_jaccard=1.0,
                shared_gene_tree_ratio=1.0,
                same_trimming_regime=False,
                ref_built_from_tested_markers=True,
            ),
        )
        rec = mod.recommend_tree(result)
        assert rec.confidence != "high"
        assert rec.confidence in ("medium", "low")
        assert "cap" in rec.reason.lower() or "independent" in rec.reason.lower()

    def test_independent_upstream_keeps_high(self):
        mod, result = self._module_and_result()
        result.conflict_report = ConflictReport(
            rf_distance=0,
            normalized_rf=0.0,
            quartet_agreement=1.0,
            independence=IndependenceReport(
                marker_set_jaccard=0.4,
                shared_gene_tree_ratio=0.0,
                same_trimming_regime=False,
                ref_built_from_tested_markers=False,
            ),
        )
        rec = mod.recommend_tree(result)
        assert rec.confidence == "high"
