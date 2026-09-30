"""Additional unit tests for HGT filtering."""

from unittest.mock import MagicMock, patch

import pytest

from markerfinder.config import HGTConfig
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.tree import Tree
from markerfinder.models.pipeline_types import PhyloStepResult
from markerfinder.modules.hgt_filter import (
    PhylogeneticHGTDetector,
    HGTDecisionEngine,
    HGTFilterModule,
    _extract_id,
    _extract_seq,
)


class TestSeqExtraction:
    def test_extract_seq_dict(self):
        assert _extract_seq({"id": "a", "seq": "ACGT"}) == "ACGT"

    def test_extract_seq_object(self):
        class Obj:
            seq = "TGCA"
        assert _extract_seq(Obj()) == "TGCA"

    def test_extract_id_defaults(self):
        assert _extract_id({"id": "a"}) == "a"
        assert _extract_id({"target": "b"}) == "b"


class TestPhylogeneticHGTDetector:
    def test_detect_identical_trees(self):
        cfg = HGTConfig()
        det = PhylogeneticHGTDetector(cfg)
        tree = Tree(newick="((A,B),(C,D));")
        result = det.detect("m1", tree, tree)
        assert result.overall_risk == 0.0
        assert result.is_suspicious is False

    def test_detect_unavailable_rf(self):
        cfg = HGTConfig()
        det = PhylogeneticHGTDetector(cfg)
        gene = Tree(newick="((A,B),(C,D));")
        species = Tree(newick="(E,F);")
        result = det.detect("m1", gene, species)
        # No common tips -> quartet consistency 1.0 but RF unavailable -> conservative
        assert 0 <= result.overall_risk <= 1

    def test_quartet_consistency_sampling(self):
        cfg = HGTConfig()
        det = PhylogeneticHGTDetector(cfg)
        # 8 tips -> more than default sample of 1000 combos? math.comb(8,4)=70
        tips = [f"t{i}" for i in range(8)]
        gene = Tree(newick="(" + ",".join(tips) + ");")
        species = Tree(newick="(" + ",".join(tips) + ");")
        q, reason = det._calculate_quartet_consistency(gene, species)
        assert q == 1.0
        assert reason == ""

    def test_unrank_quartet(self):
        det = PhylogeneticHGTDetector(HGTConfig())
        tips = [f"t{i}" for i in range(10)]
        assert det._unrank_quartet(tips, 0) == ("t0", "t1", "t2", "t3")
        # Last combination
        import math
        last = math.comb(10, 4) - 1
        assert det._unrank_quartet(tips, last) == ("t6", "t7", "t8", "t9")


class TestHGTDecisionEngine:
    def test_far_distance_active(self):
        cfg = HGTConfig(adaptive_far_thresholds=True)
        cfg._far_distance_active = True
        engine = HGTDecisionEngine(cfg)
        result = engine.evaluate_marker("m1", PhyloStepResult(overall_risk=0.7))
        assert result.level == MarkerLevel.LEVEL_2

    def test_silent_far_log(self):
        cfg = HGTConfig(adaptive_far_thresholds=True)
        cfg._far_distance_active = True
        engine = HGTDecisionEngine(cfg)
        result = engine.evaluate_marker("m1", PhyloStepResult(overall_risk=0.7), silent=True)
        assert result.level == MarkerLevel.LEVEL_2

    def test_generate_empty_report(self):
        engine = HGTDecisionEngine(HGTConfig())
        report = engine.generate_hgt_report([])
        assert report.total_markers == 0
        assert report.mean_risk == 0.0


class TestHGTFilterModule:
    def test_run_with_species_tree(self, tmp_path):
        cfg = HGTConfig(enable_phylogenetic=True, phylogenetic_threshold=0.5)
        mod = HGTFilterModule(cfg, tmp_dir=str(tmp_path))
        species = Tree(newick="((g1,g2),(g3,g4));")
        seqs = {
            "M1": [{"id": "g1", "seq": "ACDEFG" * 10}, {"id": "g2", "seq": "ACDEFG" * 10},
                   {"id": "g3", "seq": "ACDEFG" * 10}, {"id": "g4", "seq": "ACDEFG" * 10}],
        }

        def runner(cmd, **kwargs):
            if cmd[0].lower() == "mafft":
                kwargs["stdout"].write(">g1\nACDEFGACDEFG\n>g2\nACDEFGACDEFG\n>g3\nACDEFGACDEFG\n>g4\nACDEFGACDEFG\n")
            elif cmd[0] == "FastTree" or cmd[0].lower() == "fasttree":
                out = cmd[cmd.index("-out") + 1]
                from pathlib import Path
                Path(out).write_text("((g1,g2),(g3,g4));")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            report, levels, _far = mod.run(["M1"], genomes=[], marker_sequences=seqs, species_tree=species)
        assert levels["M1"] == MarkerLevel.LEVEL_1
        assert report.total_markers == 1

    def test_run_with_taxonomy_table(self, tmp_path):
        cfg = HGTConfig(enable_phylogenetic=True, phylogenetic_threshold=0.5)
        mod = HGTFilterModule(cfg, tmp_dir=str(tmp_path))
        taxonomy = {
            "g1": {"domain": "B", "phylum": "P", "class": "C1"},
            "g2": {"domain": "B", "phylum": "P", "class": "C1"},
            "g3": {"domain": "B", "phylum": "P", "class": "C2"},
            "g4": {"domain": "B", "phylum": "P", "class": "C2"},
        }
        seqs = {
            "M1": [{"id": "g1", "seq": "ACDEFG" * 10}, {"id": "g2", "seq": "ACDEFG" * 10},
                   {"id": "g3", "seq": "ACDEFG" * 10}, {"id": "g4", "seq": "ACDEFG" * 10}],
        }

        def runner(cmd, **kwargs):
            if cmd[0].lower() == "mafft":
                kwargs["stdout"].write(">g1\nACDEFGACDEFG\n>g2\nACDEFGACDEFG\n>g3\nACDEFGACDEFG\n>g4\nACDEFGACDEFG\n")
            elif cmd[0] == "FastTree" or cmd[0].lower() == "fasttree":
                out = cmd[cmd.index("-out") + 1]
                from pathlib import Path
                Path(out).write_text("((g1,g2),(g3,g4));")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=runner):
            report, levels, _far = mod.run(
                ["M1"], genomes=[], marker_sequences=seqs,
                taxonomy_map=taxonomy, monophyly_rank="class",
            )
        try:
            import ete3  # Noqa: F401 (real probe: broken installs still have a spec)
            _have_ete3 = True
        except Exception:
            _have_ete3 = False
        if _have_ete3:
            assert levels["M1"] == MarkerLevel.LEVEL_1
        else:
            # Monophyly screen needs MAD rooting (ete3); without ete3
            # The marker is NOT_MEASURABLE and lands in UNKNOWN (documented).
            assert levels["M1"] == MarkerLevel.UNKNOWN

    def test_get_gene_tree_from_prebuilt(self, tmp_path):
        cfg = HGTConfig()
        mod = HGTFilterModule(cfg, tmp_dir=str(tmp_path))
        prebuilt = {"M1": str(tmp_path / "pre.nwk")}
        (tmp_path / "pre.nwk").write_text("((A,B),(C,D));")
        tree = mod._get_gene_tree("M1", [], prebuilt)
        assert tree is not None

    def test_build_gene_tree_too_few_seqs(self, tmp_path):
        cfg = HGTConfig()
        mod = HGTFilterModule(cfg, tmp_dir=str(tmp_path))
        assert mod._build_gene_tree("M1", [{"id": "a", "seq": "ACGT"}]) is None

    def test_build_gene_tree_uses_cache(self, tmp_path):
        """A cached tree is reused only when its input stamp matches.

        The stamp (``{marker}.nwk.input_sha256``) is what makes a cache hit mean
        "this tree was built from THESE sequences", so the case that reads the
        cache has to write the stamp the current input produces. It is derived
        from the same function the module uses, not typed out here.
        """
        from markerfinder.modules.hgt_filter import _marker_input_id

        cfg = HGTConfig()
        seqs = [{"id": f"s{i}", "seq": "ACDEFG" * 10} for i in range(4)]
        mod = HGTFilterModule(cfg, tmp_dir=str(tmp_path), gene_trees_dir=str(tmp_path))
        (tmp_path / "M1.nwk").write_text("((A,B),(C,D));")
        (tmp_path / "M1.nwk.input_sha256").write_text(_marker_input_id(seqs))
        tree = mod._build_gene_tree("M1", seqs)
        assert tree is not None

    def test_cache_is_keyed_by_its_input_stamp(self, tmp_path):
        """Control: an unstamped or wrongly stamped cache is never a hit."""
        from markerfinder.modules.hgt_filter import (
            _hgt_build_gene_tree_with_source,
            _marker_input_id,
        )

        seqs = [{"id": f"s{i}", "seq": "ACDEFG" * 10} for i in range(4)]
        (tmp_path / "M1.nwk").write_text("((A,B),(C,D));")

        # No stamp at all: the tree on disk says nothing about its input.
        _tree, source = _hgt_build_gene_tree_with_source(
            "M1", seqs, str(tmp_path), 1, str(tmp_path))
        assert source != "cached", source

        # A stamp from different sequences: same refusal.
        (tmp_path / "M1.nwk.input_sha256").write_text(
            _marker_input_id([{"id": "other", "seq": "ACDEFG" * 10}]))
        _tree, source = _hgt_build_gene_tree_with_source(
            "M1", seqs, str(tmp_path), 1, str(tmp_path))
        assert source != "cached", source
