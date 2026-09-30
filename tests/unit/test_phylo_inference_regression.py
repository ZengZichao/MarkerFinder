"""Code-review regression tests for ``modules/phylogenetic_inference.py`` (#5 C29 + D38).

Locks in the ``recommend_tree`` changes:
  - C29/D38: ``_metrics_str`` folds normalized RF + quartet consistency into
             the recommendation reason string.
  - #5: an ``rf is None`` (empty/degenerate trees where RF is uncomputable)
        takes the dedicated fallback branch — no crash, and a reason string is
        produced via the None fallback.
"""

import pytest

from markerfinder.config import PhylogeneticConfig
from markerfinder.models.tree import Tree
from markerfinder.models.pipeline_types import (
    ConflictReport,
    CoalescentResult,
    PhylogeneticResult,
    SupermatrixResult,
)
from markerfinder.modules.phylogenetic_inference import (
    PhylogeneticInferenceModule,
)


class TestRecommendTreeRfNone:
    def test_rf_none_fallback_no_crash(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig())
        pr = PhylogeneticResult(
            supermatrix=SupermatrixResult(tree=Tree(newick="(A,B);")),
            coalescent=CoalescentResult(
                species_tree=Tree(newick="(A,B);"),
                species_tree_source="astral",
            ),
            conflict_report=ConflictReport(
                normalized_rf=None, quartet_agreement=None
            ),
        )
        rec = mod.recommend_tree(pr)
        # Unmeasurable conflict metrics => inconclusive and
        # NO tree handed out (previously fell back to the coalescent tree
        # With confidence "unknown", hiding the unmeasurability).
        assert rec.confidence == "inconclusive"
        assert rec.recommended_tree is None
        assert "unavailable" in rec.reason.lower()

    def test_reason_includes_metrics_when_available(self):
        mod = PhylogeneticInferenceModule(PhylogeneticConfig())
        pr = PhylogeneticResult(
            supermatrix=SupermatrixResult(tree=Tree(newick="(A,B);")),
            coalescent=CoalescentResult(
                species_tree=Tree(newick="(A,B);"),
                species_tree_source="astral",
            ),
            conflict_report=ConflictReport(
                normalized_rf=0.05, quartet_agreement=0.9
            ),
        )
        rec = mod.recommend_tree(pr)
        assert "RF=" in rec.reason
        assert rec.confidence == "high"
