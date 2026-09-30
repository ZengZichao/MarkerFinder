"""Parsimony-informative sites tests."""

import pytest

from markerfinder.utils.informative_sites import (
    effective_columns,
    gap_fraction,
    parsimony_informative_sites,
)
from markerfinder.models.marker import MarkerQualityScore

# Hand-computed reference alignment ( ②):
# Col1 A,A,G,G -> A2 G2 informative
# Col2 G,G,A,A -> G2 A2 informative
# Col3 G,G,G,G -> single residue, not informative
# Col4 C,C,A,A -> C2 A2 informative
# Col5 T,T,T,T -> not informative
# Col6 ---- all-gap, excluded
ALN = {
    "s1": "AGGCT-",
    "s2": "AGGCT-",
    "s3": "GAGAT-",
    "s4": "GAGAT-",
}


class TestPis:
    def test_known_alignment_pis(self):
        assert parsimony_informative_sites(ALN) == 3

    def test_constant_sites_not_counted(self):
        aln = {"a": "AAA", "b": "AAA", "c": "AAA", "d": "AAA"}
        assert parsimony_informative_sites(aln) == 0

    def test_singleton_variants_not_counted(self):
        # One divergent singleton does not make a column informative
        aln = {"a": "AAAA", "b": "AAAA", "c": "AAAA", "d": "GAAA"}
        assert parsimony_informative_sites(aln) == 0

    def test_two_pairs_is_informative(self):
        # All four columns hold two residues x2 -> every column informative
        aln = {"a": "GGTT", "b": "GGTT", "c": "TTGG", "d": "TTGG"}
        assert parsimony_informative_sites(aln) == 4

    def test_all_gap_column_excluded(self):
        aln = {"a": "A---", "b": "A---", "c": "G---", "d": "G---"}
        # Col1 informative, the gap columns contribute nothing
        assert parsimony_informative_sites(aln) == 1
        assert effective_columns(aln) == 1

    def test_ambiguous_treated_as_gap(self):
        aln = {"a": "AGXX", "b": "AGXX", "c": "TGXX", "d": "TGXX"}
        assert parsimony_informative_sites(aln) == 1


class TestEffectiveColumnsAndGaps:
    def test_gappy_alignment_effective_columns(self):
        assert effective_columns(ALN) == 5
        assert gap_fraction(ALN) == pytest.approx(4 / 24)

    def test_empty_alignment(self):
        assert parsimony_informative_sites({}) == 0
        assert effective_columns({}) == 0
        assert gap_fraction({}) == 0.0

    def test_uneven_lengths_raise(self):
        with pytest.raises(ValueError):
            parsimony_informative_sites({"a": "AAA", "b": "AA"})


class TestSummaryColumns:
    def test_marker_summary_has_pis_column(self, tmp_path):
        """Pis + effective_columns appended to marker_summary.tsv."""
        from markerfinder.config import ReportConfig
        from markerfinder.models.marker import (
            MarkerLevel,
            MarkerQualityScore,
            SelectedMarkerSet,
        )
        from markerfinder.models.pipeline_types import MarkerSelectionResult
        from markerfinder.modules.report_generator import PlainTextReportGenerator

        config = ReportConfig(output_dir=str(tmp_path), output_prefix="mf")
        generator = PlainTextReportGenerator(config)

        qs = MarkerQualityScore(
            marker_id="M1", occupancy=0.9,
            phylogenetic_informativeness=0.8, hgt_risk_score=0.0,
            pis=42, effective_columns=150,
        )
        marker_set = SelectedMarkerSet(
            markers=["M1"], occupancy_scores={"M1": 0.9},
            quality_scores={"M1": qs},
        )
        marker_result = MarkerSelectionResult(
            marker_set=marker_set, quality_scores={"M1": qs},
        )
        reports_dir = tmp_path / "Phase5_reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        path = generator._write_marker_summary(marker_result, reports_dir, "mf")
        content = path if isinstance(path, str) else str(path)
        lines = open(content, encoding="utf-8").read().splitlines()
        header = lines[0].split("\t")
        assert header[:4] == [
            "marker_id", "occupancy_score", "marker_quality_score",
            "marker_quality_level",
        ]
        assert header[4:] == ["pis", "effective_columns", "consistency_grade"]
        # Appended right; consistency_grade is NA because this criterion only
        # Runs under --hgt-mode consistency, and the default must not borrow a
        # Grade from the risk path.
        assert lines[1].split("\t")[4:] == ["42", "150", "NA"]

    def test_unmeasured_pis_renders_na(self):
        qs = MarkerQualityScore(marker_id="M1")
        assert qs.pis is None
        from markerfinder.modules.report_generator import _plain
        assert _plain(qs.pis) == "NA"

    def test_protein_gc_not_applicable_note(self):
        # 边界：PIS 是残基信息量；GC 指标属 且蛋白输入显式 N/A，
        # 此处锁住 pis 字段不承载 GC 语义。
        from dataclasses import fields
        names = {f.name for f in fields(MarkerQualityScore)}
        assert "pis" in names and "effective_columns" in names
        assert "gc_content" not in names
