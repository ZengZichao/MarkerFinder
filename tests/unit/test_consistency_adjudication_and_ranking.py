"""'s grading table and 's demotion of the combined score.

 ("分档规则照搬论文：两侧同侧且都过阈 -> consistent；两侧异侧 ->
inconsistent；任一侧不过阈 -> inconclusive，不是'平均一下'") and
("``combine_hgt_scores`` 在 consistency 模式下降格为'排序用'而非'裁决用'") were both
implemented in prose and code but named by no test. Together they are the reason
consistency mode exists at all: if a graded marker can be produced by averaging two
signals, the criterion has quietly become the risk score again.

Every expectation below was measured first. Two of my first guesses were wrong:
with four tips, any pair of *different* 4-taxon topologies sits at maximum RF, so a
"consistent" grade needs both references to agree -- a fixture of "gene tree matches
concat, references differ" is INCONSISTENT, not CONSISTENT, and the naive table test
would have been asserting the opposite of the rule.

Interpretation on record (flagged for the author in the closure report): when one
side clears its threshold and the other does not, the two clauses of overlap
("两侧异侧 -> inconsistent" vs "任一侧不过阈 -> inconclusive"). The implementation
grades the split case as INCONSISTENT and reserves INCONCLUSIVE for "neither side
clears" or "a side could not be measured at all".
"""

from __future__ import annotations

from pathlib import Path

import pytest

from markerfinder.config import ReportConfig
from markerfinder.models.marker import MarkerLevel, SelectedMarkerSet
from markerfinder.models.pipeline_types import (
    ConflictReport,
    HGTEvaluation,
    HGTReport,
    MarkerSelectionResult,
    PhylogeneticResult,
    SupermatrixResult,
)
from markerfinder.models.tree import Tree
from markerfinder.modules.consistency_screen import ConsistencyGrade, screen
from markerfinder.modules.report_generator import PlainTextReportGenerator

REF = "((A:1,B:1):1,(C:1,D:1):1);"
OTHER = "((A:1,C:1):1,(B:1,D:1):1);"
TINY = "(A:1,B:1);"


class _Runtime:
    duration = 1.0


def _grades(genes, concat, astral, stringency=3):
    return {
        result.marker_id: result
        for result in screen(genes, concat, astral, stringency=stringency)
    }


# ──: the four rows of the table ─────────────────────────────────────

def test_control_both_sides_clearing_gives_consistent():
    result = _grades({"m": REF}, REF, REF)["m"]
    assert result.grade is ConsistencyGrade.CONSISTENT
    assert result.reasons == []


def test_control_neither_side_clearing_is_inconclusive_not_an_average():
    """The '不是平均一下' clause: a weak-on-both-legs marker is undecided, never a
    mid-range grade produced by blending two signals."""
    result = _grades({"m": OTHER}, REF, REF)["m"]
    assert result.grade is ConsistencyGrade.INCONCLUSIVE
    assert "neither side above threshold" in result.reasons


def test_control_opposite_sides_is_inconsistent():
    """References that disagree with each other cannot both be cleared by one
    marker tree: the split is the definition of inconsistent here."""
    for gene, marker in ((REF, "matches_concat"), (OTHER, "matches_astral")):
        result = _grades({marker: gene}, REF, OTHER)[marker]
        assert result.grade is ConsistencyGrade.INCONSISTENT, marker
        assert "sides disagree" in result.reasons, marker


def test_control_unmeasurable_side_is_inconclusive():
    result = _grades({"m": TINY}, REF, REF)["m"]
    assert result.grade is ConsistencyGrade.INCONCLUSIVE
    assert any("either side unmeasurable" == reason for reason in result.reasons)
    assert any("fewer than 4 shared tips" in r for r in result.reasons)


# ──: the risk score is ranking-only, and the products say so ────────

def _summary(tmp_path, mode: str, grades: dict) -> str:
    marker = "M1"
    hgt_report = HGTReport(
        marker_evaluations=[
            HGTEvaluation(marker_id=marker, overall_risk=0.1,
                          level=MarkerLevel.LEVEL_1)
        ],
        total_markers=1,
        evidence_coverage=1.0,
    )
    phylo = PhylogeneticResult(
        supermatrix=SupermatrixResult(tree=Tree(newick=REF)),
        consistency_grades=grades,
    )
    generator = PlainTextReportGenerator(
        ReportConfig(output_dir=str(tmp_path), output_prefix="mf", hgt_mode=mode)
    )
    out_dir = tmp_path / "Phase5_reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    generator._write_pipeline_summary(
        MarkerSelectionResult(
            marker_set=SelectedMarkerSet(
                markers=[marker], occupancy_scores={marker: 0.9}, quality_scores={}
            ),
            quality_scores={},
        ),
        hgt_report,
        phylo,
        _Runtime(),
        out_dir,
        "mf",
    )
    return (out_dir / "mf.pipeline_summary.txt").read_text(encoding="utf-8")


def test_consistency_mode_declares_the_demotion(tmp_path):
    text = _summary(tmp_path, "consistency", {"M1": "consistent"})
    assert "RANKING ONLY" in text
    assert "RANKING ONLY" in text
    assert "1 marker(s) graded" in text


def test_control_risk_mode_invents_no_demotion(tmp_path):
    text = _summary(tmp_path, "risk", {})
    assert "RANKING ONLY" not in text
    assert "RANKING ONLY" not in text, (
        "the default shipped mode must not gain a sentence it has not earned"
    )


def test_consistency_mode_without_any_grade_says_the_demotion_is_off(tmp_path):
    """Loud-failure half of: if the criterion could not run, the reader must
    learn that levels are risk-based after all, not be told a comforting lie."""
    text = _summary(tmp_path, "consistency", {})
    assert "consistency demotion NOT IN EFFECT" in text
    assert "fell back to risk grading" in text


def test_every_evaluation_carries_the_ranking_only_note():
    from markerfinder.pipeline import _annotate_risk_role

    evaluations = [
        HGTEvaluation(marker_id="m1", overall_risk=0.2, level=MarkerLevel.LEVEL_2),
        HGTEvaluation(
            marker_id="m2", overall_risk=0.4, level=MarkerLevel.UNKNOWN,
            notes="already annotated",
        ),
    ]
    report = HGTReport(marker_evaluations=evaluations, total_markers=2)
    _annotate_risk_role(report, "fr59: adjudication = consistency grade")
    assert evaluations[0].notes == "fr59: adjudication = consistency grade"
    assert evaluations[1].notes.startswith("already; ") or evaluations[1].notes == (
        "already annotated; fr59: adjudication = consistency grade"
    ), evaluations[1].notes
    for evaluation in evaluations:
        assert "fr59:" in evaluation.notes


def test_report_config_carries_the_mode_from_the_cli():
    from markerfinder.cli.config_build import _build_pipeline_config
    from markerfinder.cli.parser import _build_parser

    parser = _build_parser()
    inbox = Path("")
    args = parser.parse_args(
        ["-i", str(inbox.resolve()), "--require-evidence-coverage", "0.5"]
    )
    assert _build_pipeline_config(args).report_config.hgt_mode == "risk"

# ── the excluded profile must show the legs it graded on ──────────────────

def test_profile_carries_both_leg_strengths_and_their_observation_base(tmp_path):
    """The grade alone cannot be audited from the product.

    ``concat_stance`` / ``coalescent_stance`` used to be computed, stored on the
    result, and then dropped -- ``excluded_profile.tsv`` reported only the grade,
    so a reader could not see which leg said what, nor that an agreement over one
    quartet is much weaker evidence than one over two hundred.
    """
    from markerfinder.modules.consistency_screen import write_excluded_profile

    results = screen(
        {"m_both": OTHER, "m_tiny": TINY}, REF, REF, stringency=3
    )
    path = write_excluded_profile(results, str(tmp_path), "mf")
    lines = path.read_text(encoding="utf-8").strip().split("\n")
    header = lines[0].split("\t")
    assert header[-4:] == [
        "concat_strength", "concat_quartets",
        "coalescent_strength", "coalescent_quartets",
    ], header
    rows = {line.split("\t")[0]: line.split("\t") for line in lines[1:]}
    both = rows["m_both"]
    # Graded on measurable legs -> real numbers, and one quartet for a 4-taxon set
    assert float(both[-4]) == pytest.approx(-1.0), both
    assert both[-3] == "1", both
    assert float(both[-2]) == pytest.approx(-1.0), both
    tiny = rows["m_tiny"]
    for cell in tiny[-4:]:
        assert cell == "NA", (tiny, cell)
