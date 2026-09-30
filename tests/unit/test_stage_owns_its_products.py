"""A stage computes AND writes its own artifact.

``markerfinder/pipeline.py`` is the coupling hotspot of the pipeline
(fan_out 24, ~1480 lines). It used to reach into stage internals and do their
disk I/O itself, which costs two things: a stage cannot be exercised end-to-end
without constructing ``MarkerFinderPipeline``, and the rules about *when* an
artifact exists live one layer above the module that defines them.

The two pilots chosen were ``composition`` and ``consistency_screen`` because
both already have standalone unit tests, so the move is verifiable without an
end-to-end run. The contract locked here is deliberately narrow:

  * the orchestrator must not call a stage's writer directly;
  * each stage's ``run_stage`` returns its path plus its own counts;
  * ``run_stage`` produces byte-identical output to the writer it wraps, so the
    shipped products are unchanged by the refactor;
  * an empty request is reported as "no artifact", not as an empty file.
"""
from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PIPELINE_SRC = (REPO / "markerfinder" / "pipeline.py").read_text(encoding="utf-8")

# Writers that were moved behind a stage entry point. Adding a new one here is
# The intended way to extend the policy (growth only, like STRICT_MODULES).
SUNK_WRITERS = (
    "write_composition_tsv(",
    "write_excluded_profile(",
    "write_mustpass_tsv(",
    "run_threshold_scan(",
)


@pytest.fixture()
def rows():
    return [
        {"marker_id": "M1", "rcv": 0.4, "gc_bias": 0.1, "n_sequences": 3},
        {"marker_id": "M2", "rcv": None, "gc_bias": 0.2, "n_sequences": 2},
    ]


def test_the_orchestrator_no_longer_calls_stage_writers_directly():
    offenders = [name for name in SUNK_WRITERS if name in PIPELINE_SRC]
    assert not offenders, (
        f"pipeline.py writes stage artifacts itself again: {offenders}; call the "
        "stage's run_stage() instead, or extend SUNK_WRITERS deliberately"
    )


def test_pipeline_goes_through_the_stage_entry_points():
    assert "run_composition_stage(" in PIPELINE_SRC
    assert "run_consistency_stage(" in PIPELINE_SRC
    assert "evaluate_and_write(" in PIPELINE_SRC
    assert "run_scan_stage(" in PIPELINE_SRC


def test_must_pass_stage_stamps_its_own_requirements_file(tmp_path):
    """The orchestrator used to have to remember to set the field between
    checking and writing. The stage now does it, so the record cannot be written
    without naming the file it was evaluated against."""
    from markerfinder.modules.taxonomy_mustpass import evaluate_and_write

    trees = {"M1": "((a,b),(c,d));"}
    taxonomy = {
        "a": {"domain": "Bacteria", "phylum": "P"}, "b": {"domain": "Bacteria", "phylum": "P"},
        "c": {"domain": "Bacteria", "phylum": "Q"}, "d": {"domain": "Bacteria", "phylum": "Q"},
    }
    requirements = [{"rank": "phylum", "child": "P", "parent": "domain:Bacteria"}]
    outcome = evaluate_and_write(
        trees, taxonomy, requirements,
        requirements_file=tmp_path / "baseline.yaml",
        output_dir=str(tmp_path), prefix="mf",
    )
    assert outcome.report.requirements_file == str(tmp_path / "baseline.yaml")
    assert outcome.path.name == "mf.taxonomy_mustpass.tsv" or outcome.path.exists()
    assert outcome.path.read_text(encoding="utf-8")


def test_threshold_scan_stage_names_its_artifact(tmp_path):
    from markerfinder.modules.hgt_scan import product_path, run_stage
    from markerfinder.models.marker import MarkerLevel
    from markerfinder.models.pipeline_types import HGTEvaluation

    evaluations = [
        HGTEvaluation(marker_id=f"M{i}", overall_risk=risk, level=MarkerLevel.LEVEL_1)
        for i, risk in enumerate((0.05, 0.3, 0.7))
    ]
    outcome = run_stage(evaluations, str(tmp_path), "mf")
    assert outcome.path == product_path(str(tmp_path), "mf")
    assert outcome.path.name == "mf.threshold_scan.tsv"
    assert len(outcome.results) >= 1
    assert outcome.path.read_text(encoding="utf-8").startswith("band_level1_max")


def test_composition_stage_writes_byte_identical_output(rows, tmp_path):
    from markerfinder.modules.composition import run_stage, write_composition_tsv

    direct = write_composition_tsv(rows, str(tmp_path / "a"), "mf", warn_threshold=0.15)
    outcome = run_stage(
        rows, str(tmp_path / "b"), "mf", warn_threshold=0.15
    )
    assert outcome.path is not None and outcome.rows_written == 2
    assert outcome.path.read_text(encoding="utf-8") == direct.read_text(encoding="utf-8")
    assert outcome.path.name == "mf.composition.tsv"


def test_composition_stage_reports_no_artifact_for_an_empty_request(tmp_path):
    outcome = __import__(
        "markerfinder.modules.composition", fromlist=["run_stage"]
    ).run_stage([], str(tmp_path), "mf", warn_threshold=None)
    assert outcome.path is None, "an empty screen must not leave a header-only file"
    assert outcome.rows_written == 0
    assert not (tmp_path / "Phase5_reports").exists()


def _mk(cs, marker_id: str, grade):
    """One screened marker; the two stances are irrelevant to grade tallying."""
    return cs.MarkerConsistency(
        marker_id=marker_id,
        grade=grade,
        concat_stance=None,
        coalescent_stance=None,
        stringency=1,
    )


def test_consistency_stage_counts_its_own_profile(tmp_path, monkeypatch):
    import markerfinder.modules.consistency_screen as cs

    def _fake_writer(results, output_dir, prefix, *, pis_map=None, category_map=None):
        path = Path(output_dir) / "profile.tsv"
        path.write_text("stub\n", encoding="utf-8")
        return path

    monkeypatch.setattr(cs, "write_excluded_profile", _fake_writer)
    consistent = _mk(cs, "A", cs.ConsistencyGrade.CONSISTENT)
    dropped = _mk(cs, "B", cs.ConsistencyGrade.INCONSISTENT)
    unknown = _mk(cs, "C", cs.ConsistencyGrade.INCONCLUSIVE)

    outcome = cs.run_stage([consistent, dropped, unknown], str(tmp_path), "mf")

    assert outcome.profile.rows_written == 2, "only non-CONSISTENT markers are profiled"
    assert outcome.grades == {
        "consistent": 1, "inconsistent": 1, "inconclusive": 1,
    }
    assert outcome.profile.path.name == "profile.tsv"


def test_summarize_grades_always_carries_all_three_keys():
    import markerfinder.modules.consistency_screen as cs

    assert cs.summarize_grades([]) == {
        "consistent": 0, "inconsistent": 0, "inconclusive": 0,
    }
    only = cs.summarize_grades([_mk(cs, "A", cs.ConsistencyGrade.CONSISTENT)])
    assert only["consistent"] == 1 and only["inconsistent"] == 0


def test_control_the_writer_scan_is_not_vacuous():
    """The lock in this file is a source scan; prove the scan can see a call."""
    planted = "def run_infer(self):\n    return write_composition_tsv(rows, out, pfx)\n"
    assert any(name in planted for name in SUNK_WRITERS)
