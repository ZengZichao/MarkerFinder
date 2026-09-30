"""Per-marker decision card tests."""

import json

import pytest

from markerfinder.config import HGTConfig
from markerfinder.models.evidence import MeasureState
from markerfinder.models.marker import MarkerLevel
from markerfinder.models.pipeline_types import HGTEvaluation, PhyloStepResult
from markerfinder.models.tree import Tree
from markerfinder.modules.hgt_filter import HGTDecisionEngine, PhylogeneticHGTDetector
from markerfinder.utils.hgt_utils import combine_hgt_scores

CARD_KEYS = {
    "schema_version", "marker_id", "detector", "gene_tree_source",
    "trimming_regime", "rf", "quartet", "monophyly",
    "rank_used", "n_total", "n_mono", "weights", "n_signals_used",
    "thresholds_in_effect", "far_active", "risk_basis",
    "overall_risk", "level", "confidence", "assertion_ids_fired", "notes",
}


class TestDecisionCard:
    # The complete field list;: the same card in both carriers
    # (hgt_evaluation.tsv column + Phase5_evidence/decision_<n>.json).
    def test_card_contains_all_required_keys(self):
        engine = HGTDecisionEngine(HGTConfig())
        res = PhyloStepResult(overall_risk=0.1)
        ev = engine.evaluate_marker("m1", res)
        assert CARD_KEYS <= set(ev.decision_card.keys())
        assert ev.decision_card["schema_version"] == 2

    def test_unknown_card_has_reason(self):
        engine = HGTDecisionEngine(HGTConfig())
        res = PhyloStepResult(gene_id="m")
        res.rf_state = MeasureState.NOT_MEASURABLE
        ev = engine.evaluate_marker("m", res)
        assert ev.decision_card["detector"] == "none"
        assert ev.decision_card["unknown_reason"]

    def test_card_json_serializable(self):
        engine = HGTDecisionEngine(HGTConfig())
        ev = engine.evaluate_marker("m1", PhyloStepResult(overall_risk=0.4))
        blob = json.dumps(ev.decision_card, ensure_ascii=False, default=str)
        assert "m1" in blob

    def test_weights_and_thresholds_in_card(self):
        engine = HGTDecisionEngine(HGTConfig())
        ev = engine.evaluate_marker("m1", PhyloStepResult(overall_risk=0.4))
        assert ev.decision_card["weights"] == {"rf": 0.5, "quartet": 0.5}
        assert ev.decision_card["thresholds_in_effect"] == {
            "level1_max": 0.25, "level2_max": 0.60,
        }


class TestNotesAndTsv:
    def test_non_unknown_notes_never_empty(self):
        engine = HGTDecisionEngine(HGTConfig())
        for risk in (0.0, 0.1, 0.3, 0.4, 0.9):
            ev = engine.evaluate_marker("m", PhyloStepResult(overall_risk=risk))
            assert ev.notes.strip(), f"notes empty for risk={risk}"

    def test_tab_newline_sanitized_in_tsv_row(self):
        # 375 precedent — free text must not break the TSV shape.
        engine = HGTDecisionEngine(HGTConfig())
        ev = engine.evaluate_marker("m", PhyloStepResult(overall_risk=0.4))
        assert "\t" not in ev.notes and "\n" not in ev.notes


class TestStateSchema:
    def test_schema_version_is_2(self):
        from markerfinder.utils.state_codec import SCHEMA_VERSION
        assert SCHEMA_VERSION == 2

    def test_resume_with_old_state_gives_actionable_error(self, tmp_path):
        from markerfinder.utils.state_codec import StateSchemaError
        from markerfinder.utils.state_codec import (
            STATE_DIR_NAME, STATE_FILE_NAME, read_completed_steps,
        )

        state_dir = tmp_path / STATE_DIR_NAME
        state_dir.mkdir()
        (state_dir / STATE_FILE_NAME).write_text(
            json.dumps({"schema_version": 1, "completed_steps": ["scan"]}),
            encoding="utf-8",
        )
        with pytest.raises(StateSchemaError) as excinfo:
            read_completed_steps(str(tmp_path))
        msg = str(excinfo.value)
        assert "schema_version" in msg
        assert "--redo" in msg  # Actionable rerun hint

    def test_state_mismatch_not_silent(self):
        # Version mismatch must raise, never silently return None.
        import inspect
        from markerfinder.utils import state_codec

        source = inspect.getsource(state_codec.read_completed_steps)
        assert "return None  # no state" not in source
        assert "StateSchemaError" in source


class TestCombineStillReporting:
    def test_n_used_counted(self):
        _, n = combine_hgt_scores({"rf": 1.0, "quartet": 0.0},
                                  {"rf": 0.5, "quartet": 0.5})
        assert n == 2
