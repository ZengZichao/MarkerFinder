"""Phase numbering unification tests."""

import re
from pathlib import Path

import pytest

from markerfinder.phases import HISTORICAL_ALIASES, PHASES, canonical_log_tag, get_phase


class TestPhasesModule:
    def test_all_canonical_phases_registered(self):
        numbers = [p.number for p in PHASES]
        assert numbers == ["0", "0.1", "0.2", "1", "1.5", "2", "3", "4"]

    def test_phase_1_5_records_unwired_status(self):
        p15 = get_phase("1.5")
        assert "NOT wired in" in p15.title

    def test_historical_aliases_cover_directories(self):
        # Directory names are unchanged but registered as aliases.
        for alias in ("Phase4_trees", "Phase5_reports", "Phase5_metadata"):
            assert alias in HISTORICAL_ALIASES

    def test_canonical_log_tag(self):
        assert canonical_log_tag("0.1") == "[Phase 0.1]"
        with pytest.raises(KeyError):
            canonical_log_tag("9")


class TestNoStalePhaseReferences:
    def test_pipeline_comments_no_longer_call_report_phase5(self):
        pipeline = Path(__file__).parent.parent.parent / "markerfinder" / "pipeline.py"
        src = pipeline.read_text(encoding="utf-8")
        assert "Step 4: report (Phase 5)" not in src
        assert "Step 4: report (Phase 4)" not in src  # Unified wording

    def test_hgt_filter_comment_names_canonical_phase(self):
        hgt = Path(__file__).parent.parent.parent / "markerfinder" / "modules" / "hgt_filter.py"
        src = hgt.read_text(encoding="utf-8")
        assert "供 Phase 4 (CoalescentInference)" not in src
        assert "canonical phase 3" in src

    def test_directory_names_unchanged(self):
        """The integration test asserting Phase5_reports must still hold."""
        integration = (
            Path(__file__).parent.parent / "integration" / "test_pipeline.py"
        )
        src = integration.read_text(encoding="utf-8")
        assert "Phase5_reports" in src  # Untouched by
