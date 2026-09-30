"""Unit tests for stepwise subcommand state inspection helpers."""

import sys
from io import StringIO
from pathlib import Path

import pytest

from markerfinder.__main__ import (
    _handle_step_state,
    _load_step_state,
    _step_status,
    _truncate_state_for_redo,
)
from markerfinder.models.report import PhaseContext
from markerfinder.utils.state_codec import save_pipeline_state


def _make_state(output_dir: Path, completed_steps):
    state = {
        "completed_steps": completed_steps,
        "context": PhaseContext(),
        "genomes": [],
        "marker_selection": None,
        "marker_sequences": {},
        "rank_map": {},
    }
    save_pipeline_state(str(output_dir), state)


class TestLoadStepState:
    def test_missing_returns_none(self, tmp_path):
        assert _load_step_state(str(tmp_path)) is None

    def test_loads_existing_state(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter"])
        state = _load_step_state(str(tmp_path))
        assert state is not None
        assert state["completed_steps"] == ["scan", "filter"]


class TestStepStatus:
    def test_no_state(self, tmp_path):
        status = _step_status("scan", str(tmp_path))
        assert not status["state_exists"]
        assert status["completed_steps"] == []

    def test_scan_not_completed(self, tmp_path):
        _make_state(tmp_path, [])
        status = _step_status("scan", str(tmp_path))
        assert status["state_exists"]
        assert not status["already_completed"]
        assert status["prerequisite_ok"]

    def test_scan_completed(self, tmp_path):
        _make_state(tmp_path, ["scan"])
        status = _step_status("scan", str(tmp_path))
        assert status["already_completed"]
        assert status["prerequisite_ok"]

    def test_filter_missing_scan(self, tmp_path):
        _make_state(tmp_path, [])
        status = _step_status("filter", str(tmp_path))
        assert not status["already_completed"]
        assert not status["prerequisite_ok"]
        assert status["missing_prerequisites"] == ["scan"]
        assert status["next_needed_step"] == "scan"

    def test_filter_ready(self, tmp_path):
        _make_state(tmp_path, ["scan"])
        status = _step_status("filter", str(tmp_path))
        assert not status["already_completed"]
        assert status["prerequisite_ok"]

    def test_filter_completed(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter"])
        status = _step_status("filter", str(tmp_path))
        assert status["already_completed"]
        assert status["prerequisite_ok"]


class TestTruncateStateForRedo:
    def test_truncate_filter(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter", "infer", "report"])
        _truncate_state_for_redo(str(tmp_path), "filter")
        state = _load_step_state(str(tmp_path))
        assert state["completed_steps"] == ["scan"]

    def test_truncate_scan(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter"])
        _truncate_state_for_redo(str(tmp_path), "scan")
        state = _load_step_state(str(tmp_path))
        assert state["completed_steps"] == []


class FakeTTY(StringIO):
    def isatty(self):
        return True


class TestHandleStepState:
    def _args(self, tmp_path, **kwargs):
        class Args:
            output = str(tmp_path)
            redo = False
            resume = False

        for k, v in kwargs.items():
            setattr(Args, k, v)
        return Args()

    def test_no_state_runs_requested(self, tmp_path):
        args = self._args(tmp_path)
        cmd, code = _handle_step_state("scan", args)
        assert cmd == "scan"
        assert code == 0

    def test_completed_step_non_tty_warns_and_exits(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter"])
        args = self._args(tmp_path)
        cmd, code = _handle_step_state("filter", args)
        assert cmd is None
        assert code == 0

    def test_completed_step_with_redo_runs(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter"])
        args = self._args(tmp_path, redo=True)
        cmd, code = _handle_step_state("filter", args)
        assert cmd == "filter"
        assert code == 0

    def test_missing_prerequisite_non_tty_warns_and_errors(self, tmp_path):
        _make_state(tmp_path, ["scan"])
        args = self._args(tmp_path)
        cmd, code = _handle_step_state("infer", args)
        assert cmd is None
        assert code == 2

    def test_missing_prerequisite_with_resume_runs(self, tmp_path):
        _make_state(tmp_path, ["scan"])
        args = self._args(tmp_path, resume=True)
        cmd, code = _handle_step_state("infer", args)
        assert cmd == "infer"
        assert code == 0

    def test_completed_step_tty_choose_redo(self, tmp_path, monkeypatch):
        _make_state(tmp_path, ["scan", "filter"])
        args = self._args(tmp_path)
        monkeypatch.setattr(sys, "stdin", FakeTTY("r\n"))
        cmd, code = _handle_step_state("filter", args)
        assert cmd == "filter"
        assert args.redo is True
        assert code == 0

    def test_completed_step_tty_choose_skip(self, tmp_path, monkeypatch):
        _make_state(tmp_path, ["scan", "filter"])
        args = self._args(tmp_path)
        monkeypatch.setattr(sys, "stdin", FakeTTY("s\n"))
        cmd, code = _handle_step_state("filter", args)
        assert cmd is None
        assert code == 0

    def test_missing_prerequisite_tty_choose_continue(self, tmp_path, monkeypatch):
        _make_state(tmp_path, ["scan"])
        args = self._args(tmp_path)
        monkeypatch.setattr(sys, "stdin", FakeTTY("c\n"))
        cmd, code = _handle_step_state("infer", args)
        assert cmd == "infer"
        assert args.resume is True
        assert code == 0

    def test_resume_from_completed_step_advances(self, tmp_path):
        _make_state(tmp_path, ["scan", "filter"])
        args = self._args(tmp_path, resume=True)
        cmd, code = _handle_step_state("filter", args)
        assert cmd == "infer"
        assert code == 0
