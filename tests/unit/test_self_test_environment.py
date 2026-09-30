""": ``--check`` interpreter version gate must-fail control."""

import pytest

from markerfinder.cli.self_test import interpreter_in_range


class TestInterpreterRange:
    def test_in_range_lower_bound(self):
        assert interpreter_in_range((3, 10, 0)) is True

    def test_in_range_upper_bound(self):
        assert interpreter_in_range((3, 12, 7)) is True

    def test_below_range_fails(self):
        assert interpreter_in_range((3, 9, 18)) is False

    def test_out_of_range_upper_fails(self):
        # The dev machine's situation: 3.14 is outside >=3.10,<3.13.
        assert interpreter_in_range((3, 14, 6)) is False

    def test_major_jump_fails(self):
        assert interpreter_in_range((4, 0, 0)) is False

    def test_custom_bounds(self):
        assert interpreter_in_range((2, 7), min_version=(2, 7), max_exclusive=(3, 0))
        assert not interpreter_in_range((3, 0), min_version=(2, 7), max_exclusive=(3, 0))

    def test_test_dependencies_reports_version_fail_on_bad_interpreter(self, monkeypatch):
        """Must-fail: a fake out-of-range interpreter must turn the check red."""
        import markerfinder.cli.self_test as st

        monkeypatch.setattr(st.sys, "version_info", (3, 14, 6))
        results = dict((name, (status, detail)) for name, status, detail in st._test_dependencies())
        version_rows = [k for k in results if k.startswith("Interpreter 3.14")]
        assert version_rows, "interpreter version check entry missing"
        assert all(results[k][0] == "FAIL" for k in version_rows)
