""": ``--check`` interpreter version gate must-fail control."""

from markerfinder.cli.self_test import (
    SUPPORTED_PY_MAX_EXCLUSIVE,
    _declared_range_label,
    interpreter_in_range,
)


class TestInterpreterRange:
    def test_in_range_lower_bound(self):
        assert interpreter_in_range((3, 10, 0)) is True

    def test_in_range_upper_bound(self):
        assert interpreter_in_range((3, 12, 7)) is True

    def test_below_range_fails(self):
        assert interpreter_in_range((3, 9, 18)) is False

    def test_there_is_no_upper_bound(self):
        # 3.14 is supported: ``markerfinder._cgi_compat`` restored the stdlib
        # ``cgi`` module ete3 still imports, so the former ``<3.13`` ceiling was
        # retired and pyproject declares an unbounded ``>=3.10``. Anything at or
        # above the floor must therefore pass -- this is the assertion that stops
        # the retired ceiling from creeping back in.
        assert SUPPORTED_PY_MAX_EXCLUSIVE is None
        for version in ((3, 13, 0), (3, 14, 6), (3, 15, 0), (4, 0, 0)):
            assert interpreter_in_range(version) is True, version

    def test_declared_label_tracks_the_constants(self):
        # The label used to be a hardcoded ">=3.10,<3.13" that pyproject had
        # already stopped declaring. It must now be derived, not restated.
        assert _declared_range_label() == ">=3.10"

    def test_major_jump_passes_now_that_there_is_no_ceiling(self):
        # Documented consequence of the unbounded range, not an oversight.
        assert interpreter_in_range((4, 0, 0)) is True

    def test_custom_bounds(self):
        assert interpreter_in_range((2, 7), min_version=(2, 7), max_exclusive=(3, 0))
        assert not interpreter_in_range((3, 0), min_version=(2, 7), max_exclusive=(3, 0))

    def test_test_dependencies_reports_version_fail_on_bad_interpreter(self, monkeypatch):
        """Must-fail: a fake below-floor interpreter must turn the check red."""
        import markerfinder.cli.self_test as st

        monkeypatch.setattr(st.sys, "version_info", (3, 9, 18))
        results = dict((name, (status, detail)) for name, status, detail in st._test_dependencies())
        version_rows = [k for k in results if k.startswith("Interpreter 3.9")]
        assert version_rows, "interpreter version check entry missing"
        assert all(results[k][0] == "FAIL" for k in version_rows)
        # The failure text must name the range the manifest actually declares.
        assert ">=3.10" in results[version_rows[0]][1]
