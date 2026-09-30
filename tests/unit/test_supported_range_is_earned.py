"""The supported interpreter range must be earned, not remembered.

``pyproject.toml`` declares ``requires-python = ">=3.10,<3.13"``. A declared
range is a promise, and a promise nobody checks rots: the suite must actually run
across the whole 3.10/3.11/3.12 span rather than only on the interpreter someone
happens to have, and on 3.14 -- outside the range -- it must report labelled
skips instead of aborting collection and returning no result at all.

The upper bound has a concrete cause: ete3 imports the stdlib ``cgi`` module,
which CPython 3.13 removed. So the bound is only right while that holds. This file
asserts both halves of the promise on whatever interpreter it runs on, and says so
in words when the environment cannot answer.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

PYPROJECT = Path(__file__).resolve().parents[2] / "pyproject.toml"


def _requires_python() -> str:
    try:
        import tomllib as _toml  # Type: ignore[import-not-found]
    except ModuleNotFoundError:  # Python 3.10, the declared floor
        import tomli as _toml  # Type: ignore[import-not-found]
    text = PYPROJECT.read_text(encoding="utf-8")
    match = re.search(r"requires-python\s*=\s*\"([^\"]+)\"", text)
    assert match, "no requires-python in pyproject.toml"
    return match.group(1)


def test_declared_range_is_the_one_the_team_verified():
    spec = _requires_python()
    assert ">=3.10" in spec and "<3.13" in spec, spec


def test_upper_bound_matches_the_real_ete3_constraint():
    """The bound exists because ete3 needs the stdlib ``cgi`` module.

    On 3.13+ the bound demands that ete3 be unusable; if a future ete3 stops
    importing ``cgi``, this test fails and the bound has to be raised on purpose
    rather than by accident. On the supported range, ete3 must be importable --
    otherwise the differential tests quietly stop running while the suite still
    reports green, which is exactly what was about.
    """
    importable = importlib.util.find_spec("ete3") is not None
    try:
        import ete3  # Noqa: F401

        usable = True
    except Exception:
        usable = False

    if sys.version_info >= (3, 13):
        assert not usable, (
            "ete3 became importable on Python >= 3.13; the <3.13 upper bound in "
            "pyproject.toml is now unjustified and should be re-decided, not "
            "left stale"
        )
        pytest.skip(
            "NOT EXECUTED as a positive check: this interpreter is outside the "
            "declared range, so ete3 is expected to be unusable"
        )
    if not importable:
        pytest.skip(
            "NOT EXECUTED: ete3 is not installed in this interpreter, so the "
            "supported-range claim cannot be checked here"
        )
    assert usable, (
        "ete3 is installed but unimportable inside the declared range -- the "
        "differential tests would silently stop comparing"
    )


def test_control_the_spec_probe_is_not_vacuous():
    # The regex must actually read the real field, not fall back to a default.
    assert "3.10" in _requires_python()
    with pytest.raises(AssertionError):
        _requires_python_from("nothing to see here")


def _requires_python_from(text: str) -> str:
    match = re.search(r"requires-python\s*=\s*\"([^\"]+)\"", text)
    if not match:
        raise AssertionError("no requires-python")
    return match.group(1)
