"""The supported interpreter range must be earned, not remembered.

``pyproject.toml`` declares ``requires-python = ">=3.10"``. A declared
range is a promise, and a promise nobody checks rots: the suite must actually run
across the whole span rather than only on the interpreter someone
happens to have.

The range used to be ``>=3.10,<3.13``. That upper bound had a concrete cause:
ete3 imports the stdlib ``cgi`` module, which CPython 3.13 removed (PEP 594),
which made ete3 — and with it every MAD / monophyly measurement — silently
unusable on newer interpreters. ``markerfinder._cgi_compat`` now installs a
minimal ``cgi`` stand-in, so ete3 is importable on 3.13+ and there is no reason to
exclude those interpreters. This file asserts the promise that replaced the
retired one: ete3 must be importable on **every** supported interpreter, because
the pure-fallback path must never be reached merely because the stdlib moved on.
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
    assert ">=3.10" in spec, spec
    assert "<3.13" not in spec, (
        "the <3.13 upper bound was retired once _cgi_compat made ete3 usable on "
        f"3.13+; found a stale bound in {spec!r}"
    )


def test_ete3_is_usable_on_every_supported_interpreter():
    """ete3 must import on every interpreter inside the declared range.

    The pure-Python split-set fallback exists for trees too large for ete3, not
    as a way to absorb a broken import. If ete3 stops importing, the differential
    tests quietly stop comparing while the suite still reports green -- the
    failure mode this assertion exists to prevent.

    On 3.13+ this only holds because ``markerfinder._cgi_compat`` restores the
    removed ``cgi`` module, hence the deliberate package import first.

    The failure text carries the underlying exception. It used to swallow it,
    which cost real diagnosis time: an undeclared transitive dependency of ete3
    (it does ``import six.moves.cPickle`` while declaring no requirements at
    all) made every ete3-dependent test fail in CI with a message that named
    only the symptom, "ete3 is installed but unimportable". A green local
    checkout hid it, because some unrelated package happened to have six
    installed. An assertion that cannot say why cannot be acted on.
    """
    importable = importlib.util.find_spec("ete3") is not None
    failure = ""
    try:
        import markerfinder  # noqa: F401  -- installs the cgi shim when needed

        import ete3  # Noqa: F401

        usable = True
    except Exception as exc:
        usable = False
        failure = f"{type(exc).__name__}: {exc}"

    if not importable:
        pytest.skip(
            "NOT EXECUTED: ete3 is not installed in this interpreter, so the "
            "supported-range claim cannot be checked here"
        )
    assert usable, (
        "ete3 is installed but unimportable on Python "
        f"{sys.version_info.major}.{sys.version_info.minor}, which is inside the "
        "declared range -- the differential tests would silently stop comparing. "
        f"Underlying import error: {failure or 'none captured'}"
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
