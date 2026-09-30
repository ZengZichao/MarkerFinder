"""Strict static checking of the new/rewritten modules, as a gate.

The design scope limits mypy strict to the rewritten modules ("既有代码不追溯").
Without this test the check is never actually RUN — the acceptance item was
carried on the assumption that the annotations looked right. Running it found a
real latent crash: ``_stance`` computed ``agreement - norm_rf`` after checking
only ``rf is None``, so a provider returning ``(0, None)`` would raise TypeError
inside the run path that wiring had just made reachable.

The gate reports NOT EXECUTED rather than passing when mypy is not installed:
a skipped static check must never be read as a clean one.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

# The eight modules names, plus 's must-pass gate added later.
STRICT_MODULES = (
    "markerfinder/utils/etree.py",
    "markerfinder/assertions.py",
    "markerfinder/utils/reference.py",
    "markerfinder/models/evidence.py",
    "markerfinder/utils/informative_sites.py",
    "markerfinder/modules/hgt_scan.py",
    "markerfinder/modules/consistency_screen.py",
    "markerfinder/modules/composition.py",
    "markerfinder/modules/taxonomy_mustpass.py",
    # Growth-only addition (the anti-shrink lock stays).
    # Models/tree.py is where a None-vs-float comparison and ete3's default
    # Support label both survived every gate, precisely because it was not
    # Checked here. See test_models_tree_not_measurable.py.
    "markerfinder/models/tree.py",
    # Ledger: new module, so applies to it from day one.
    "markerfinder/utils/unmeasured.py",
    # Growth-only addition: the two modules every stage reads from and
    # Writes to. Both are plain dataclasses, so the cost is annotation-only, and
    # The payoff is that the input/output contract of all five stages is now
    # Type-checked: 6 bare ``List``/``Dict``/``tuple`` declarations were hiding
    # Here (untyped generics are an error under --strict, not a warning).
    "markerfinder/config.py",
    "markerfinder/models/pipeline_types.py",
)


def test_all_strict_modules_still_exist():
    # A renamed or deleted file must not silently shrink the checked set.
    missing = [m for m in STRICT_MODULES if not (REPO / m).exists()]
    assert not missing, f"strict-check list points at absent files: {missing}"


def test_mypy_strict_reports_no_errors():
    if importlib.util.find_spec("mypy") is None:
        pytest.skip(
            "NOT EXECUTED: mypy is not installed in this interpreter "
            "(pip install -e '.[dev]' / conda env). The gate is unproven here."
        )
    result = subprocess.run(
        [sys.executable, "-m", "mypy", "--strict",
         "--follow-imports=silent", "--ignore-missing-imports",
         *STRICT_MODULES],
        cwd=str(REPO), capture_output=True, text=True,
    )
    output = (result.stdout + result.stderr).strip()
    assert result.returncode == 0, (
        "mypy --strict failed on the checked module set:\n" + output[-4000:]
    )


def test_dev_extra_declares_mypy():
    """The gate is only real if CI can install it."""
    try:
        import tomllib
    except ModuleNotFoundError:  # Python 3.10, the declared floor
        import tomli as tomllib
    data = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    dev = " ".join(data["project"]["optional-dependencies"]["dev"]).lower()
    assert "mypy" in dev, "mypy missing from the dev extra"
    runtime = " ".join(data["project"]["dependencies"]).lower()
    # Mypy must stay out of the runtime dependency set.
    assert "mypy" not in runtime
