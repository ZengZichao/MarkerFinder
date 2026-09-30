"""Exit codes and shared CLI constants for the MarkerFinder command line.

These mirror the historical constants defined in ``markerfinder.__main__`` so
that the split ``cli`` sub-package can reference them without creating a
circular import back into ``__main__``.
"""

from __future__ import annotations

from pathlib import Path

# Where the shipped must-pass baseline lives inside a source checkout.
_MUSTPASS_RELPATH = "tests/benchmark/expected/taxonomy_mustpass.yaml"


def recommendation_is_inconclusive(recommendation) -> bool:
    """Did the recommendation layer finish without being able to name a tree?

    Exit code 5 (``EXIT_INCONCLUSIVE``, / ) exists for exactly this
    state: the pipeline ran to completion, conflict metrics were unavailable,
    and ``recommend_tree`` returned ``recommended_tree=None`` /
    ``confidence="inconclusive"``. Both invocation routes — the one-shot run and
    the ``infer``/``report`` subcommands — must answer with the same code, or the
    same scientific state exits 0 on one path and 5 on the other.
    """
    if recommendation is None:
        return False
    if str(getattr(recommendation, "confidence", "")) == "inconclusive":
        return True
    return bool(getattr(recommendation, "recommended_tree", None) is None
                 and getattr(recommendation, "reason", None))


def default_taxonomy_mustpass_path() -> str:
    """Resolve the bare ``--taxonomy-mustpass`` default to an absolute path.

    A repo-relative default would only work when the user happens to run from
    the checkout root, and the gate would then abort for a packaging reason
    rather than a biological one.
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / _MUSTPASS_RELPATH
        if candidate.exists():
            return str(candidate)
    return _MUSTPASS_RELPATH


EXIT_SUCCESS = 0
EXIT_RUNTIME_ERROR = 1
EXIT_ARG_ERROR = 2
EXIT_DATA_ERROR = 3
# Output-numeric assertion failures
# And taxonomy must-pass failures exit with this code so
# "output unreasonable" is distinguishable from a runtime error.
EXIT_ASSERTION_FAILED = 4
# The pipeline completed but the recommendation layer
# Could not hand out a tree (conflict metrics unmeasurable).
EXIT_INCONCLUSIVE = 5
EXIT_INTERRUPT = 130
