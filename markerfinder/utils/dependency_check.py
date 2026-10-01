"""Dependency pre-checks for MarkerFinder.

Called at pipeline initialization to detect missing critical dependencies
early, before any computation begins. This prevents confusing downstream
failures (e.g. ete3 ImportError inside an HGT filter worker thread).
"""

from __future__ import annotations

import logging
import shutil
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


# External CLI tools and their purpose. Each entry is
# (binary_name, purpose, is_critical).
# Critical tools raise an error if missing; non-critical produce a warning.
#
# The table used to stop at the six tools the original design assumed. The
# Package really invokes more of them -- ``mag_optimization._run_checkm`` runs
# ``checkm lineage_wf`` (mag_optimization.py:121-133) and silently falls back to
# Default quality estimates when CheckM is absent, so the run "succeeded" with
# Estimated completeness/contamination numbers and nothing in the startup
# Pre-check ever said so. Every entry below is backed by a subprocess call site.
EXTERNAL_TOOLS: List[Tuple[str, str, bool]] = [
    ("hmmsearch", "HMM marker scanning (Phase 1)", True),
    ("mafft", "multiple sequence alignment (Phase 1/2/3)", False),
    ("trimal", "alignment trimming (Phase 1/3)", False),
    ("fasttree", "fast gene-tree inference (Phase 1/2/3)", False),
    ("iqtree3", "ML gene-tree inference for GTDB-TK markers (Phase 1, "
                "optional; --gene-tree-builder iqtree)", False),
    ("astral", "coalescent species-tree inference (Phase 3, optional)", False),
    ("checkm", "MAG completeness/contamination (Phase 0, optional; missing "
               "falls back to default quality estimates)", False),
]

# Which pipeline steps invoke which binary, read off the call sites.
# ``check_external_tools(required_phases=...)`` documented this filtering for a
# Long time without implementing it; the mapping is now data, and the filter
# Honours it.
#
# Deliberately absent: ``diamond`` (modules/ortholog_resolver.py:79,84) belongs
# To the Phase 1.5 ortholog-resolution path, which the pipeline does not run --
# See phases.PHASES["1.5"] and the disclosure item in ``--check``. Listing it
# Here would tell users to install a tool that cannot change their results.
#
# Note iqtree3: it is NOT a Phase 3 tool. The only call site is
# ``utils/gtdb_tk_markers.py:287``, reached from marker selection (``scan``) via
# ``--gene-tree-builder iqtree``; the wrong table said "(Phase 3)".
TOOL_STEPS: Dict[str, List[str]] = {
    "hmmsearch": ["scan"],
    "checkm": ["scan"],
    "mafft": ["scan", "filter", "infer"],
    "trimal": ["scan", "infer"],
    "fasttree": ["scan", "filter", "infer"],
    "iqtree3": ["scan"],
    "astral": ["infer"],
}

# The four resumable steps exposed by the CLI (scan/filter/infer/report).
PIPELINE_STEPS: List[str] = ["scan", "filter", "infer", "report"]

# Alternative binary names for tools with platform-dependent naming.
_TOOL_ALIASES: Dict[str, List[str]] = {
    "fasttree": ["fasttree", "FastTree"],
}


def check_python_dependencies(strict: bool = False) -> List[str]:
    """Verify critical Python package imports at startup.

    Returns a list of warning/error messages. If strict=True, raises
    ImportError on critical failures.
    """
    messages: List[str] = []

    # Ete3 is required for HGT monophyly screening and tree operations
    try:
        # Single entry point, and its message carries the REAL reason.
        from markerfinder.utils.etree import require_ete3

        require_ete3()
    except Exception as e:  # Uninstalled, or installed but unimportable
        # "not installed" used to be asserted whenever the import failed, which
        # is wrong when ete3 IS installed but cannot be imported. The usual cause
        # was Python 3.13+, where the stdlib `cgi` module ete3 imports at
        # package-import time no longer exists; `markerfinder._cgi_compat`
        # normally papers over that, so a failure here means something else is
        # wrong. Report what actually happened (name the cause, never mislabel
        # it) and do not send the user downgrading an interpreter that is fine.
        msg = (
            f"ete3 is unusable in this interpreter ({e}). The phylogenetic HGT "
            "screening (MAD rooting + monophyly proportion) will not be "
            "available. Every supported interpreter (Python >=3.10) can import "
            "ete3 once `markerfinder` has been imported — so this is a broken "
            "installation, not an unsupported interpreter. Reinstall with "
            "`conda install -c etetoolkit ete3` or `pip install ete3 six` "
            "(the PyPI wheel does not declare six)."
        )
        if strict:
            raise ImportError(msg)
        messages.append(msg)
        logger.warning(msg)

    # Biopython is critical for sequence/tree parsing
    try:
        import Bio  # Noqa: F401
    except ImportError:
        msg = "Biopython is not installed. This is a critical dependency."
        if strict:
            raise ImportError(msg)
        messages.append(msg)
        logger.error(msg)

    return messages


def check_external_tools(
    required_phases: List[str] | None = None,
) -> Dict[str, bool]:
    """Check availability of external CLI tools on PATH.

    Args:
        required_phases: If provided, only check tools needed for these steps.
            E.g. ["scan"] only checks hmmsearch/checkm; ["filter", "infer"]
            checks alignment/tree tools too. Names come from ``PIPELINE_STEPS``.

    Returns:
        Dict mapping tool name to availability (True/False).

    Raises:
        ValueError: if ``required_phases`` names a step that does not exist.
            Silently ignoring a typo'd step once meant "check nothing" read as
            "everything is present".
    """
    if required_phases is not None:
        unknown = sorted(set(required_phases) - set(PIPELINE_STEPS))
        if unknown:
            raise ValueError(
                f"unknown step(s) {unknown}; valid steps are {PIPELINE_STEPS}"
            )

    results: Dict[str, bool] = {}

    for tool_name, purpose, is_critical in EXTERNAL_TOOLS:
        if required_phases is not None and not any(
            step in required_phases for step in TOOL_STEPS.get(tool_name, [])
        ):
            continue
        aliases = _TOOL_ALIASES.get(tool_name, [tool_name])
        found = any(shutil.which(alias) for alias in aliases)
        results[tool_name] = found

        if not found:
            level = "ERROR" if is_critical else "WARNING"
            logger.log(
                logging.ERROR if is_critical else logging.WARNING,
                f"[Dependency Check] {level}: '{tool_name}' not found on PATH. "
                f"Purpose: {purpose}. "
                f"Tried: {', '.join(aliases)}.",
            )

    return results


def precheck_all(strict: bool = False) -> None:
    """Run all dependency pre-checks.

    Called once at pipeline startup. Logs all issues; if strict=True,
    raises on critical failures.

     requires one distinguishable report entry per outage CLASS. The
    external-tool summary is therefore emitted whether or not anything is
    missing: a run on a machine that has every tool must still show, in the
    same place, that the tool class was checked and came back clean. Scoring
    the class only when it fails made "no line at all" ambiguous between
    "all tools present" and "the check never ran".
    """
    logger.debug("Running dependency pre-checks...")
    py_issues = check_python_dependencies(strict=strict)
    tool_status = check_external_tools()

    n_missing = sum(1 for v in tool_status.values() if not v)
    missing_names = sorted(t for t, ok in tool_status.items() if not ok)
    if n_missing:
        logger.info(
            f"[Dependency Check] external tools: {n_missing}/{len(tool_status)} "
            f"not found ({', '.join(missing_names)}). Some pipeline phases may be "
            "unavailable. This entry names the TOOL class only; a missing "
            "Python package is reported separately above."
        )
    else:
        logger.info(
            f"[Dependency Check] external tools: {len(tool_status)}/{len(tool_status)} "
            f"found on PATH ({', '.join(sorted(tool_status))})."
        )
    if not n_missing and not py_issues:
        logger.debug("[Dependency Check] All dependencies satisfied.")
