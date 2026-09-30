"""Canonical phase numbering for MarkerFinder.

This module is the SINGLE source of truth for pipeline phase numbers and the
historical output-directory aliases. Directory names are deliberately NOT
renamed (user scripts and an integration test assert on
``Phase4_trees/`` / ``Phase5_reports/``); only log labels and comments are
unified, with each alias registered here so tooling can map old→new.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


@dataclass(frozen=True)
class Phase:
    number: str          # Canonical label, e.g. "0.1"
    title: str
    log_tag: str         # Exact tag used in log lines, e.g. "[Phase 0.1]"
    directory: str       # Output directory associated with this phase


PHASES: List[Phase] = [
    Phase("0",   "Quality preprocessing (CheckM / defaults)",        "[Phase 0]",   ""),
    Phase("0.1", "Reference species tree intake",                    "[Phase 0.1]", ""),
    Phase("0.2", "Taxonomy table intake",                            "[Phase 0.2]", ""),
    Phase("1",   "Marker selection (GTDB-TK / HMM)",                 "[Phase 1]",   ""),
    Phase("1.5", "Marker sequence extraction (ortholog resolution NOT wired in — F-2)", "[Phase 1.5]", ""),
    Phase("2",   "HGT-aware marker filtering (gene trees + grading)", "[Phase 2]",   "Phase4_trees/gene_trees"),
    Phase("3",   "Phylogenetic inference (concatenation + coalescent)", "[Phase 3]", "Phase4_trees"),
    Phase("4",   "Reports, metadata and evidence artifacts",          "[Phase 4]",   "Phase5_reports"),
]

# Historical aliases: legacy names that STILL appear in output paths (
# Directories unchanged) mapped to their canonical phase. Any comment or log
# Line referring to "Phase 4" as the coalescent step or "Phase 5" as reports
# Is a historical alias of phase 3 / 4 respectively.
HISTORICAL_ALIASES: Dict[str, str] = {
    "Phase4_trees": "3",          # Directory name kept; belongs to phase 3 (and 2 caches)
    # Added by the ratchet (tests/unit/test_help_phase_labels_are_consistent.py),
    # Which scans every PhaseN_dir token in the package and fails on any that is
    # Not registered here. It immediately found Phase4_alignments: created by
    # Report_generator.py:343, listed as reserved in validation.py:41 and
    # Documented in the `report` subcommand help -- a real product directory that
    # The alias table had simply never heard of, i.e. the exact "指向错误"
    # Is about. Per-marker alignments/partitions are produced by the inference
    # Step, so the canonical phase is 3.
    "Phase4_alignments": "3",
    "Phase4_intermediate": "3",
    "Phase5_reports": "4",
    "Phase5_metadata": "4",
    "Phase5_evidence": "4",
    "Phase 4 (CoalescentInference)": "3",  # Pre- comment wording
    "Step 4: report (Phase 5)": "4",       # Pre- comment wording
}


def get_phase(number: str) -> Phase:
    for phase in PHASES:
        if phase.number == number:
            return phase
    raise KeyError(f"unknown phase {number!r}")


def canonical_log_tag(number: str) -> str:
    return get_phase(number).log_tag
