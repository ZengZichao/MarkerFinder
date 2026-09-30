"""Four-state evidence model.

``None`` alone cannot distinguish "could not be measured" from "does not
apply on this path"; both need different user action. Every measurement that
feeds an HGT decision carries a:class:`MeasureState` and, when it was not
measured, a machine-readable reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class MeasureState(Enum):
    MEASURED = "measured"          # Actually measured
    NOT_MEASURABLE = "not_measurable"  # Dependency/tool/data missing
    NOT_APPLICABLE = "not_applicable"  # This decision path has no such signal
    REJECTED = "rejected"          # Measured, but reference legality failed


class UnknownReason(Enum):
    """Distinguishable sources of ``MarkerLevel.UNKNOWN``."""

    NO_REFERENCE = "no_reference"              # No species tree and no taxonomy table
    TREE_BUILD_FAILED = "tree_build_failed"    # Gene tree could not be built
    NO_TESTABLE_RANK = "no_testable_rank"      # No rank has >=2 representatives
    RF_UNMEASURABLE = "rf_unmeasurable"        # RF could not be computed
    QUARTET_UNMEASURABLE = "quartet_unmeasurable"  # No quartet could be evaluated
    REFERENCE_ILLEGAL = "reference_illegal"    # Reference tree failed legality checks
    SINGLE_SIGNAL_ONLY = "single_signal_only"  # Synthesis needs >=2 real signals


@dataclass
class Measurement:
    """One measured (or unmeasured) quantity with its evidence state."""

    value: Optional[float]
    state: MeasureState
    method: str = ""     # E.g. "species_tree_rf" | "taxonomy_monophyly" | "splits-python"
    reason: str = ""     # Required when state is NOT_MEASURABLE / REJECTED
    n_obs: int = 0       # Observation base (quartet count, n_total,...)

    @property
    def measured(self) -> bool:
        return self.state is MeasureState.MEASURED and self.value is not None
