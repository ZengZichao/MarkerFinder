"""Threshold scan.

Re-grades the already-computed per-marker risks under several threshold
bands — **no external tool is re-run**: band comparison is
O(bands x markers) over the stored ``overall_risk`` values.

Output: ``Phase5_reports/{prefix}.threshold_scan.tsv`` with one row per band
``(level1_max, level2_max) -> (n_L1, n_L2, n_L3, n_UNKNOWN, n_selected,
n_at_boundary, mean_risk)`` plus a cross-band Jaccard matrix of the selected
marker sets. ``n_at_boundary`` exists so readers do not mistake boundary
effects for sensitivity (strict-less-than grading, hgt_filter.py).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from markerfinder.models.marker import MarkerLevel
from markerfinder.models.pipeline_types import HGTEvaluation

logger = logging.getLogger(__name__)

# E.2: default scan grid, far band included.
DEFAULT_BANDS: Tuple[Tuple[float, float], ...] = (
    (0.15, 0.60),
    (0.25, 0.60),
    (0.40, 0.60),
    (0.25, 0.95),  # Far band
)


@dataclass
class BandResult:
    level1_max: float
    level2_max: float
    n_L1: int = 0
    n_L2: int = 0
    n_L3: int = 0
    n_unknown: int = 0
    n_selected: int = 0
    n_at_boundary: int = 0
    mean_risk: float = 0.0
    selected: Set[str] = field(default_factory=set)


def _grade(risk: float, l1: float, l2: float) -> MarkerLevel:
    # Strict-less-than, mirroring HGTDecisionEngine.evaluate_marker
    if risk < l1:
        return MarkerLevel.LEVEL_1
    if risk < l2:
        return MarkerLevel.LEVEL_2
    return MarkerLevel.LEVEL_3


def run_threshold_scan(
    evaluations: List[HGTEvaluation],
    output_dir: str,
    prefix: str,
    *,
    bands: Optional[List[Tuple[float, float]]] = None,
    stability_min: float = 0.6,
) -> List[BandResult]:
    """Re-grade stored risks under each band and write
    ``Phase5_reports/{prefix}.threshold_scan.tsv``."""
    graded = [e for e in evaluations if e.level is not MarkerLevel.UNKNOWN]
    if bands is None:
        bands = list(DEFAULT_BANDS)

    results: List[BandResult] = []
    for l1, l2 in bands:
        band = BandResult(level1_max=l1, level2_max=l2)
        risk_values = []
        for ev in graded:
            risk = ev.overall_risk
            level = _grade(risk, l1, l2)
            if level is MarkerLevel.LEVEL_1:
                band.n_L1 += 1
                band.selected.add(ev.marker_id)
            elif level is MarkerLevel.LEVEL_2:
                band.n_L2 += 1
                band.selected.add(ev.marker_id)
            else:
                band.n_L3 += 1
            if risk in (l1, l2):
                band.n_at_boundary += 1
            risk_values.append(risk)
        band.n_unknown = len(evaluations) - len(graded)
        band.n_selected = band.n_L1 + band.n_L2
        band.mean_risk = (
            sum(risk_values) / len(risk_values) if risk_values else 0.0
        )
        results.append(band)

    _write_tsv(results, output_dir, prefix)
    _warn_on_instability(results, floor=stability_min)
    return results


def _jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    return len(a & b) / len(union) if union else 1.0


def _warn_on_instability(results: List[BandResult], floor: float = 0.6) -> None:
    """Flag bands whose selected set is unstable across the scan."""
    if len(results) < 2:
        return
    worst = 1.0
    for i in range(len(results)):
        for j in range(i + 1, len(results)):
            jacc = _jaccard(results[i].selected, results[j].selected)
            worst = min(worst, jacc)
    if worst < floor:
        logger.warning(
            f"  Threshold scan: minimum cross-band Jaccard {worst:.2f} < "
            f"{floor:.2f} — the selected-marker set is unstable across bands; "
            f"treat band-specific conclusions as inconclusive."
        )


def product_path(output_dir: str, prefix: str) -> Path:
    """Where the scan's artifact lands -- one definition, used by writer and caller."""
    return Path(output_dir) / "Phase5_reports" / f"{prefix}.threshold_scan.tsv"


@dataclass(frozen=True)
class ThresholdScanOutcome:
    """Bands re-graded under each threshold, plus the file that records them."""

    results: List[BandResult]
    path: Path


def run_stage(
    evaluations: List[HGTEvaluation],
    output_dir: str,
    prefix: str,
    *,
    bands: Optional[List[Tuple[float, float]]] = None,
    stability_min: float = 0.6,
) -> ThresholdScanOutcome:
    """``run_threshold_scan`` plus the artifact it wrote, so callers can name it."""
    results = run_threshold_scan(
        evaluations, output_dir, prefix, bands=bands, stability_min=stability_min
    )
    return ThresholdScanOutcome(results=results, path=product_path(output_dir, prefix))


def _write_tsv(results: List[BandResult], output_dir: str, prefix: str) -> Path:
    out = Path(output_dir) / "Phase5_reports"
    out.mkdir(parents=True, exist_ok=True)
    path = product_path(output_dir, prefix)

    lines = [
        "band_level1_max\tband_level2_max\tn_L1\tn_L2\tn_L3\tn_UNKNOWN\t"
        "n_selected\tn_at_boundary\tmean_risk\n"
    ]
    for band in results:
        lines.append(
            f"{band.level1_max:.2f}\t{band.level2_max:.2f}\t{band.n_L1}\t"
            f"{band.n_L2}\t{band.n_L3}\t{band.n_unknown}\t{band.n_selected}\t"
            f"{band.n_at_boundary}\t{band.mean_risk:.4f}\n"
        )

    # Cross-band Jaccard matrix of the selected sets.
    lines.append("\n")
    labels = [f"({b.level1_max:.2f},{b.level2_max:.2f})" for b in results]
    lines.append("jaccard\t" + "\t".join(labels) + "\n")
    for i, bi in enumerate(results):
        row = [labels[i]]
        for j, bj in enumerate(results):
            row.append(f"{_jaccard(bi.selected, bj.selected):.3f}")
        lines.append("\t".join(row) + "\n")

    target = path.resolve()
    target.write_text("".join(lines), encoding="utf-8", newline="\n")
    logger.info(f"  Threshold scan written: {target}")
    return target
