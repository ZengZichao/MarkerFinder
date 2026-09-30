"""Composition / GC bias diagnostics.

Default OFF (``enable_composition_screen=False``). These metrics are
PARALLEL EVIDENCE COLUMNS: they must never be merged into ``overall_risk`` or
``overall_score`` (locked by test). Protein inputs have no GC — the
GC metric is explicitly NOT_APPLICABLE for.faa (never 0).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

_AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def aa_composition(seq: str) -> Dict[str, float]:
    """Normalised amino-acid frequencies (gaps/ambiguous excluded)."""
    counts: Dict[str, int] = {}
    total = 0
    for ch in seq.upper():
        if ch in _AA20:
            counts[ch] = counts.get(ch, 0) + 1
            total += 1
    if total == 0:
        return {}
    return {aa: c / total for aa, c in counts.items()}


def rcv(seqs: Sequence[str], group_labels: Optional[List[str]] = None) -> float:
    """Relative Compositional Variability (per Otu et al.-style RCV).

     implement RCV or an equivalent chi-square/Euclidean
    composition-outlier metric. With no group labels, measures global
    compositional variability across sequences; with labels, compares each
    group against the pooled composition and returns the worst deviation.
    """
    if not seqs:
        return 0.0
    compositions = [aa_composition(s) for s in seqs]
    compositions = [c for c in compositions if c]
    if not compositions:
        return 0.0

    aas = sorted({aa for comp in compositions for aa in comp})
    n = len(compositions)

    def _deviation(comp: Dict[str, float], pooled: Dict[str, float]) -> float:
        return sum(abs(comp.get(aa, 0.0) - pooled.get(aa, 0.0)) for aa in aas) / 2.0

    pooled: Dict[str, float] = {}
    for comp in compositions:
        for aa, freq in comp.items():
            pooled[aa] = pooled.get(aa, 0.0) + freq / n

    if not group_labels:
        return max(_deviation(c, pooled) for c in compositions)

    # Group-wise: worst between-group deviation from the pooled composition
    worst = 0.0
    seen_groups = sorted(set(group_labels))
    for group in seen_groups:
        members = [compositions[i] for i, g in enumerate(group_labels) if g == group]
        if not members:
            continue
        group_comp: Dict[str, float] = {}
        for comp in members:
            for aa, freq in comp.items():
                group_comp[aa] = group_comp.get(aa, 0.0) + freq / len(members)
        worst = max(worst, _deviation(group_comp, pooled))
    return worst


def gc_content(seq: str, protein: bool = False) -> Optional[float]:
    """GC fraction of a NUCLEOTIDE sequence.

    Protein input (``protein=True``) is NOT_APPLICABLE — returns None, never
    0 (0 would read as "no GC bias"). Protein/nucleotide alphabets share
    A/C/G/T letters, so the caller must declare the molecule type.
    """
    if protein:
        return None
    seq_u = seq.upper()
    total = sum(1 for ch in seq_u if ch in "ACGTUN")
    if total == 0:
        return None
    gc = sum(1 for ch in seq_u if ch in "GC")
    return gc / total


def gc_or_codon_bias(seqs: Sequence[str], protein: bool = True) -> Optional[float]:
    """Spread of GC content across sequences.

    Protein input (default,.faa) => NOT_APPLICABLE (None): GC diagnostics do
    not apply to amino-acid data and must never be reported as 0.
    """
    if protein:
        return None
    raw_values = [gc_content(s) for s in seqs]
    values = [v for v in raw_values if v is not None]
    if len(values) < 2:
        return None
    mean = sum(values) / len(values)
    variance = sum((v - mean) ** 2 for v in values) / len(values)
    return math.sqrt(variance)


def chi2_vs_group_mean(seq: str, group_seqs: Sequence[str]) -> Optional[float]:
    """Chi-square-like statistic of one sequence's composition vs group mean.

    Returns None when the composition cannot be computed (empty inputs).
    """
    comp = aa_composition(seq)
    if not comp or not group_seqs:
        return None
    members = [aa_composition(s) for s in group_seqs]
    members = [m for m in members if m]
    if not members:
        return None
    aas = sorted({aa for m in members for aa in m} | set(comp))
    n = len(members)
    pooled = {aa: sum(m.get(aa, 0.0) for m in members) / n for aa in aas}
    chi2 = 0.0
    for aa in aas:
        expected = max(pooled[aa], 1e-6)
        observed = comp.get(aa, 0.0)
        chi2 += (observed - expected) ** 2 / expected
    return chi2


def marker_composition_metrics(
    aln_map: Dict[str, str], protein: bool = True
) -> Dict[str, Any]:
    """Per-marker composition diagnostics from a ``{taxon: aligned_seq}`` map.

    Produces the parallel evidence columns promised: ``rcv``
    and ``gc_bias``. ``gc_bias`` is ``None`` for protein input — GC does not
    apply to amino-acid data and must never be reported as 0 (that is the
    "unmeasurable encoded as a number" defect this whole refactor exists to
    remove). One pass, pure stdlib.
    """
    seqs = [s for s in (aln_map or {}).values() if s]
    if len(seqs) < 2:
        return {"rcv": None, "gc_bias": None, "n_sequences": len(seqs)}
    return {
        "rcv": rcv(seqs),
        "gc_bias": gc_or_codon_bias(seqs, protein=protein),
        "n_sequences": len(seqs),
    }


def write_composition_tsv(rows: Sequence[Dict[str, Any]], output_dir: str,
                          prefix: str,
                          warn_threshold: Optional[float] = None) -> Path:
    """Write ``{prefix}.composition.tsv`` — parallel columns only.

    Sorted worst-outlier first so a planted composition outlier lands at the
    head of the list. Nothing here is combined into ``overall_risk``
    or ``overall_score``: composition is shown *beside* the HGT verdict so a
    reader can see both, which is the whole point of.
    """
    out = Path(output_dir) / "Phase5_reports"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{prefix}.composition.tsv"

    def _key(row: Dict[str, Any]) -> Tuple[float, str]:
        value = row.get("rcv")
        return (-(value if value is not None else -1.0), str(row.get("marker_id", "")))

    lines = ["marker_id\trcv\tgc_bias\tn_sequences\toutlier_flag\n"]
    for row in sorted(rows, key=_key):
        rcv_value = row.get("rcv")
        gc_value = row.get("gc_bias")
        flag = "NA"
        if (warn_threshold is not None and rcv_value is not None
                and rcv_value >= warn_threshold):
            flag = "warn"
        lines.append(
            "{mid}\t{rcv}\t{gc}\t{n}\t{flag}\n".format(
                mid=row.get("marker_id", ""),
                rcv="NA" if rcv_value is None else f"{rcv_value:.6f}",
                gc="NA" if gc_value is None else f"{gc_value:.6f}",
                n=row.get("n_sequences", "NA"),
                flag=flag,
            )
        )
    path.resolve().write_text("".join(lines), encoding="utf-8", newline="\n")
    return path


@dataclass(frozen=True)
class CompositionOutcome:
    """What the composition stage produced.

    a stage returns its own products, so the orchestrator does not have to
    know the file name, the row count or the empty-input rule to report them.
    """

    path: Optional[Path]
    rows_written: int


def run_stage(rows: Sequence[Dict[str, Any]], output_dir: str, prefix: str, *,
              warn_threshold: Optional[float] = None) -> CompositionOutcome:
    """The whole composition product, in one call.

    Empty input is neither an error nor an artifact: the stage says so by
    returning ``path=None``, which keeps the caller's log line truthful instead
    of letting it claim "diagnostics written" over an empty dict.
    """
    materialised = list(rows)
    if not materialised:
        return CompositionOutcome(path=None, rows_written=0)
    return CompositionOutcome(
        path=write_composition_tsv(
            materialised, output_dir, prefix, warn_threshold=warn_threshold
        ),
        rows_written=len(materialised),
    )
