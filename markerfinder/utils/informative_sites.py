"""Per-marker phylogenetic informativeness.

PIS (parsimony-informative sites, Farris definition, protein alphabet): a
column is informative when at least two different residues each occur at
least twice. Gaps and ambiguous characters ('-', 'X', '?', 'B', 'Z') do not
count as residues. Pure stdlib; O(sites x taxa).
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple

_GAP_CHARS = set("-.~")
_AMBIGUOUS = set("XBZJUO?*")  # Treated like gaps for informativeness purposes
_NON_RESIDUE = _GAP_CHARS | _AMBIGUOUS


def _as_columns(alignment: Dict[str, str]) -> Tuple[List[str], int]:
    """Transpose {id: seq} into column strings; validates equal lengths."""
    if not alignment:
        return [], 0
    seqs = [str(s).upper() for s in alignment.values()]
    length = len(seqs[0])
    for seq_id, seq in alignment.items():
        if len(seq) != length:
            raise ValueError(
                f"unaligned sequence {seq_id!r}: length {len(seq)} != {length}"
            )
    columns = ["".join(seq[i] for seq in seqs) for i in range(length)]
    return columns, length


def parsimony_informative_sites(alignment: Dict[str, str]) -> int:
    """Count columns where >=2 distinct residues each appear >=2 times."""
    columns, _ = _as_columns(alignment)
    informative = 0
    for column in columns:
        counts: Dict[str, int] = {}
        for ch in column:
            if ch in _NON_RESIDUE or ch.isspace():
                continue
            counts[ch] = counts.get(ch, 0) + 1
        residues_with_min_two = sum(1 for c in counts.values() if c >= 2)
        if residues_with_min_two >= 2:
            informative += 1
    return informative


def effective_columns(alignment: Dict[str, str]) -> int:
    """Columns kept after removing entirely-gapped columns."""
    columns, _ = _as_columns(alignment)
    kept = 0
    for column in columns:
        if any(ch not in _GAP_CHARS and not ch.isspace() for ch in column):
            kept += 1
    return kept


def gap_fraction(alignment: Dict[str, str]) -> float:
    """Fraction of gap characters over the whole alignment."""
    columns, length = _as_columns(alignment)
    if length == 0 or not alignment:
        return 0.0
    total = length * len(alignment)
    gaps = sum(
        1
        for column in columns
        for ch in column
        if ch in _GAP_CHARS or ch.isspace()
    )
    return gaps / total if total else 0.0


def apply_pis_floor(evaluations: Iterable[Any], pis_by_marker: Dict[str, int],
                    floor: int) -> List[str]:
    """A marker with too few parsimony-informative sites is inconclusive.

    "Support from too little information is not support" — a sequence-level
    analogue of the ``|dlnL| > 2`` rule of the retracted sponge paper
    (Steenwyk & King 2025, Science, doi:10.1126/science.adw9456). A marker whose PIS is below
    ``floor`` gets its decision card and notes stamped ``inconclusive`` even
    when its risk score is low, so the report cannot present it as clean.

    ``floor <= 0`` (the shipped default) is a deliberate no-op: the threshold
    is provisional until supplies real data, and requires
    the default behaviour to equal the baseline. Markers whose PIS was never
    measured are left alone rather than guessed at.

    Returns the demoted marker ids. Grading itself (``ev.level``) is NOT
    touched — exclusion from the usable set is a publication decision for the
    author, not something this layer invents.
    """
    from markerfinder.models.marker import MarkerLevel

    if not evaluations or not floor or int(floor) <= 0:
        return []
    floor = int(floor)
    demoted: List[str] = []
    for ev in evaluations:
        if getattr(ev, "level", None) is MarkerLevel.UNKNOWN:
            continue
        pis = pis_by_marker.get(getattr(ev, "marker_id", ""))
        if pis is None:
            continue
        if pis >= floor:
            continue
        card = getattr(ev, "decision_card", None)
        if isinstance(card, dict):
            card["pis_grade"] = "inconclusive"
            card["pis_floor"] = floor
            card["pis_measured"] = pis
        message = (
            f"PIS {pis} < floor {floor}: inconclusive (provisional "
            f"threshold, not yet calibrated on real data)"
        )
        existing = getattr(ev, "notes", "") or ""
        ev.notes = f"{existing} {message}".strip()
        demoted.append(ev.marker_id)
    return demoted
