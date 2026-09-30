"""Reference-tree legality validation.

Direct defence against "Error A" of the retracted paper — Steenwyk, J. L. &
King, N. (2025). Integrative phylogenomics positions sponges at the root of
the animal tree. Science 390(6774), 751-756, doi:10.1126/science.adw9456;
retracted 2026-02-05 (retraction notice: Thorp, H. H., 2026. Science
391(6785), 564, doi:10.1126/science.aef5589). Error A is the
merged-reference defect: a reference object that is
unresolved, too small, does not cover the gene tree's tips, or has the
focal clades already collapsed must NEVER be used as a scoring reference.
"Legality" here means availability (resolved, covering, not collapsed) —
not biological truth.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, FrozenSet, List, Optional, Set

from markerfinder.exceptions import PhyloFormatError
from markerfinder.utils import etree

# Reference legality codes.
POLYTOMY = "POLYTOMY"
TIPSET_MISMATCH = "TIPSET_MISMATCH"
TOO_FEW_TIPS = "TOO_FEW_TIPS"
MALFORMED = "MALFORMED"
CLADE_COLLAPSED = "CLADE_COLLAPSED"

# Legality checks reuse the pure-Python measurement cap.
MAX_REFERENCE_TIPS = etree.MAX_PURE_PYTHON_TIPS


@dataclass
class ReferenceVerdict:
    """Outcome of a reference-tree legality check."""

    ok: bool
    code: Optional[str] = None
    detail: str = ""


def _is_resolved(newick: str) -> bool:
    """True iff every internal node has unrooted degree <= 3 (no polytomy).

    For a Newick writing, an internal node's unrooted degree is its child
    count (root) or child count + 1 (non-root). Any degree > 3 is a polytomy.
    """
    root = etree._parse(newick)  # Noqa: SLF001 — same-package internal reuse

    def walk(node: Any) -> int:
        """Return subtree size; raise on polytomy."""
        if node.is_leaf:
            return 1
        degree = len(node.children) + (0 if node is root else 1)
        if degree > 3:
            raise PhyloFormatError(f"polytomy: internal node with degree {degree}")
        total = 0
        for child in node.children:
            total += walk(child)
        return total

    walk(root)
    return True


def _clade_collapsed(
    ref_splits: Set[FrozenSet[str]],
    all_tips: Set[str],
    focal_groups: List[Set[str]],
) -> Optional[str]:
    """Detect Error-A-style collapse among focal groups.

    If some split of the reference has one side exactly equal to the union
    of two focal groups (or its complement), those two groups — which the
    competing hypotheses require to be independent — are already merged
    into a single reference edge, and the reference cannot arbitrate
    between them.
    """
    for a, b in combinations(focal_groups, 2):
        union = frozenset(set(a) | set(b))
        complement = frozenset(all_tips - set(union))
        for split in ref_splits:
            if split == union or split == complement:
                return (
                    f"focal groups {sorted(a)} and {sorted(b)} are merged on a "
                    f"single reference edge (split side {sorted(split)})"
                )
    return None


def validate_reference_tree(
    ref_newick: str,
    gene_newick: str,
    *,
    target_clades: Optional[List[Set[str]]] = None,
) -> ReferenceVerdict:
    """Validate a reference tree before it may participate in scoring.

    Checks:
      1. the reference parses and passes structural sanity (MALFORMED);
      2. it is fully resolved — no polytomy (POLYTOMY);
      3. tipset(ref) ⊇ tipset(gene) on this marker's tip set (TIPSET_MISMATCH);
      4. the reference has >= 4 tips (TOO_FEW_TIPS);
      5. optional: the focal clades of the competing hypotheses are pairwise
         not collapsed onto a single reference edge (CLADE_COLLAPSED).

    A failed check yields ``ok=False`` with a distinguishable ``code``; the
    caller must route the marker to UNKNOWN/REJECTED rather than score.
    """
    # 1. Malformed / unparseable.
    try:
        ref_tips = etree.tip_set(ref_newick)
        etree.tip_set(gene_newick)
        splits, all_ref_tips = etree.splits_of(ref_newick)
    except etree.TreeMeasureError as e:
        return ReferenceVerdict(False, MALFORMED, f"reference newick unparseable: {e}")
    except Exception as e:  # Defensive: any parser error is MALFORMED
        return ReferenceVerdict(False, MALFORMED, f"reference newick unparseable: {e}")

    # 3 (ordering: cheap count first for a precise message) — tip coverage.
    gene_tips = etree.tip_set(gene_newick)
    missing = gene_tips - ref_tips
    if missing:
        return ReferenceVerdict(
            False,
            TIPSET_MISMATCH,
            f"reference lacks {len(missing)} gene-tree tip(s): {sorted(missing)[:5]}",
        )

    # 4. Minimum tip count.
    if len(ref_tips) < 4:
        return ReferenceVerdict(
            False, TOO_FEW_TIPS, f"reference has {len(ref_tips)} tips (<4)"
        )

    # 2. Full resolution. Checked after cheap counting so degenerate inputs
    # Report the more specific code first.
    try:
        _is_resolved(ref_newick)
    except PhyloFormatError as e:
        return ReferenceVerdict(False, POLYTOMY, str(e))

    # 5. Optional focal-clade collapse check (Error A).
    if target_clades:
        if len(ref_tips) > MAX_REFERENCE_TIPS:
            # Beyond the pure-Python cap the split set is not measurable
            #; collapse checking is refused, not fabricated.
            return ReferenceVerdict(
                False,
                "TIP_CAP_EXCEEDED",
                f"{len(ref_tips)} reference tips exceed cap {MAX_REFERENCE_TIPS}",
            )
        collapsed = _clade_collapsed(splits, all_ref_tips, target_clades)
        if collapsed:
            return ReferenceVerdict(False, CLADE_COLLAPSED, collapsed)

    return ReferenceVerdict(True, None, "reference tree is legal")
