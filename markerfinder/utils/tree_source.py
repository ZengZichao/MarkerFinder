"""Species-tree source priority resolution.

This module is the single source of truth for choosing the "primary" species
tree source among the available phylogenetic results. The logic used to be
duplicated in three places:

* ``markerfinder/pipeline.py`` (the availability signal that feeds
  ``PipelineResult.species_tree_source``), and
* ``markerfinder/modules/report_generator.py`` (the Interactive and PlainText
  report generators both resolved the same priority).

All three sites now call:func:`prioritize_tree_source` so the priority ranking
cannot silently drift between the pipeline result and the reports.
"""

from typing import Optional

# Priority ranking of coalescent species-tree sources (higher = more reliable).
# Astral >= consensus > concat (supermatrix) > first_gene_tree >= none.
# 'first_gene_tree' is a *degenerate* fallback (a single gene tree used when a
# Proper coalescent merge such as ASTRAL-III fails); its reliability is lower
# Than the IQ-TREE3 concatenation tree, so its priority value sits below
# 'concat'.
TREE_SOURCE_PRIORITY = {
    "astral": 4,
    "consensus": 3,
    "concat": 2,            # Concatenation / supermatrix
    "first_gene_tree": 1,   # Degenerate single-gene-tree fallback
    "none": 0,
}

# Default display label returned when the concatenation tree is the chosen
# Source. The reports use this human-readable label; callers that need a bare
# Machine-readable key (e.g. ``PipelineResult.species_tree_source``) pass
# ``concat_label="concat"``.
DEFAULT_CONCAT_LABEL = "concat (supermatrix)"


def prioritize_tree_source(
    coalescent_source: Optional[str],
    has_concat_tree: bool,
    concat_label: str = DEFAULT_CONCAT_LABEL,
) -> str:
    """Resolve the primary species-tree source.

    Args:
        coalescent_source: the coalescent result's ``species_tree_source``
            (e.g. ``"astral"``, ``"consensus"``, ``"first_gene_tree"``,
            ``"none"`` or ``None``). Callers that have already determined the
            coalescent tree is absent should pass ``None``.
        has_concat_tree: whether a concatenation / supermatrix tree is
            available.
        concat_label: the string returned when the concatenation tree is the
            chosen source. Defaults to ``"concat (supermatrix)"``; pass
            ``"concat"`` when a bare key is required.

    Returns:
        The resolved primary source string.
    """
    source = coalescent_source or "none"
    rank = TREE_SOURCE_PRIORITY.get(source, TREE_SOURCE_PRIORITY["none"])
    if rank >= TREE_SOURCE_PRIORITY["concat"]:
        # Astral / consensus / concat: a proper merged tree wins outright.
        return source
    if has_concat_tree:
        # No proper merged tree, but a concatenation tree exists: prefer it
        # Over a degenerate single-gene-tree fallback.
        return concat_label
    if rank >= TREE_SOURCE_PRIORITY["first_gene_tree"]:
        # No concatenation tree: fall back to the degenerate single gene tree
        # If one is available.
        return source
    return "none"


def coalescent_is_degenerate(coalescent_source: Optional[str]) -> bool:
    """Return ``True`` when the coalescent source is the degenerate single-gene-tree fallback."""
    return coalescent_source in ("first_gene_tree",)


def has_usable_tree(tree: Optional[object]) -> bool:
    """Return ``True`` only for a *usable* inferred tree (present and >= 4 tips).

    Dataclass instances are always truthy, so a bare ``bool(tree)``
    previously counted a placeholder ``;`` tree (0 tips) as a valid
    concatenation tree and let a degraded run report success. All
    availability checks (pipeline usability signal and both report
    generators) share this helper so the threshold cannot drift.
    """
    if tree is None:
        return False
    n_tips = getattr(tree, "n_tips", None)
    if n_tips is None:
        return True
    try:
        return int(n_tips) >= 4
    except (TypeError, ValueError):
        return True
