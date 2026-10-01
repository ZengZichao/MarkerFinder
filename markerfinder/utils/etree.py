"""Single ete3 import entry point + pure-Python tree measurement fallback.

 of:

*:func:`require_ete3` is the ONLY sanctioned place to import ete3. Scattered
  ``try: from ete3...`` blocks that swallow failures into placeholder numbers
  are the root cause of the zero-discrimination defect.
* When ete3 is unavailable (Python ≥ 3.13 removed ``cgi`` which ete3 still
  imports), the pure-Python split-set path below provides REAL measurements
  for RF distance, quartet topology and clade monophyly on trees of up to
:data:`MAX_PURE_PYTHON_TIPS` shared tips. These are measurements, not
  placeholders — callers report ``method="splits-python"`` when they are used.

The Newick parser here is deliberately minimal: it understands unquoted and
quoted tip labels (numeric labels such as ``"123"`` / ``"12.5"`` are kept as
labels, never parsed as support values), branch lengths, internal node labels
and NHX/ASTRAL bracket annotations (stripped defensively, same rules as
``utils.tree_utils.strip_nhx_annotations``).
"""

from __future__ import annotations

import re
from typing import Any, FrozenSet, List, Optional, Set, Tuple, Union

from markerfinder.exceptions import PhyloToolUnavailable

# Pure-Python measurement ceiling. Default max_markers=60, so tip sets
# Never exceed this in the default pipeline; going past it is a deliberate
# "not measurable" (reason=tip_cap_exceeded), never a silently degraded value.
MAX_PURE_PYTHON_TIPS = 64


class TreeMeasureError(Exception):
    """Base error for pure-Python tree measurement failures."""


class TipCapExceeded(TreeMeasureError):
    """Shared tip count exceeds MAX_PURE_PYTHON_TIPS (not measurable)."""


class NewickParseError(TreeMeasureError):
    """The Newick string could not be parsed by the minimal parser."""


def require_ete3() -> Any:
    """Return the ete3 module, or raise PhyloToolUnavailable.

    Every ete3 use in MarkerFinder must go through this function so that
    "dependency missing" has exactly one observable failure type.

    On Python 3.13+ the stdlib ``cgi`` module no longer exists and ete3's
    import-time ``from .webplugin.webapp import *`` therefore raises
    ``ModuleNotFoundError``. :mod:`markerfinder._cgi_compat` installs a minimal
    stand-in for it, restoring ete3 (and with it the real MAD / monophyly
    measurements) instead of letting every affected marker degrade to
    "unmeasured". On 3.12 and earlier the shim is a no-op.
    """
    try:
        from markerfinder import _cgi_compat  # noqa: F401, PLC0415

        import ete3  # Noqa: PLC0415 — deferred on purpose: import must be lazy

        return ete3
    except Exception as e:  # ImportError and broken ete3 installs alike
        raise PhyloToolUnavailable(
            f"ete3 is not importable in this interpreter: {e}"
        ) from e


# ---------------------------------------------------------------------------
# Minimal Newick parsing
# ---------------------------------------------------------------------------

class _Node:
    __slots__ = ("name", "label", "length", "children")

    def __init__(self) -> None:
        self.name: Optional[str] = None      # Leaf label (may be "" for empty)
        self.label: Optional[str] = None     # Internal node label (e.g. support)
        self.length: Optional[str] = None    # Branch length as written
        self.children: List["_Node"] = []

    @property
    def is_leaf(self) -> bool:
        return not self.children


_NAME_STOP = set(",():;")

def _strip_annotations(newick: str) -> str:
    if not newick:
        return newick
    stripped = re.sub(r"'[^']*'", "", newick)
    stripped = re.sub(r"\[.*?\]", "", stripped)
    return stripped.strip()


def _parse(newick: str) -> _Node:
    s = _strip_annotations(newick)
    if not s:
        raise NewickParseError("empty newick string")
    if not s.endswith(";"):
        s = s + ";"
    pos = 0
    n = len(s)

    def skip_ws() -> None:
        nonlocal pos
        while pos < n and s[pos].isspace():
            pos += 1

    def parse_name() -> str:
        nonlocal pos
        skip_ws()
        if pos < n and s[pos] in ("'", '"'):
            quote = s[pos]
            end = s.find(quote, pos + 1)
            if end < 0:
                raise NewickParseError(f"unterminated quoted label at {pos}")
            name = s[pos + 1:end]
            pos = end + 1
            return name
        start = pos
        while pos < n and s[pos] not in _NAME_STOP and not s[pos].isspace():
            pos += 1
        return s[start:pos]

    def parse_length() -> Optional[str]:
        nonlocal pos
        if pos >= n or s[pos] != ":":
            return None
        pos += 1
        start = pos
        while pos < n and s[pos] not in _NAME_STOP and s[pos] != "(":
            pos += 1
        return s[start:pos].strip()

    def parse_subtree() -> _Node:
        nonlocal pos
        skip_ws()
        if pos >= n:
            raise NewickParseError("unexpected end of newick")
        if s[pos] == "(":
            node = _Node()
            pos += 1
            while True:
                node.children.append(parse_subtree())
                skip_ws()
                if pos < n and s[pos] == ",":
                    pos += 1
                    continue
                break
            skip_ws()
            if pos >= n or s[pos] != ")":
                raise NewickParseError(f"unbalanced parentheses at {pos}")
            pos += 1
            # Optional internal label (e.g. bootstrap) then optional length.
            start = pos
            while pos < n and s[pos] not in _NAME_STOP and not s[pos].isspace():
                pos += 1
            tag = s[start:pos]
            node.length = parse_length()
            if tag:
                node.label = tag
        else:
            leaf = _Node()
            leaf.name = parse_name()
            skip_ws()
            leaf.length = parse_length()
            node = leaf
        return node

    root = parse_subtree()
    skip_ws()
    if pos < n and s[pos] == ";":
        pos += 1
    skip_ws()
    if pos != n:
        raise NewickParseError(f"trailing content at {pos}: {s[pos:pos + 20]!r}")
    return root


def _collect_leaf_sets(root: _Node) -> Tuple[Set[str], List[FrozenSet[str]]]:
    """Return (all tip names, leaf-set of every internal node incl. root)."""
    tips: Set[str] = set()
    clades: List[FrozenSet[str]] = []

    def walk(node: _Node) -> FrozenSet[str]:
        if node.is_leaf:
            if node.name:
                tips.add(node.name)
            return frozenset([node.name] if node.name else [])
        child_sets = [walk(c) for c in node.children]
        merged: FrozenSet[str] = frozenset().union(*child_sets)
        clades.append(merged)
        return merged

    walk(root)
    return tips, clades


# ---------------------------------------------------------------------------
# Public pure-Python measurements
# ---------------------------------------------------------------------------

def tip_set(newick: str) -> Set[str]:
    """All (non-empty) tip labels of the tree."""
    tips, _ = _collect_leaf_sets(_parse(newick))
    return tips


def clades_of(newick: str) -> Set[FrozenSet[str]]:
    """Leaf set of every internal node (root included) of the tree.

    A tip subset S is monophyletic (clade) iff ``frozenset(S)`` is in the
    returned set — same semantics as the ete3 MRCA-descendants check used by
    ``taxonomy.monophyly_proportion``.
    """
    _, clades = _collect_leaf_sets(_parse(newick))
    return set(clades)


def is_clade_monophyletic(newick: str, tip_subset: Set[str]) -> bool:
    """True iff ``tip_subset`` is exactly the leaf set of some internal node."""
    return frozenset(tip_subset) in clades_of(newick)


def splits_of(newick: str) -> Tuple[Set[FrozenSet[str]], Set[str]]:
    """Non-trivial (unrooted) splits of the tree and its tip set.

    For each internal node, one side of its edge is the node's leaf set S;
    the bipartition {S, T\\S} is kept when both sides have >= 2 tips (pendant
    edges carry no split) and normalized to a canonical representative so the
    two sides of the same edge dedupe (root children in particular).
    """
    root = _parse(newick)
    tips, clades = _collect_leaf_sets(root)
    splits: Set[FrozenSet[str]] = set()
    for clade in clades[:-1] if clades else []:
        # The last appended clade is the root's own leaf set == all tips;
        # Its edge (to a nonexistent parent) is not a split. Root children
        # Were appended before it and are kept.
        comp = tips - clade
        if len(clade) < 2 or len(comp) < 2:
            continue
        splits.add(_canonical_side(clade, comp))
    return splits, tips


def _canonical_side(
    side_a: Union[FrozenSet[str], Set[str]],
    side_b: Union[FrozenSet[str], Set[str]],
) -> FrozenSet[str]:
    sa, sb = sorted(side_a), sorted(side_b)
    return frozenset(side_a) if sa <= sb else frozenset(side_b)


def _restrict_splits(
    splits: Set[FrozenSet[str]], shared: Set[str]
) -> Set[FrozenSet[str]]:
    """Restrict splits to the shared tip set, dropping trivial bipartitions."""
    out: Set[FrozenSet[str]] = set()
    for split in splits:
        s = split & shared
        comp = shared - s
        if len(s) < 2 or len(comp) < 2:
            continue
        out.add(_canonical_side(s, comp))
    return out


def rf_distance(tree1_newick: str, tree2_newick: str) -> Tuple[int, float]:
    """Raw and normalized RF distance via pure-Python split sets.

    Only splits over the shared tip set are compared (matching ete3's
    ``robinson_foulds(unrooted_trees=True)`` over common leaves).
    ``max_rf = 2 * (n_shared - 3)``.

    Raises:
        TipCapExceeded: shared tips > MAX_PURE_PYTHON_TIPS.
        NewickParseError: either tree unparseable.
    """
    splits1, tips1 = splits_of(tree1_newick)
    splits2, tips2 = splits_of(tree2_newick)
    shared = tips1 & tips2
    if len(tips1) < 1 or len(tips2) < 1:
        raise NewickParseError("tree without tip labels")
    if len(shared) > MAX_PURE_PYTHON_TIPS:
        raise TipCapExceeded(
            f"{len(shared)} shared tips exceeds pure-Python cap "
            f"{MAX_PURE_PYTHON_TIPS} (reason=tip_cap_exceeded)"
        )
    if len(shared) < 3:
        # Fewer than 3 shared tips: every split is trivial, trees are
        # Uninformative; ete3 would report max_rf=0. Keep the (0, 0.0)
        # "identical on shared tips" convention? No — 0 shared internal
        # Splits means the comparison carries no information, but RF here is
        # Still well-defined as 0/0. Mirror ete3: max_rf=0 -> norm 0.0.
        return 0, 0.0
    r1 = _restrict_splits(splits1, shared)
    r2 = _restrict_splits(splits2, shared)
    rf = len(r1 ^ r2)
    max_rf = 2 * (len(shared) - 3)
    return rf, (rf / max_rf if max_rf > 0 else 0.0)


def quartet_topology(newick: str, quartet: Tuple[str, str, str, str]) -> Optional[str]:
    """Canonical induced topology of ``quartet`` on the tree.

    Returns ``((a,b),(c,d))`` — pairs sorted internally and between pairs —
    for a resolved 2-2 split, or ``(a,b,c,d)`` (sorted) for an unresolved
    (star) induced quartet. Rotation/order-of-writing invariant, so the
    ete3 path and this path produce identical strings and are directly
    comparable ( contract).

    Returns ``None`` when the tree cannot be parsed or a tip is missing.
    """
    try:
        splits, tips = splits_of(newick)
    except NewickParseError:
        return None
    q = frozenset(quartet)
    if len(q) != 4 or not q <= tips:
        return None
    ind = _restrict_splits(splits, set(q))
    for split in ind:
        if len(split) == 2:
            s1, s2 = sorted(split), sorted(q - split)
            pair_a, pair_b = sorted([tuple(s1), tuple(s2)])
            return f"(({pair_a[0]},{pair_a[1]}),({pair_b[0]},{pair_b[1]}))"
    names = ",".join(sorted(q))
    return f"({names})"


def canonical_small_newick(newick: str) -> Optional[str]:
    """Canonicalize an already-pruned 4-tip newick (ete3 helper).

    Feeds ete3's ``write(format=9)`` quartet output through the same
    canonical form as:func:`quartet_topology` so both measurement paths are
    comparable. Returns ``None`` when unparseable.
    """
    try:
        splits, tips = splits_of(newick)
    except NewickParseError:
        return None
    if len(tips) != 4:
        return None
    for split in splits:
        if len(split) == 2:
            s1, s2 = sorted(split), sorted(tips - split)
            pair_a, pair_b = sorted([tuple(s1), tuple(s2)])
            return f"(({pair_a[0]},{pair_a[1]}),({pair_b[0]},{pair_b[1]}))"
    return "(" + ",".join(sorted(tips)) + ")"
