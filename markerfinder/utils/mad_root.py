"""MAD rooting for gene trees (Minimal Ancestor Deviation).

Implements the original MAD (Minimal Ancestor Deviation) criterion of Tria
et al. (2017, *Nature Ecology & Evolution*) — note this is the original MAD
formulation, not a separate "2.0" variant. The optimal root position minimizes
the sum of squared
"ancestor deviations" over all internal nodes v:

    δ(v) = dist(root, v) − mean(dist(root, children of v))
    score(root) = Σ_v δ(v)²

Because δ(v) is linear in the root position along any single edge, the score
along an edge is a quadratic function of the interior root position, so the
unconstrained minimum is available in closed form and clamped to the edge.
The global optimum is taken over every internal node and every edge of the
(unrooted) tree.

MAD only needs branch lengths; when they are absent (or the tree is too large
to root cheaply) the tree is midpoint-rooted instead. The root position does
not change monophyly, but rooting with MAD is the method MarkerFinder uses to
prepare per-marker gene trees before the single-rank monophyly-proportion HGT
screen.
"""

from __future__ import annotations

import logging
from typing import List, Optional, Tuple

from markerfinder.exceptions import PhyloToolUnavailable

logger = logging.getLogger(__name__)

# Above this many tips, MAD rooting is skipped (O(n^2) candidates) and we fall
# Back to midpoint rooting, which is sufficient for the root-invariant
# Monophyly screen downstream.
_MAX_MAD_TIPS = 600


def mad_root(tree_newick: str) -> str:
    """Root a (unrooted) gene tree using the MAD criterion.

    Args:
        tree_newick: Newick string (branch lengths optional).

    Returns:
        Rerooted Newick string.

    Raises:
        PhyloToolUnavailable: when ete3 cannot be imported. MAD rooting has no
        pure-Python equivalent ( deliberately covers only the monophyly
        *test*, not the rooting), and the previous behaviour was to hand back the
        input tree unchanged with a single warning line — a caller that then
        measured clades on it had no way to learn that it was measuring an
        unrooted tree. "Cannot root" is a different statement from "this is the
        root", so it is now an exception ('s one failure type).

    Note:
        Two fallbacks remain and both are logged, but neither is a MAD result:
        oversized trees (> ``_MAX_MAD_TIPS``) and trees without branch lengths
        fall back to midpoint rooting. On a lengthless tree ete3 supplies unit
        branch lengths while parsing, so the "midpoint" it finds is arbitrary and
        the returned Newick carries lengths the input never had. See the closure
        report (open defect O).
    """
    from markerfinder.utils.etree import require_ete3

    try:
        EteTree = require_ete3().Tree
    except Exception as e:
        raise PhyloToolUnavailable(
            f"MAD rooting requires ete3, which is unusable here ({e}); "
            "there is no pure-Python rooting fallback."
        ) from e

    try:
        t = EteTree(tree_newick, format=1)
    except Exception as e:  # Malformed input
        logger.warning(f"mad_root: failed to parse tree: {e}")
        return tree_newick

    leaves = [n for n in t.get_leaves() if n.name]
    if len(leaves) < 3:
        return tree_newick

    if len(leaves) > _MAX_MAD_TIPS:
        logger.warning(
            f"mad_root: {len(leaves)} tips > {_MAX_MAD_TIPS}; MAD search "
            "skipped, rooting method = midpoint (not MAD)."
        )
        try:
            t.set_outgroup(t.get_midpoint_outgroup())
        except Exception:
            pass
        return t.write(format=1)

    if not _has_branch_lengths(t):
        logger.warning(
            "mad_root: tree carries no branch lengths, so MAD is undefined; "
            "falling back to an ARBITRARY midpoint root (ete3 invents unit "
            "branch lengths during parsing). Rooting method = midpoint/"
            "arbitrary, not MAD."
        )
        try:
            t.set_outgroup(t.get_midpoint_outgroup())
        except Exception:
            pass
        return t.write(format=1)

    nodes = list(t.traverse())
    internal = [n for n in nodes if not n.is_leaf()]

    # Pre-compute the leaf-clade sets for every original edge. These let us
    # Re-locate each edge on the (possibly re-rooted) working tree and on a
    # Fresh copy of the original tree by topology alone, eliminating the
    # Per-candidate ``t.copy`` that made this O(n^2) in memory.
    edge_clades = []
    for (u, v) in _get_edges(t):
        cu = frozenset(n.name for n in u.get_leaves() if n.name)
        cv = frozenset(n.name for n in v.get_leaves() if n.name)
        edge_clades.append((cu, cv))

    # Single pristine copy of the original topology, used only to materialize
    # The best rooting at the end (O(n) memory instead of O(n^2)).
    orig_t = t.copy()

    best_score: Optional[float] = None
    best_is_node = False
    best_node_clade: Optional[frozenset] = None
    best_edge_u: Optional[frozenset] = None
    best_edge_v: Optional[frozenset] = None
    best_t_star: float = 0.0

    # Candidate 1: root at each internal node (move the root on the working
    # Tree instead of copying it).
    for r in internal:
        try:
            t.set_outgroup(r)
        except Exception:
            continue
        s = _score_rooted(t)
        if best_score is None or s < best_score - 1e-12:
            best_score = s
            best_is_node = True
            # Record the outgroup as the FROZEN SET of its descendant leaf names
            # (not the node object, which is invalid on the fresh ``orig_t`` copy
            # Used below). ete3 3.1.3's ``set_outgroup`` rejects a bare list of
            # Leaf-name strings; the matching node is re-located by topology via
            # ``_node_for_clade``.
            best_node_clade = frozenset(n.name for n in r.get_leaves() if n.name)

    # Candidate 2: interior optimum along each edge (closed-form t*).
    for (cu, cv) in edge_clades:
        u_node, v_node = _find_edge_by_clades(t, cu, cv)
        if u_node is None or v_node is None:
            continue
        L = u_node.get_distance(v_node)
        if L <= 0:
            continue
        deg_u = len(_neighbors(u_node))
        deg_v = len(_neighbors(v_node))
        if deg_u == 0 or deg_v == 0:
            continue
        S = sum(u_node.get_distance(c) for c in _neighbors(u_node) if c is not v_node)
        T = sum(v_node.get_distance(c) for c in _neighbors(v_node) if c is not u_node)
        p1 = float(deg_u * deg_u)
        p2 = float(deg_v * deg_v)
        if p1 + p2 == 0:
            continue
        A = L + S
        B = T
        # Unconstrained minimum of the quadratic score(t) along the edge.
        t_star = (p1 * (L - B) + p2 * A) / (2.0 * (p1 + p2))
        if t_star <= 0.0 or t_star >= L:
            # Optimum is at an endpoint -> already covered by node candidates
            continue
        try:
            _root_on_edge(t, u_node, v_node, t_star)
        except Exception:
            continue
        s = _score_rooted(t)
        if best_score is None or s < best_score - 1e-12:
            best_score = s
            best_is_node = False
            best_edge_u = cu
            best_edge_v = cv
            best_t_star = t_star

    if best_score is None:
        return t.write(format=1)

    # Materialize the single best rooting on a fresh copy of the original tree
    # (its edges are intact, so the closed-form placement is exact). This is the
    # Only extra copy — O(n) memory instead of O(n^2) deep copies.
    best_tree = orig_t.copy()
    if best_is_node:
        # ``best_tree`` is a *copy* of ``orig_t``; node identities differ from the
        # Working tree ``t``, so re-locate the outgroup node by its leaf-clade
        # Topology instead of reusing the old reference ``r``.
        clade_node = _node_for_clade(best_tree, best_node_clade)
        if clade_node is not None:
            best_tree.set_outgroup(clade_node)
    else:
        u_node, v_node = _find_edge_by_clades(best_tree, best_edge_u, best_edge_v)
        if u_node is not None and v_node is not None:
            _root_on_edge(best_tree, u_node, v_node, best_t_star)
    return best_tree.write(format=1)


# ── helpers ──

def _neighbors(node) -> List:
    """All adjacent nodes of *node* in the unrooted sense (parent + children)."""
    ns = []
    if node.up is not None:
        ns.append(node.up)
    ns.extend(node.get_children())
    return ns


def _get_edges(tree) -> List[Tuple]:
    """All undirected edges of the tree as (parent, child) pairs (each once)."""
    edges = []
    for node in tree.traverse():
        for ch in node.get_children():
            edges.append((node, ch))
    return edges


def _find_edge_by_clades(tree, clade_u: frozenset, clade_v: frozenset):
    """Locate the (parent, child) node pair whose two sides match clade_u/clade_v.

    Returns ``(parent, child)`` where ``child``'s descendant leaves equal one of
    the clades and the rest of the tree equals the other. Used to re-locate an
    original edge on a (possibly re-rooted) tree so the closed-form rooting can
    be replayed without re-copying the whole tree.
    """
    total = frozenset(n.name for n in tree.get_leaves() if n.name)
    for node in tree.traverse():
        for ch in node.get_children():
            desc = frozenset(n.name for n in ch.get_leaves() if n.name)
            rest = total - desc
            if desc == clade_u and rest == clade_v:
                return node, ch
            if desc == clade_v and rest == clade_u:
                return node, ch
    return None, None


def _node_for_clade(tree, clade_set: frozenset):
    """Locate the node in *tree* whose descendant leaf-name set equals *clade_set*.

    Used to re-position a stored outgroup target (recorded as a leaf-name
    frozenset) onto a *fresh copy* of the original tree, where node identities
    from the working tree are no longer valid. Returns the matching ete3
    TreeNode, or ``None`` when no node's leaf clade matches.

    A clade set defines a unique MRCA in a tree, so this is equivalent to
    ``tree.get_common_ancestor(list(clade_set))`` but safe across tree copies.
    """
    if not clade_set:
        return None
    for node in tree.traverse():
        desc = frozenset(n.name for n in node.get_leaves() if n.name)
        if desc == clade_set:
            return node
    return None


def _has_branch_lengths(tree) -> bool:
    """True if at least one branch has a positive length."""
    for node in tree.traverse():
        if node.up is not None and node.dist:
            try:
                if float(node.dist) > 0:
                    return True
            except (TypeError, ValueError):
                continue
    return False


def _score_rooted(root_node) -> float:
    """Σ δ(v)² for a rooted ete3 tree (root_node is the tree root).

    δ(v) = dist(root, v) − mean(dist(root, children of v)).
    """
    dist_from_root: dict = {}
    stack = [(root_node, 0.0)]
    order: List = []
    while stack:
        node, d = stack.pop()
        if id(node) in dist_from_root:
            continue
        dist_from_root[id(node)] = d
        order.append(node)
        for ch in node.get_children():
            bl = ch.dist if (ch.dist is not None) else 0.0
            try:
                bl = float(bl)
            except (TypeError, ValueError):
                bl = 0.0
            stack.append((ch, d + bl))

    total = 0.0
    for node in order:
        if node.is_leaf():
            continue
        children = node.get_children()
        if not children:
            continue
        d_node = dist_from_root[id(node)]
        child_dists = [dist_from_root[id(c)] for c in children]
        mean_c = sum(child_dists) / len(child_dists)
        delta = d_node - mean_c
        total += delta * delta
    return total


def _root_on_edge(tree, u, v, t: float) -> None:
    """Reroot *tree* at a point at distance *t* from *u* along edge (u, v)."""
    from markerfinder.utils.etree import require_ete3

    EteTree = require_ete3().Tree

    L = u.get_distance(v)
    if t <= 0.0:
        tree.set_outgroup(u)
        return
    if t >= L:
        tree.set_outgroup(v)
        return
    v.detach()
    mid = EteTree()
    u.add_child(mid, dist=t)
    mid.add_child(v, dist=L - t)
    tree.set_outgroup(mid)
