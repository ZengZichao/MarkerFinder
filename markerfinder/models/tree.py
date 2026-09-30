"""Tree model -- a data shape that also carries derived measurements.

 (read this before treating ``Tree`` as a DTO): the methods below compute
topology facts (RF distance, quartet agreement, MAD rooting, support reading),
so a behaviour change in ``markerfinder.utils.etree`` or ``utils.tree_utils``
changes what this *model* reports, not just some service above it.

The reverse import inside method bodies (``utils.tree_utils`` importing this
module at module scope, this module importing them per call) is deliberate and
is the only bidirectional pair in the package: it keeps the import graph acyclic
at import time and keeps a single ete3 entry point, instead of leaving a
weaker second implementation of the same measurements inside the model.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


@dataclass
class Tree:
    """A gene/species tree plus the measurements derived from its topology."""

    newick: str = ""
    support_values: Dict[str, float] = field(default_factory=dict)
    # Reserved: no producer and no consumer anywhere in the package; kept as
    # declared API, annotated honestly.
    quartet_support: Optional[Dict[str, float]] = None
    # Which parser produced the latest get_tips result —
    # "ete3" | "regex". Numeric tip labels are preserved by both paths, but
    # The regex path never filtered support-like tokens, so decision cards and
    # Assertions need to know which path served the tips.
    tip_parse_mode: Optional[str] = None

    def get_tips(self) -> List[str]:
        # Prefer ete3 parsing so numeric taxon labels (e.g. "123", "12.5") are
        # Treated as legitimate tips rather than being dropped as if they were
        # Bootstrap/support values. The regex fallback is used only when ete3
        # Is unavailable and deliberately does NOT filter numeric names.
        try:
            # Single sanctioned ete3 entry point.
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            t = EteTree(self.newick, format=1)
            tips = [n.name for n in t.get_leaves() if n.name]
            self.tip_parse_mode = "ete3"
        except Exception:
            raw = re.findall(r"([A-Za-z0-9_.\-]+)(?=[,:)\]])", self.newick)
            tips = [t.strip() for t in raw if t.strip()]
            self.tip_parse_mode = "regex"
        # Warn (per call) when the source carries duplicate tip labels. The
        # De-duplication below silently drops them, which can mask a malformed
        # Input tree (e.g. the same taxon listed twice).
        if len(tips) != len(set(tips)):
            logger.warning(
                f"Tree.get_tips: input tree contains duplicate tip names "
                f"(kept {len(set(tips))} unique of {len(tips)}); the tree may "
                f"be malformed."
            )
        # Preserve order while de-duplicating.
        tips = list(dict.fromkeys(tips))
        return tips

    def to_newick(self) -> str:
        """Return the current Newick string backing this tree.

        Semantics are immutable: this method always returns ``self.newick``
        verbatim and does NOT re-serialize the in-memory tree state. The
        derived-measurement helpers (``get_induced_tree``,
        ``get_quartet_topology``) return a fresh value built from a copy, and
        ``None`` when that measurement is not available. Callers that mutate
        ``self`` must write the result back themselves; ``to_newick`` will not
        reflect unsaved edits.
        """
        return self.newick

    def get_induced_tree(self, species_set: Set[str]) -> Optional[Tree]:
        """Induced subtree on ``species_set``, or ``None`` when not measurable.

        the previous fallback returned ``self.newick`` verbatim, i.e. a
        FULL tree handed back as if it were the induced one. No caller could
        tell that apart from a real pruning result, and every downstream
        measurement on it would have looked plausible while being wrong.
        """
        try:
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            t = EteTree(self.newick)
            t.prune(list(species_set), preserve_branch_length=True)
            return Tree(newick=t.write(format=9))
        except Exception as e:
            logger.warning(
                f"Induced tree on {len(species_set)} taxon(a) not measurable: {e}"
            )
            return None

    def get_quartet_topology(
        self, quartet: Tuple[str, str, str, str]
    ) -> Optional[str]:
        """Induced quartet topology in canonical form, or ``None``.

        Delegates to:func:`markerfinder.utils.tree_utils.get_quartet_topology`
        (ete3 first, then the pure-Python split-set measurement) instead of
        keeping a second, weaker copy of the same measurement here. The old
        inline version returned ``""`` on any failure, which a caller cannot
        distinguish from "this quartet has no topology".
        """
        from markerfinder.utils.tree_utils import get_quartet_topology

        return get_quartet_topology(self.newick, quartet)

    def contains_bipartition(self, left: Set[str], right: Set[str]) -> Optional[bool]:
        """Check whether a bipartition exists in this tree.

        Returns ``None`` when the tree cannot be read: "not measured" is not the
        same statement as "the bipartition is absent".

.. note:: Reserved API — not called by the main pipeline. Kept for
           ad-hoc tree-set analysis and future use.
        """
        try:
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            t = EteTree(self.newick)
            for node in t.traverse("postorder"):
                if not node.is_leaf():
                    desc = set(n.name for n in node.get_leaves())
                    all_tips = set(n.name for n in t.get_leaves())
                    if desc == left or (all_tips - desc) == left:
                        comp = all_tips - desc
                        if comp == right or desc == right:
                            return True
            return False
        except Exception as e:
            logger.warning(f"Bipartition check not measurable: {e}")
            return None

    def get_mrca_depth(self, species_set: Set[str]) -> Optional[int]:
        """Depth of the MRCA of a species set from the root.

.. note:: Reserved API — not called by the main pipeline. Kept for
           ad-hoc tree-set analysis and future use.
        """
        try:
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            t = EteTree(self.newick)
            mrca = t.get_common_ancestor(list(species_set))
            depth = 0
            node = mrca
            while node.up is not None:
                depth += 1
                node = node.up
            return depth
        except Exception:
            # An MRCA depth we could not compute is NOT_MEASURABLE, not
            # Depth 0 (which would silently read as "the two taxa meet at the
            # Root", a real and very different biological statement).
            return None

    def are_sister_groups(self, left: Set[str], right: Set[str]) -> Optional[bool]:
        """Check whether two clades are sister groups.

        ``None`` means the question could not be asked of this tree; ``False``
        means it was asked and the clades are not sisters.

.. note:: Reserved API — not called by the main pipeline. Kept for
           ad-hoc tree-set analysis and future use.
        """
        try:
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            t = EteTree(self.newick)
            mrca_left = t.get_common_ancestor(list(left))
            mrca_right = t.get_common_ancestor(list(right))
            return mrca_left.up is not None and mrca_left.up is mrca_right.up
        except Exception as e:
            logger.warning(f"Sister-group check not measurable: {e}")
            return None

    @property
    def internal_branches(self) -> Optional[List[Dict[str, Any]]]:
        """Internal branches with support values and clade tip sets, or ``None``.

        An empty list asserts "this tree has no internal branches"; a tree that
        could not be read yields ``None`` instead.

.. note:: Reserved API — not called by the main pipeline. Kept for
           ad-hoc tree-set analysis and future use.
        """
        try:
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            t = EteTree(self.newick, format=1)
            branches = []
            for node in t.traverse("postorder"):
                if not node.is_leaf() and node.up is not None:
                    label = (node.name or "").strip()
                    try:
                        support: Optional[float] = float(label)
                    except ValueError:
                        # No numeric label means no support value. Under the
                        # Default ete3 parse this branch would have reported
                        # 0.0 or the fabricated 1.0 default instead.
                        support = None
                    branches.append({
                        "support": support,
                        "clade": set(n.name for n in node.get_leaves()),
                    })
            return branches
        except Exception as e:
            logger.warning(f"Internal-branch listing not measurable: {e}")
            return None

    @property
    def n_tips(self) -> int:
        return len(self.get_tips())

    def support_labels(self) -> List[float]:
        """Internal-node support labels actually present in the Newick string.

        Root and leaves are excluded, and — the part that matters — a node with
        **no label contributes nothing**. ete3's ``TreeNode.support`` attribute
        defaults to ``1.0`` for unlabelled nodes, so reading ``node.support``
        turns "this tree carries no support information" into "every branch has
        support 1.0", i.e. the maximum under the 0-1 convention FastTree and
        ASTRAL use. Measured with ete3 3.1.3::

            ((A,B),(C,D)); -> node.support == 1.0 (fabricated)
            ((A:1,B:1)98:1,(C:1,D:1)76:1) -> 98.0 / 76.0 (native percent)

        Parsing with ``format=1`` keeps the label text in ``node.name`` instead
        of substituting the default, so presence is decided from the evidence
. Values are returned in whatever scale the tree itself
        uses; callers that need a 0-1 fraction must detect the scale explicitly
        and refuse when it is ambiguous.
        """
        try:
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            tree = EteTree(self.newick, format=1)
        except Exception as e:
            logger.warning(f"Support labels not measurable: {e}")
            return []
        labels: List[float] = []
        for node in tree.traverse("postorder"):
            if node.is_leaf() or node.up is None:
                continue
            name = (node.name or "").strip()
            if not name:
                continue
            try:
                labels.append(float(name))
            except ValueError:
                # Clade identifiers such as "Node1" are not support values.
                continue
        return labels

    def get_average_support(self) -> Optional[float]:
        if self.support_values:
            vals = list(self.support_values.values())
            return sum(vals) / len(vals) if vals else 0.0
        labels = self.support_labels()
        # A tree that carries no readable internal support values has an
        # Undefined mean support, not a mean support of 0.0 (which reads as
        # "every branch is unsupported" — ) and not ete3's default 1.0
        # ("maximum support"). Native scale; see support_labels.
        return sum(labels) / len(labels) if labels else None

    @staticmethod
    def read(path: str) -> Tree:
        with open(path, encoding="utf-8", newline="") as f:
            return Tree(newick=f.read().strip())


@dataclass
class TreeRecommendation:
    recommended_tree: Optional[Tree] = None
    confidence: str = "unknown"
    reason: str = ""
    warnings: List[str] = field(default_factory=list)
    flagged_branches: List[Any] = field(default_factory=list)


@dataclass
class TreeVisualization:
    trees: Dict[str, Any] = field(default_factory=dict)
    default_layout: str = "rectangular"
    metadata_fields: List[str] = field(default_factory=list)
