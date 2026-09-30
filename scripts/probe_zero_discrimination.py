"""Zero-discrimination probe: does the output move when the input moves?

Feeds two 4-taxon trees — one maximally conflicted with the reference
``((A,B),(C,D));``, one identical to it — through
``PhylogeneticHGTDetector.detect`` + ``HGTDecisionEngine.evaluate_marker``.

Zero discrimination means the two inputs produce byte-identical output
(``normalized_rf=0.5``, ``quartet_consistency=1.0``, ``overall_risk=0.25``,
``level_2``, ``confidence=high``, empty notes). With the split-set measurement
and the unmeasurable→UNKNOWN / provenance-notes path, the two inputs MUST
produce different risks and non-empty notes. This script is also the
implementation body of assertion A-12.

Exit code 0 = discriminating; 1 = zero discrimination / missing provenance.
"""

import logging
import sys

from markerfinder.config import HGTConfig
from markerfinder.models.tree import Tree
from markerfinder.modules.hgt_filter import (
    HGTDecisionEngine,
    PhylogeneticHGTDetector,
)

logging.disable(logging.CRITICAL)


def probe(gene_tree: Tree, ref_tree: Tree, label: str):
    cfg = HGTConfig()
    detector = PhylogeneticHGTDetector(cfg)
    engine = HGTDecisionEngine(cfg)
    res = detector.detect(f"M_{label}", gene_tree, ref_tree)
    ev = engine.evaluate_marker(f"M_{label}", res)
    return (
        res.rf_distance,
        res.normalized_rf,
        res.quartet_consistency,
        res.overall_risk,
        ev.level.value,
        ev.confidence,
        ev.notes,
    )


def main() -> int:
    ref = Tree(newick="((A,B),(C,D));")
    clash = Tree(newick="((A,C),(B,D));")  # Maximally conflicting with ref

    a = probe(clash, ref, "clash")  # Expect: high risk / LEVEL_3 / notes
    b = probe(ref, ref, "same")     # Expect: low risk / LEVEL_1 / notes

    print(f"clash: rf={a[0]} norm_rf={a[1]} q={a[2]} risk={a[3]} "
          f"level={a[4]} conf={a[5]} notes={a[6]!r}")
    print(f"same : rf={b[0]} norm_rf={b[1]} q={b[2]} risk={b[3]} "
          f"level={b[4]} conf={b[5]} notes={b[6]!r}")

    ok = True
    if a[3] == b[3]:
        print("FAIL: zero discrimination — identical overall_risk for "
              "conflicting vs identical inputs (E0, §1.3)")
        ok = False
    if not str(a[6]).strip() or not str(b[6]).strip():
        print("FAIL: provenance missing — notes empty for a graded marker "
              "(the assertion layer's own probe)")
        ok = False
    if ok:
        print("PASS: outputs discriminate (A-12 / §8.4)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
