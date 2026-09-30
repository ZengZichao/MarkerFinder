"""The consistency strength metric, and the two things it must not do.

 is short -- "强度指标用已有的 ``calculate_rf_distance`` 与
``_calculate_quartet_consistency`` 取绝对值差形式，不引入似然计算、不新增外部工具" --
but every clause is a ban, and a ban nobody asserts
gradually stops holding. The formula comment was in the code; no test named
or checked the metric's shape.

Measured properties here: identical trees must score exactly 1.0, the metric must
be symmetric (that is what "absolute difference" buys), a disagreeing pair must
score strictly lower (a metric that cannot discriminate is zero-discrimination again), and the two
prohibitions are locked by a scanner that carries its own positive control, so a
silently-empty scan cannot pass.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from markerfinder.modules.consistency_screen import _stance

MODULE_ROOT = Path(__file__).resolve().parents[2] / "markerfinder"

IDENTICAL = "((A:1,B:1):1,(C:1,D:1):1);"
DISAGREEING = "((A:1,C:1):1,(B:1,D:1):1);"
TWO_TIP = "(A:1,B:1);"
DUPLICATE_TIP = "((A:1,A:1):1,(C:1,D:1):1);"


def test_control_identical_trees_score_exactly_one():
    agrees, measurement = _stance(IDENTICAL, IDENTICAL, 0.3, 0.7)
    assert measurement.reason == ""
    assert agrees is True
    assert measurement.value == pytest.approx(1.0), measurement
    # The observation base survives: a 4-taxon pair has exactly one quartet, and
    # The product must be able to say so instead of a bare agreement ratio.
    assert measurement.n_obs == 1, measurement


def test_control_disagreement_scores_strictly_lower():
    """The discrimination property: a constant metric would pass every other test."""
    _, same = _stance(IDENTICAL, IDENTICAL, 0.3, 0.7)
    _, different = _stance(IDENTICAL, DISAGREEING, 0.3, 0.7)
    assert different.value < same.value, (different, same)


def test_strength_is_symmetric_in_the_two_trees():
    """Absolute-difference form: swapping gene and reference cannot change it."""
    _, forward = _stance(IDENTICAL, DISAGREEING, 0.3, 0.7)
    _, backward = _stance(DISAGREEING, IDENTICAL, 0.3, 0.7)
    assert forward.value == pytest.approx(backward.value)


def test_bounds_decide_agreement_in_both_directions():
    agrees_loose, _ = _stance(IDENTICAL, DISAGREEING, 1.0, 0.0)
    agrees_strict, _ = _stance(IDENTICAL, DISAGREEING, 0.0, 1.0)
    assert agrees_loose is True
    assert agrees_strict is False


def test_unmeasurable_pairs_report_a_reason_instead_of_a_zero_strength():
    cases = (
        ((TWO_TIP, TWO_TIP), "reference_illegal"),
        # A duplicated tip name passes the reference legality check but shrinks
        # The shared tip set below four, so the quartet leg cannot be evaluated
        # At all -- measured, not assumed: the validator lets it through and the
        # Abstention happens one step later.
        ((DUPLICATE_TIP, IDENTICAL), "fewer than 4 shared tips"),
    )
    for (gene, ref), expected in cases:
        agrees, measurement = _stance(gene, ref, 0.3, 0.7)
        assert agrees is None, (gene, ref)
        # NOT_MEASURABLE carries no number: strength stays None rather than 0.0
        assert measurement.value is None, measurement
        assert expected in measurement.reason, measurement
    # Two-sided control: a measurable leg does produce a number
    agrees, measurable = _stance(IDENTICAL, IDENTICAL, 0.3, 0.7)
    assert agrees is True and measurable.value is not None


def _identifiers(path: Path) -> set:
    """Names the code actually uses: imports, call targets, attribute lookups.

    Deliberately NOT a substring search over the file: the module's own comment
    says "no likelihood computation", so grepping the raw text would either fail
    on the prohibition it is enforcing or need an exception list.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            names.add(node.id.lower())
        elif isinstance(node, ast.Attribute):
            names.add(node.attr.lower())
        elif isinstance(node, ast.Import):
            names.update(a.name.lower() for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.lower())
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            names.add(node.func.id.lower())
    return names


def test_module_never_uses_likelihood_or_a_new_tool():
    """AD's two prohibitions, with a positive control for the probe.

    A banned word must be absent here and a sought word must be present in a
    module that really does shell out; otherwise "I grepped and found nothing" is
    indistinguishable from a broken probe.
    """
    names = _identifiers(MODULE_ROOT / "modules" / "consistency_screen.py")
    for banned in ("lnl", "likelihood", "subprocess", "popen", "iqtree", "raxml"):
        assert not any(banned in name for name in names), (banned, sorted(names))

    control = _identifiers(MODULE_ROOT / "modules" / "phylogenetic_inference.py")
    assert any("subprocess" in name for name in control), (
        "the probe found no subprocess anywhere, so the negative lock above would "
        "be vacuous -- fix the probe before trusting it"
    )
