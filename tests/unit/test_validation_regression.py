"""Code-review regression tests for ``validation.py`` (#3 A8/A9/A10).

Locks in three changed behaviours of ``validate_newick_string``:
  - A8: scientific-notation negative branches (``:-1.5e-3``) are flagged.
  - A9: quoted labels (``'Taxon (A)'``) parse without false "unbalanced
        quotes" / empty-tip errors; duplicate internal names are warnings.
  - A10: the function returns a 2-tuple ``(errors, warnings)`` instead of a
         single list.
"""

import pytest

from markerfinder.validation import validate_newick_string


class TestScientificNotationNegativeBranch:
    def test_neg_e_notation_flags(self):
        errors, _ = validate_newick_string("((A,(B):-1.5e-3),C);")
        assert any("Negative" in e for e in errors)

    def test_neg_e_notation_internal(self):
        errors, _ = validate_newick_string("((A,B:-2.0e1),C);")
        assert any("Negative" in e for e in errors)

    def test_plain_negative_still_flags(self):
        errors, _ = validate_newick_string("((A:-1,B),(C,D));")
        assert any("Negative" in e for e in errors)


class TestQuotedLabels:
    def test_quoted_tips_no_false_error(self):
        errors, warnings = validate_newick_string("('Taxon (A)','Node 1');")
        assert errors == []
        assert warnings == []

    def test_quoted_internal_label(self):
        errors, _ = validate_newick_string("((A,B)'Internal (x)',C);")
        assert errors == []

    def test_quoted_tips_not_seen_as_empty(self):
        errors, _ = validate_newick_string("('Taxon (A)','Node 1');")
        assert not any("Empty tip" in e for e in errors)


class TestReturnSignature:
    def test_returns_two_tuple(self):
        res = validate_newick_string("((A,B),(C,D));")
        assert isinstance(res, tuple) and len(res) == 2
        errors, warnings = res
        assert isinstance(errors, list)
        assert isinstance(warnings, list)

    def test_empty_returns_two_tuple(self):
        errors, warnings = validate_newick_string("")
        assert isinstance(errors, list) and isinstance(warnings, list)
        assert len(errors) > 0


class TestDuplicateInternalWarning:
    def test_duplicate_internal_is_warning_not_error(self):
        errors, warnings = validate_newick_string("((A,B)X,(C,D)X);")
        assert not any("Duplicate internal" in e for e in errors)
        assert any("Duplicate internal" in w for w in warnings)
