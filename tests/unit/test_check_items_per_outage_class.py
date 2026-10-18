""" +: failures must be loud, distinguishable, and self-checked.

"依赖缺失、工具缺失、参照非法必须产生可区分的退出码与报告条目". Three
different outages must not collapse into one generic failure, and each must leave
an artifact a reader can act on.

 additionally asked for a ``--check`` item for the composition diagnostics.
That item did not exist: ``grep composition markerfinder/cli/self_test.py``
returned nothing, while the requirement names it. It is implemented in
``self_test._test_composition`` and the controls below make sure it is not a
rubber stamp -- flattening the metric must turn the item red.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from markerfinder.cli import constants
from markerfinder.cli import self_test


# ──: distinguishable failure surface ───────────────────────────────

def test_control_exit_codes_are_pairwise_distinct():
    codes = {
        name: value
        for name, value in vars(constants).items()
        if name.startswith("EXIT_") and isinstance(value, int)
    }
    assert codes, "no exit-code constants found"
    assert len(set(codes.values())) == len(codes), codes
    assert codes["EXIT_ASSERTION_FAILED"] == 4
    assert codes["EXIT_INCONCLUSIVE"] == 5
    assert codes["EXIT_DATA_ERROR"] == 3


def test_missing_external_tools_are_named_individually(monkeypatch):
    """A missing aligner must be named tool-by-tool, never a generic failure."""
    import shutil as shutil_module

    monkeypatch.setattr(shutil_module, "which", lambda name: None)
    from markerfinder.utils.dependency_check import check_external_tools

    status = check_external_tools()
    assert status, "no tools were checked at all"
    assert not any(status.values()), status
    for tool in ("mafft", "trimal"):
        assert tool in status, sorted(status)


def test_precheck_logs_the_tool_gap_as_its_own_entry(caplog, monkeypatch):
    """ Wants a report entry per outage class: tools missing must be
    distinguishable in the log from a broken python dependency.

    The assertion must not depend on what happens to be installed on the
    machine running the suite, so the tool-absent case is produced by
    simulating an empty PATH rather than by hoping the machine is bare.
    """
    import logging
    import shutil as shutil_module

    from markerfinder.utils.dependency_check import precheck_all

    monkeypatch.setattr(shutil_module, "which", lambda _name: None)
    with caplog.at_level(logging.DEBUG, logger="markerfinder.utils.dependency_check"):
        precheck_all()
    text = caplog.text.lower()
    assert "external tools" in text, text[-800:]
    assert "not found" in text, text[-800:]


def test_precheck_names_the_tool_class_even_when_everything_is_present(caplog):
    """'No log line at all' used to be ambiguous between 'all tools present'
    and 'the tool check never ran'. The class must be reported either way."""
    import logging

    from markerfinder.utils.dependency_check import precheck_all

    with caplog.at_level(logging.DEBUG, logger="markerfinder.utils.dependency_check"):
        precheck_all()
    assert "external tools" in caplog.text.lower(), caplog.text[-800:]


def test_tool_gap_and_python_gap_are_not_the_same_entry(caplog, monkeypatch):
    """The two outage classes must be separable by a reader of the log."""
    import logging
    import shutil as shutil_module

    import markerfinder.utils.etree as etree_module
    from markerfinder.utils.dependency_check import precheck_all

    def _broken():
        raise RuntimeError("No module named 'cgi'")

    monkeypatch.setattr(etree_module, "require_ete3", _broken)
    monkeypatch.setattr(shutil_module, "which", lambda _name: None)
    with caplog.at_level(logging.DEBUG, logger="markerfinder.utils.dependency_check"):
        precheck_all()
    text = caplog.text.lower()
    assert "external tools" in text, text[-800:]
    assert "ete3 is unusable" in text, text[-800:]
    # One generic catch-all line would make the two outages indistinguishable.
    assert "dependency problem" not in text, text[-800:]


def test_broken_dependency_names_the_real_cause_not_a_false_diagnosis(monkeypatch):
    """`ete3 is not installed` was printed when ete3 WAS installed but could not
    be imported. Those are different repairs, and the message must name the
    real one."""
    import markerfinder.utils.etree as etree_module
    from markerfinder.exceptions import PhyloToolUnavailable

    def _broken():
        raise PhyloToolUnavailable(
            "ete3 is not importable in this interpreter: No module named 'cgi'"
        )

    monkeypatch.setattr(etree_module, "require_ete3", _broken)
    from markerfinder.utils.dependency_check import check_python_dependencies

    messages = check_python_dependencies()
    text = "\n".join(messages)
    assert "cgi" in text, messages
    assert "is not installed" not in text, (
        "the message must not assert a diagnosis the interpreter contradicts"
    )


def test_illegal_reference_fails_as_a_named_state_not_a_generic_error():
    """The third outage class: an unusable reference must carry its own
    code, so a reader can tell 'reference merged the focal clades' from
    'reference has three tips' without re-running anything."""
    from pathlib import Path as _P

    from markerfinder.utils.reference import (
        CLADE_COLLAPSED,
        MALFORMED,
        TOO_FEW_TIPS,
        validate_reference_tree,
    )

    fixtures = _P(__file__).resolve().parents[1] / "fixtures" / "reference"

    def _newick(name: str) -> str:
        # Validate_reference_tree takes Newick text, not a path: passing a path
        # Silently produced TOO_FEW_TIPS on the first run of this test.
        return (fixtures / name).read_text(encoding="utf-8").strip()

    collapsed = validate_reference_tree(
        _newick("collapsed_sponge_other.nwk"),
        _newick("collapsed_sponge_other.nwk"),
        target_clades=[{"S1", "S2"}, {"O1"}],
    )
    assert not collapsed.ok
    assert collapsed.code == CLADE_COLLAPSED

    malformed = validate_reference_tree("((A,B;", "(A,B);")
    assert not malformed.ok and malformed.code == MALFORMED

    tiny = validate_reference_tree(_newick("three_tip_ref.nwk"), "(A,B);")
    assert not tiny.ok and tiny.code == TOO_FEW_TIPS

    assert len({collapsed.code, malformed.code, tiny.code}) == 3, (
        "three different outages collapsed onto one code"
    )


def _item(results, prefix):
    matches = [row for row in results if row[0].startswith(prefix)]
    assert matches, (prefix, [r[0] for r in results])
    return matches[0]


def test_control_composition_self_test_passes_on_the_real_metric():
    results = self_test._test_composition()
    assert len(results) == 3, results
    assert _item(results, "Composition RCV discriminates")[1] == "PASS"
    assert _item(results, "Composition: protein GC bias is NA")[1] == "PASS"
    assert _item(results, "Composition: empty input")[1] == "PASS"


def test_composition_self_test_turns_red_when_the_metric_goes_flat(monkeypatch):
    """Must-fail control: if the check cannot see a broken metric it is a stamp."""
    import markerfinder.modules.composition as composition

    monkeypatch.setattr(composition, "rcv", lambda seqs, group_labels=None: 0.0)
    results = self_test._test_composition()
    name, status, detail = _item(results, "Composition RCV discriminates")
    assert status == "FAIL", (name, status, detail)
    assert "flat" in detail, detail


def test_composition_self_test_turns_red_when_na_becomes_a_number(monkeypatch):
    import markerfinder.modules.composition as composition

    monkeypatch.setattr(
        composition, "gc_or_codon_bias", lambda seqs, protein=True: 0.0
    )
    results = self_test._test_composition()
    name, status, detail = _item(results, "Composition: protein GC bias is NA")
    assert status == "FAIL", (name, status, detail)
    assert "placeholder" in detail, detail


def test_composition_self_test_is_registered_in_the_check_table():
    source = Path(self_test.__file__).read_text(encoding="utf-8")
    assert "results += _test_composition()" in source, (
        "an unregistered self-test is invisible to --check, which is how the "
        "second half went unimplemented for a whole round"
    )
