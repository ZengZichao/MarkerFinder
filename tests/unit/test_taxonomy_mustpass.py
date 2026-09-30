""" +: the taxonomy must-pass gate must actually gate.

The must-pass baseline existed only as a YAML file with no
reader anywhere in the package, so "the pipeline is rejected and aborts with
exit 4 when a required relationship fails" could never happen. These tests
plant a real violation and require the abort, and require that an
unloadable / unsupported rule is reported rather than silently passing.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from markerfinder.cli.constants import EXIT_ASSERTION_FAILED
from markerfinder.exceptions import AssertionFailureError
from markerfinder.modules.taxonomy_mustpass import (
    check_taxonomy_mustpass,
    enforce_mustpass,
    load_mustpass,
    write_mustpass_tsv,
)

TAXONOMY = {
    "A1": {"genus": "Alpha", "phylum": "P1"},
    "A2": {"genus": "Alpha", "phylum": "P1"},
    "B1": {"genus": "Beta", "phylum": "P1"},
    "B2": {"genus": "Beta", "phylum": "P1"},
}
RULES = [
    {"relation": "monophyletic", "level": "genus", "min_taxa": 2},
    {"relation": "monophyletic", "level": "phylum", "min_taxa": 2},
]
MONOPHYLETIC = "((A1,A2),(B1,B2));"
# Alpha and Beta each scattered across both sides: taxon monophyly is broken.
POLYPHYLETIC = "((A1,B1),(A2,B2));"


class TestGateBehaviour:
    def test_clean_tree_passes(self):
        report = check_taxonomy_mustpass(
            {"M1": MONOPHYLETIC}, TAXONOMY, RULES
        )
        assert report.passed
        assert report.groups_tested == 3  # Alpha, Beta, P1

    def test_planted_polyphyly_aborts_with_exit_4_semantics(self):
        report = check_taxonomy_mustpass(
            {"M1": POLYPHYLETIC}, TAXONOMY, RULES
        )
        assert not report.passed
        assert {v.taxon for v in report.violations} == {"Alpha", "Beta"}
        with pytest.raises(AssertionFailureError) as exc:
            enforce_mustpass(report)
        assert "must-pass" in str(exc.value)
        assert EXIT_ASSERTION_FAILED == 4

    def test_no_taxonomy_map_is_not_a_pass(self):
        report = check_taxonomy_mustpass({"M1": POLYPHYLETIC}, {}, RULES)
        assert report.not_checked
        assert "no taxonomy map supplied" in report.not_checked[0]

    def test_unsupported_relation_is_never_silently_passed(self):
        report = check_taxonomy_mustpass(
            {"M1": MONOPHYLETIC}, TAXONOMY,
            [{"relation": "sister_of", "level": "genus"}],
        )
        assert any("unsupported relation" in line for line in report.not_checked)

    def test_partial_sampling_is_reported_not_invented(self):
        trees = {"M1": "((A1,B1),(A2,B2));"}
        partial = {k: v for k, v in TAXONOMY.items()}
        partial["ghost"] = {"genus": "Gamma"}  # A taxon absent from the tree
        report = check_taxonomy_mustpass(trees, partial, RULES)
        # Gamma never appears in the tree, so it must not be counted as tested.
        assert all(v.taxon != "Gamma" for v in report.violations)

    def test_unreadable_tree_is_recorded(self):
        report = check_taxonomy_mustpass({"M1": "((A1,A2"}, TAXONOMY, RULES)
        assert any("unparseable" in line for line in report.not_checked)


class TestBaselineLoading:
    def test_shipped_skeleton_loads_two_rules(self):
        rules = load_mustpass(
            "tests/benchmark/expected/taxonomy_mustpass.yaml"
        )
        assert len(rules) == 2
        assert {r["level"] for r in rules} == {"genus", "phylum"}

    def test_missing_file_yields_no_rules(self):
        assert load_mustpass("nope/does-not-exist.yaml") == []
        assert load_mustpass(None) == []


class TestProduct:
    def test_verdict_is_persisted(self, tmp_path):
        report = check_taxonomy_mustpass({"M1": POLYPHYLETIC}, TAXONOMY, RULES)
        path = write_mustpass_tsv(report, str(tmp_path), "mf")
        assert path == tmp_path / "Phase5_reports" / "mf.mustpass.tsv"
        body = path.read_text(encoding="utf-8")
        assert body.splitlines()[0].startswith("kind\tmarker_id")
        assert "summary" in body and "FAIL" in body
        assert "violation" in body and "Alpha" in body


class TestWiringLocks:
    def test_cli_exposes_the_flag(self):
        import markerfinder.cli.parser as parser_module

        src = inspect.getsource(parser_module)
        assert "--taxonomy-mustpass" in src

    def test_config_and_pipeline_consumed_it(self):
        import markerfinder.cli.config_build as config_build
        import markerfinder.pipeline as pipeline_module

        assert "taxonomy_mustpass=" in inspect.getsource(config_build)
        pipeline_src = inspect.getsource(pipeline_module)
        # The pipeline now drives the gate through the stage's own entry
        # Point (check + requirements-file stamping + TSV write in one call), so
        # The lock follows the call it actually makes. Intent is unchanged: the
        # Gate must run, must record which file it read, and must be able to abort.
        assert "evaluate_and_write(" in pipeline_src, (
            "must-pass gate is never evaluated by the run"
        )
        assert "write_mustpass_tsv(" not in pipeline_src, (
            "the orchestrator writes the gate's artifact itself again -- see "
            "tests/unit/test_stage_owns_its_products.py"
        )
        assert "enforce_mustpass(" in pipeline_src, (
            "gate computes a verdict but never aborts"
        )
        assert "must-pass gate NOT ENABLED" in pipeline_src, (
            "an unset gate must say so, not stay silent"
        )

    def test_repo_relative_default_resolves(self):
        # A bare --taxonomy-mustpass must not depend on the caller's cwd, or
        # The gate would abort for a packaging reason instead of a biological
        # One.
        from markerfinder.cli.constants import default_taxonomy_mustpass_path

        import markerfinder.cli.parser as parser_module

        resolved = Path(default_taxonomy_mustpass_path())
        assert resolved.is_absolute() and resolved.exists(), resolved
        assert "default_taxonomy_mustpass_path()" in inspect.getsource(
            parser_module
        )


class TestBaselineScope:
    def test_unfilled_skeleton_means_no_restriction(self):
        from markerfinder.modules.taxonomy_mustpass import mustpass_marker_scope

        assert mustpass_marker_scope(
            "tests/benchmark/expected/taxonomy_mustpass.yaml"
        ) is None
        assert mustpass_marker_scope(None) is None

    def test_filled_ids_are_honoured(self, tmp_path):
        from markerfinder.modules.taxonomy_mustpass import mustpass_marker_scope

        f = tmp_path / "mp.yaml"
        f.write_text(
            "markers:\n  ids:\n    - UBA_P07473\n    - UBA_P07483\n"
            "must_pass:\n  - relation: monophyletic\n    level: genus\n",
            encoding="utf-8",
        )
        assert mustpass_marker_scope(str(f)) == ["UBA_P07473", "UBA_P07483"]
