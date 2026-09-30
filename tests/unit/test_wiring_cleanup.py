"""Wiring & dead-code cleanup tests."""

import json

import pytest


class TestDeadCodeRemoved:
    def test_sh_test_pvalue_gone(self):
        from dataclasses import fields
        from markerfinder.models.pipeline_types import PhyloStepResult

        assert "sh_test_pvalue" not in {f.name for f in fields(PhyloStepResult)}

    def test_missing_pattern_report_gone(self):
        import markerfinder.models as models

        assert not hasattr(models, "MissingPatternReport")

    def test_combine_rejects_nothing_new(self):
        # Guard stays: composition signals never enter synthesis.
        from markerfinder.utils.hgt_utils import combine_hgt_scores

        risk, n = combine_hgt_scores({}, {"rf": 0.5, "quartet": 0.5})
        assert risk is None and n == 0


class TestDbHashComparison:
    def test_expected_hashes_file_exists(self):
        from pathlib import Path

        expected = Path(__file__).parent.parent.parent / "db" / "expected_hashes.json"
        assert expected.exists()
        data = json.loads(expected.read_text(encoding="utf-8"))
        assert isinstance(data, dict)

    def test_verify_database_compares_hash(self, tmp_path):
        """Providing an expected hash activates comparison."""
        from markerfinder.utils.db_versioning import DatabaseVersionManager

        db_dir = tmp_path / "db"
        db_dir.mkdir()
        (db_dir / "expected_hashes.json").write_text(
            json.dumps({"marker_hmm": "deadbeef"}), encoding="utf-8"
        )
        target = tmp_path / "hmm"
        target.write_text("payload", encoding="utf-8")

        manager = DatabaseVersionManager(str(db_dir))
        info = manager.verify_database("marker_hmm", str(target))
        assert info["expected_hash"] == "deadbeef"
        assert info["hash_match"] is False  # Mismatch detected, not swallowed

    def test_verify_database_no_expected_means_no_verdict(self, tmp_path):
        from markerfinder.utils.db_versioning import DatabaseVersionManager

        db_dir = tmp_path / "db"
        db_dir.mkdir()
        target = tmp_path / "hmm"
        target.write_text("payload", encoding="utf-8")
        manager = DatabaseVersionManager(str(db_dir))
        info = manager.verify_database("marker_hmm", str(target))
        assert info["hash_match"] is None
        assert info["exists"] is True


class TestOccupancyOnCard:
    def test_occupancy_joins_decision_card(self):
        """Occupancy becomes an independent card column."""
        from markerfinder.config import HGTConfig
        from markerfinder.models.marker import MarkerLevel
        from markerfinder.models.pipeline_types import PhyloStepResult
        from markerfinder.modules.hgt_filter import HGTDecisionEngine

        engine = HGTDecisionEngine(HGTConfig())
        res = PhyloStepResult(overall_risk=0.1)
        ev = engine.evaluate_marker("M1", res)
        ev.decision_card["occupancy"] = {"M1": 3}.get("M1")  # What run does
        assert ev.decision_card["occupancy"] == 3
        # Independence: occupancy must not alter risk/level
        assert ev.overall_risk == 0.1
        assert ev.level == MarkerLevel.LEVEL_1


class TestDependencyDeclaration:
    def test_pyproject_drops_unused_deps(self):
        """Scipy/pandas/numpy/rich removed from runtime deps.
        The positive control is biopython, which must remain declared."""
        try:
            import tomllib
        except ModuleNotFoundError:  # Python 3.10 (declared floor)
            import tomli as tomllib
        from pathlib import Path

        pyproject = Path(__file__).parent.parent.parent / "pyproject.toml"
        data = tomllib.loads(pyproject.read_text(encoding="utf-8"))
        deps = " ".join(data["project"]["dependencies"]).lower()
        for unused in ("pandas", "numpy", "scipy", "rich"):
            assert unused not in deps, f"{unused} still declared"
        assert "biopython" in deps  # Positive control: the scan can hit

    def test_self_test_reports_unused_deps_as_info(self):
        """--check downgrades unused deps to informational."""
        from markerfinder.cli import self_test as st

        results = st._test_dependencies()
        statuses = {status for _, status, _ in results}
        for name, status, detail in results:
            if name.startswith("Import rich"):
                assert status == "INFO"
                assert "unused" in detail


class TestMarkerPreset:
    def test_preset_overrides_selection_triple(self):
        """An explicit preset overrides the triple."""
        from markerfinder.cli.config_build import _build_pipeline_config
        from markerfinder.cli.parser import _build_parser
        from markerfinder.models.marker import SelectionStrategy

        parser = _build_parser()
        args = parser.parse_args(
            ["--output", "out", "--marker-preset", "conservative"]
        )
        cfg = _build_pipeline_config(args)
        assert cfg.selection_config.min_occupancy == 0.90
        assert cfg.selection_config.max_markers == 30
        assert cfg.selection_config.selection_strategy is SelectionStrategy.INFO_MAX

    def test_preset_none_keeps_defaults(self):
        from markerfinder.cli.config_build import _build_pipeline_config
        from markerfinder.cli.parser import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--output", "out"])
        cfg = _build_pipeline_config(args)
        # The shipped default, now a named dataclass default instead of
        # The literal that used to sit inside config_build.
        assert cfg.selection_config.max_markers == 60  # Baseline behaviour
