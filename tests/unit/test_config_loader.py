"""Unit tests for configuration file loading and merging."""

import json
from pathlib import Path

import pytest

from markerfinder.config_loader import (
    _load_json,
    _load_toml,
    _load_yaml,
    load_config_file,
    merge_config_with_args,
)
from markerfinder.exceptions import ConfigError


class TestLoadYaml:
    def test_load_simple_yaml(self, tmp_path):
        path = tmp_path / "cfg.yaml"
        path.write_text("threads: 4\noutput: /tmp/out\n")
        data = _load_yaml(str(path))
        assert data["threads"] == 4
        assert data["output"] == "/tmp/out"

    def test_load_empty_yaml(self, tmp_path):
        path = tmp_path / "empty.yaml"
        path.write_text("")
        assert _load_yaml(str(path)) == {}

    def test_load_non_mapping_yaml(self, tmp_path):
        path = tmp_path / "bad.yaml"
        path.write_text("- a\n- b\n")
        with pytest.raises(ConfigError):
            _load_yaml(str(path))


class TestLoadToml:
    def test_load_simple_toml(self, tmp_path):
        path = tmp_path / "cfg.toml"
        path.write_text('threads = 4\noutput = "/tmp/out"\n')
        data = _load_toml(str(path))
        assert data["threads"] == 4
        assert data["output"] == "/tmp/out"


class TestLoadJson:
    def test_load_simple_json(self, tmp_path):
        path = tmp_path / "cfg.json"
        path.write_text(json.dumps({"threads": 4, "output": "/tmp/out"}))
        data = _load_json(str(path))
        assert data["threads"] == 4

    def test_load_non_mapping_json(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(["a", "b"]))
        with pytest.raises(ConfigError):
            _load_json(str(path))


class TestLoadConfigFile:
    def test_load_yaml_by_extension(self, tmp_path):
        path = tmp_path / "cfg.yaml"
        path.write_text("threads: 4\n")
        assert load_config_file(str(path))["threads"] == 4

    def test_load_unknown_extension_tries_all_formats(self, tmp_path):
        path = tmp_path / "cfg"
        path.write_text("threads: 4\n")
        assert load_config_file(str(path))["threads"] == 4

    def test_load_unparseable_file_raises(self, tmp_path):
        path = tmp_path / "cfg"
        path.write_text("not valid yaml or toml or json")
        with pytest.raises(Exception):
            load_config_file(str(path))

    def test_load_missing_file_raises(self, tmp_path):
        with pytest.raises(ConfigError):
            load_config_file(str(tmp_path / "missing.yaml"))


class TestMergeConfigWithArgs:
    def test_cli_overrides_config(self):
        defaults = {"threads": 1, "output": "./out"}
        config = {"threads": 4}
        args = {"threads": 8, "output": "./cli"}
        merged = merge_config_with_args(config, args, defaults)
        assert merged["threads"] == 8
        assert merged["output"] == "./cli"

    def test_config_overrides_defaults(self):
        defaults = {"threads": 1, "output": "./out"}
        config = {"threads": 4}
        args = {}
        merged = merge_config_with_args(config, args, defaults)
        assert merged["threads"] == 4
        assert merged["output"] == "./out"

    def test_arg_equal_to_default_is_ignored(self):
        defaults = {"threads": 1}
        config = {"threads": 4}
        args = {"threads": 1}
        merged = merge_config_with_args(config, args, defaults)
        assert merged["threads"] == 4

    def test_hyphen_keys_normalized(self):
        defaults = {"output_dir": "./out"}
        config = {"output-dir": "/tmp/out"}
        merged = merge_config_with_args(config, {}, defaults)
        assert merged["output_dir"] == "/tmp/out"


class TestRecordedSnapshotReplay:
    """``--config <previous run>/Phase5_metadata/run_config.json`` must work.

    The README's reproducibility promise is that the file the software wrote is
    enough to reproduce the run. The snapshot nests the per-module dataclasses
    and renames a few fields, so reading it verbatim left the values that decide
    the run — marker directory, taxonomy table, requested marker budget — behind
    names the command line does not have, and the replay silently ran on
    defaults while reporting the recorded numbers back.
    """

    SNAPSHOT = {
        "markerfinder_version": "0.1.1",
        "timestamp": "2026-09-30T00:00:00",
        "run_duration_seconds": 12.5,
        "database_versions": {"gtdb_markers": {"hash": "ab"}},
        "parameters": {
            "input_dir": "/run/in",
            "output_dir": "/run/out",
            "tmp_dir": "/run/tmp/scratch",
            "tmp_dir_auto": True,
            "cpus": 2,
            "db_dir": "/install/db",
            "selection_config": {
                "marker_mode": "gtdb_tk",
                "gtdb_markers_dir": "/run/markers",
                "max_markers": 150,
                "user_max_markers": 4,
                "min_occupancy": 0.75,
                "user_min_occupancy": None,
            },
            "taxonomy_config": {"taxonomy_table": "/run/tax.tsv"},
            "hgt_config": {
                "level_thresholds": {"level1_max": 0.2, "level2_max": 0.7},
                "monophyly_rank": "order",
                "adaptive_far_thresholds": True,
            },
        },
    }

    @pytest.fixture
    def snapshot_path(self, tmp_path):
        meta = tmp_path / "out" / "Phase5_metadata"
        meta.mkdir(parents=True)
        path = meta / "run_config.json"
        path.write_text(json.dumps(self.SNAPSHOT), encoding="utf-8")
        return path

    def test_nested_blocks_are_flattened_onto_command_line_names(
            self, snapshot_path):
        data = load_config_file(str(snapshot_path))
        assert data["gtdb_markers_dir"] == "/run/markers"
        assert data["taxonomy_table"] == "/run/tax.tsv"
        assert data["monophyly_rank"] == "order"
        # Renamed fields come back under the names the parser owns.
        assert data["input"] == "/run/in"
        assert data["output"] == "/run/out"
        assert data["threads"] == 2
        assert data["hgt_threshold_l1l2"] == 0.2
        assert data["hgt_threshold_l2l3"] == 0.7
        assert data["hgt_adaptive_thresholds"] is True

    def test_the_requested_budget_wins_over_the_adapted_one(self, snapshot_path):
        data = load_config_file(str(snapshot_path))
        assert data["max_markers"] == 4, (
            "a snapshot recorded with --max-markers 4 replayed as the adapted "
            "budget of 150"
        )
        # No user intent recorded for the occupancy floor: the value in the
        # Snapshot is what Phase 0 derived for that dataset, and replaying it as
        # A request would pin an adaptive default. It is dropped instead.
        assert "min_occupancy" not in data, data.get("min_occupancy")
        assert "user_max_markers" not in data
        assert "user_min_occupancy" not in data

    def test_per_run_scratch_and_recording_envelope_are_not_replayed(
            self, snapshot_path):
        data = load_config_file(str(snapshot_path))
        # The temporary directory of the recorded run is not this run's scratch.
        assert "tmp_dir" not in data
        for envelope in ("timestamp", "run_duration_seconds",
                         "database_versions", "markerfinder_version",
                         "tmp_dir_auto", "parameters"):
            assert envelope not in data, envelope

    def test_a_hand_written_config_is_left_alone(self, tmp_path):
        """Only snapshots are flattened: a user's own nesting stays as written."""
        path = tmp_path / "cfg.json"
        path.write_text(json.dumps({"selection_config": {"max_markers": 9}}),
                        encoding="utf-8")
        data = load_config_file(str(path))
        assert data == {"selection_config": {"max_markers": 9}}
