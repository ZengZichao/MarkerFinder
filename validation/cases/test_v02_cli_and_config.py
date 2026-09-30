"""V-02 — command line, configuration files and the parameter priority.

Covers ``--config`` (YAML / TOML / JSON), the documented CLI > config >
default precedence, the unknown-key warning, ``--threads`` bounds, ``-v``
verbosity, ``--log-file`` persistence, ``--output`` creation rules
(``--force`` / ``--no-clobber``) and the legacy no-subcommand invocation.
"""

from __future__ import annotations

import json

import pytest

from markerfinder.cli import constants

# The defaults `mf` puts on the command line, which a config-file case must
# Suppress so that the value can only have come from the file.
CONFIG_OWNED = ("--marker-mode", "--gtdb-markers-dir", "--taxonomy-table",
                "--max-markers", "--monophyly-rank", "--gene-tree-builder",
                "--coalescent-mode", "--skip-checkm")


@pytest.mark.capability("config", "input", "threads",
                        "subcommand:aggregate-run", "exit:EXIT_SUCCESS")
def test_yaml_config_drives_a_whole_run(mf, data_dir, record_metric):
    """A config file must be able to run the pipeline on its own.

    Only ``-i`` and ``-o`` stay on the command line (the input directory is a
    per-worker runtime path, and the output location is what lets the case find
    the products). Everything else in the assertions below can only have come
    from the file, because the matching command-line defaults are omitted.
    """
    cfg = data_dir / "configs" / "quad4.yaml"
    run = mf(omit=CONFIG_OWNED, extra=["--config", str(cfg)])
    run.assert_ok("--config YAML run")
    assert run.recorded("selection_config", "marker_mode") == "gtdb_tk"
    assert run.recorded("selection_config", "max_markers") == 4, (
        run.recorded("selection_config", "max_markers"))
    assert run.recorded("", "cpus") == 2, run.recorded("", "cpus")
    assert "quad4_core" in str(run.recorded("selection_config",
                                            "gtdb_markers_dir")), (
        run.recorded("selection_config", "gtdb_markers_dir"))
    assert run.recorded("phylo_config", "coalescent_mode") == "off"
    assert run.recorded("mag_config", "skip_checkm") is True
    record_metric("v02_config", "yaml_run_rc", run.rc)


@pytest.mark.capability("config")
def test_toml_config_drives_a_whole_run(mf, data_dir):
    cfg = data_dir / "configs" / "quad4.toml"
    run = mf(omit=CONFIG_OWNED, extra=["--config", str(cfg)])
    run.assert_ok("--config TOML run")
    assert run.recorded("selection_config", "max_markers") == 4
    assert run.recorded("taxonomy_config", "taxonomy_table").endswith(
        "taxonomy_quad4.tsv"), run.recorded("taxonomy_config", "taxonomy_table")


@pytest.mark.capability("config")
def test_json_config_with_a_parameters_block_is_flattened(mf, data_dir):
    """``--config`` accepts the JSON shape a run records (a top-level
    ``parameters`` mapping), because replaying a recorded run is the documented
    reproducibility route."""
    cfg = data_dir / "configs" / "quad4_flat_params.json"
    payload = json.loads(cfg.read_text(encoding="utf-8"))
    intended = payload["parameters"]
    assert "marker_mode" in intended, intended

    run = mf(omit=CONFIG_OWNED, extra=["--config", str(cfg)])
    run.assert_ok("JSON --config run")
    assert run.recorded("selection_config", "marker_mode") == "gtdb_tk"
    assert run.recorded("selection_config", "max_markers") == intended["max_markers"]


@pytest.mark.capability("config", "workflow:determinism")
def test_a_recorded_run_config_can_be_replayed(mf, record_metric):
    """Feed run A's own ``run_config.json`` back through ``--config``.

    This is the strongest form of the reproducibility claim: the file the
    software wrote must be enough to reproduce the same configuration, not
    merely a log of it.
    """
    first = mf(extra=["--force"])
    first.assert_ok("baseline run")
    recorded = first.product("Phase5_metadata/run_config.json")
    replay = mf(omit=CONFIG_OWNED, extra=["--config", str(recorded), "--force"])
    replay.assert_ok("replay from a recorded run_config.json")

    a = first.recorded_parameters()
    b = replay.recorded_parameters()
    # The two runs may differ only where the case itself moved them: output
    # Directory, temp directory, input path. Everything else must agree.
    volatile = {"input_dir", "output_dir", "tmp_dir", "keep_tmp", "timestamp",
                "run_duration_seconds"}

    def flatten(node, prefix=""):
        if isinstance(node, dict):
            out = {}
            for k, v in node.items():
                out.update(flatten(v, f"{prefix}{k}."))
            return out
        if isinstance(node, list):
            return {prefix.rstrip("."): sorted(map(str, node))}
        return {prefix.rstrip("."): node}

    fa, fb = flatten(a), flatten(b)
    diffs = {}
    for key in set(fa) | set(fb):
        leaf = key.rsplit(".", 1)[-1]
        if leaf in volatile:
            continue
        if fa.get(key) != fb.get(key):
            diffs[key] = (fa.get(key), fb.get(key))
    record_metric("v02_replay", "differing_parameters", sorted(diffs))
    assert not diffs, (
        f"replay did not reproduce the recorded configuration: {diffs}"
    )


@pytest.mark.capability("config", "workflow:failure-loudness")
def test_unknown_config_key_is_warned_about_not_swallowed(mf, case_workdir):
    """A key the software does not have must never look like it took effect."""
    cfg = case_workdir / "bogus.yaml"
    cfg.write_text("adaptive_far_thresholds: true\nhgt_mode: nonsense\n",
                   encoding="utf-8")
    run = mf(extra=["--config", str(cfg)])
    run.assert_ok("run with unknown config keys")
    assert "unrecognized keys" in run.text.lower(), (
        f"the unknown keys were silently ignored:\n{run.tail()}")
    assert "adaptive_far_thresholds" in run.text, run.text[:400]


@pytest.mark.capability("config", "mode")
def test_cli_overrides_the_config_file(mf, case_workdir):
    """Same key from two sources: the command line wins; a key only the file
    supplies still reaches the pipeline."""
    cfg = case_workdir / "mode.yaml"
    cfg.write_text("mode: conservative\nmax_markers: 6\n", encoding="utf-8")
    # --max-markers is omitted from the command line so the file's value is the
    # Only possible source for it; --mode is passed explicitly against the file.
    # The budget is 6 rather than 3 because a budget that leaves no marker
    # Passing the HGT screen is refused (exit 3), and then this case could not
    # Read the record it is about.
    run = mf(omit=("--max-markers",),
             extra=["--config", str(cfg), "--mode", "expanded"])
    run.assert_ok()
    assert run.recorded("", "mode") == "expanded", run.recorded("", "mode")
    assert run.recorded("selection_config", "max_markers") == 6, (
        "a key that only the config file supplied was not applied: "
        f"{run.recorded('selection_config', 'max_markers')}"
    )


@pytest.mark.capability("threads", "exit:EXIT_ARG_ERROR")
def test_thread_count_below_one_is_an_argument_error(run_markerfinder,
                                                    genome_input, tmp_path):
    result = run_markerfinder(["-i", str(genome_input("quad4")),
                               "-o", str(tmp_path / "o"), "-t", "0"])
    assert result.rc == constants.EXIT_ARG_ERROR, result.tail()
    assert "below minimum" in result.text.lower(), result.tail()


@pytest.mark.capability("verbose")
def test_verbose_flag_raises_the_output_verbosity(mf, record_metric):
    """``-v`` must actually add log lines, not just be accepted."""
    quiet = mf(extra=["--force"])
    quiet.assert_ok("quiet baseline run")
    verbose = mf(extra=["-v", "--force"])
    verbose.assert_ok("verbose run")
    n_quiet = len(quiet.text.splitlines())
    n_verbose = len(verbose.text.splitlines())
    assert n_verbose > n_quiet, (
        f"-v produced {n_verbose} lines against {n_quiet} without it"
    )
    record_metric("v02_verbosity", "default_lines", n_quiet)
    record_metric("v02_verbosity", "verbose_lines", n_verbose)


@pytest.mark.capability("log_file")
def test_log_file_is_written_in_utf8_and_holds_the_run(mf, case_workdir):
    """The log file must carry the run's own progress lines.

    At the default verbosity the logger only emits warnings, so a file written
    by a quiet run can legitimately be short; ``-v`` is what makes the phase
    narration land in it, and that is the combination a user is told to use.
    """
    log = case_workdir / "run.log"
    run = mf(extra=["--log-file", str(log), "-v", "--force"])
    run.assert_ok()
    assert log.exists(), f"--log-file {log} was never created"
    text = log.read_text(encoding="utf-8")
    assert len(text) > 100, f"the log file is a stub ({len(text)} chars)"
    assert "Dependency Check" in text or "Phase" in text, text[:400]


@pytest.mark.capability("output", "force", "no_clobber")
def test_existing_output_requires_force_or_no_clobber(mf, case_output,
                                                      record_metric):
    first = mf(out_dir=case_output, extra=["--force"])
    first.assert_ok("first run")
    second = mf(out_dir=case_output, omit=("--force",))
    assert second.rc != 0, (
        "a run into an existing output directory succeeded without --force or "
        f"--no-clobber:\n{second.tail()}")
    assert "exists" in second.text.lower() or "force" in second.text.lower(), (
        f"the refusal does not name its reason:\n{second.tail()}")
    third = mf(out_dir=case_output, omit=("--force",), extra=["--no-clobber"])
    assert third.rc == 0, f"--no-clobber should let the run finish: {third.tail()}"


@pytest.mark.capability("output", "subcommand:aggregate-run")
def test_output_directory_is_created_when_absent(mf, case_output):
    target = case_output / "nested" / "out"
    run = mf(out_dir=target, extra=["--force"])
    run.assert_ok()
    assert target.is_dir(), f"{target} was not created"
    assert (target / "Phase5_reports").is_dir(), sorted(target.iterdir())
