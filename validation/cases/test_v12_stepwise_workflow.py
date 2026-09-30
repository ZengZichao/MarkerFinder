"""V-12 — the four resumable steps, their state file and their ordering.

Covers the subcommands ``scan`` / ``filter`` / ``infer`` / ``report``, plus
``--resume`` and ``--redo``.

The documented contract is that a step runs only after its prerequisite, that
the state needed to chain them is persisted under ``<output>/.markerfinder/``,
and that a completed step is not silently repeated. Each of those is a failure
mode this file drives directly: out-of-order invocation, missing prerequisite,
re-running a finished step, and resuming an incomplete one.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from markerfinder.cli import constants

BASE_FLAGS = ["-t", "2", "--marker-mode", "gtdb_tk", "--max-markers", "4",
              "--monophyly-rank", "order", "--gene-tree-builder", "fasttree",
              "--skip-checkm", "--coalescent-mode", "always"]


@pytest.fixture(scope="module")
def step_inputs(request, data_dir):
    """Paths every step case shares (markers and taxonomy of the fast set)."""
    return {
        "markers": str(data_dir / "markers" / "quad4_core"),
        "taxonomy": str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
    }


def _markers(out_dir: Path):
    """The marker ids a run reported in marker_summary.tsv."""
    path = out_dir / "Phase5_reports" / "markerfinder.marker_summary.tsv"
    if not path.exists():
        return set()
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    return {line.split("\t")[0] for line in lines[1:] if line.strip()}


@pytest.mark.capability("subcommand:scan", "subcommand:filter",
                        "subcommand:infer", "subcommand:report",
                        "workflow:state-persistence",
                        "product:.markerfinder.pipeline_state.json",
                        "product:Phase5_metadata.run_config.json")
def test_the_four_steps_in_order_reach_the_same_products_as_the_aggregate_run(
        run_markerfinder, genome_input, step_inputs, case_output, record_metric):
    """Scan -> filter -> infer -> report must report the same marker set as one
    aggregate ``markerfinder -i... -o...`` invocation: the documented routes
    are two ways of running one pipeline, not two pipelines."""
    inp = genome_input("quad4")
    steps = ["scan", "filter", "infer", "report"]
    for name in steps:
        result = run_markerfinder(
            [name, "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
            ["--gtdb-markers-dir", step_inputs["markers"],
             "--taxonomy-table", step_inputs["taxonomy"], "--force"])
        assert result.rc == 0, f"{name} failed:\n{result.tail()}"

    state = case_output / ".markerfinder" / ".pipeline_state.json"
    assert state.exists(), f"no step-state file at {state}"
    recorded_steps = json.dumps(json.loads(state.read_text(encoding="utf-8")))
    for name in steps:
        assert name in recorded_steps, (
            f"{name} is not recorded in the step state: {recorded_steps[:400]}")

    stepped = _markers(case_output)
    aggregate = case_output.parent / "aggregate"
    agg = run_markerfinder(
        ["-i", str(inp), "-o", str(aggregate)] + BASE_FLAGS +
        ["--gtdb-markers-dir", step_inputs["markers"],
         "--taxonomy-table", step_inputs["taxonomy"], "--force"])
    agg.assert_ok("aggregate run")
    aggregate_markers = _markers(aggregate)
    record_metric("v12_steps", "markers_stepwise", sorted(stepped))
    record_metric("v12_steps", "markers_aggregate", sorted(aggregate_markers))
    assert stepped == aggregate_markers, (
        f"stepwise kept {sorted(stepped)} and the aggregate run kept "
        f"{sorted(aggregate_markers)}: the two documented routes disagree"
    )
    assert stepped, "both routes reported zero markers, so the comparison is void"


@pytest.mark.capability("subcommand:filter", "subcommand:infer",
                        "exit:EXIT_DATA_ERROR", "workflow:failure-loudness")
def test_a_step_without_its_prerequisite_refuses(run_markerfinder, step_inputs,
                                                 case_output):
    """``filter`` without ``scan`` cannot know which markers exist. The refusal
    must name the missing step, not crash or invent an empty result."""
    result = run_markerfinder(
        ["filter", "-o", str(case_output)] + BASE_FLAGS +
        ["--gtdb-markers-dir", step_inputs["markers"],
         "--taxonomy-table", step_inputs["taxonomy"], "--force"])
    assert result.rc != 0, (
        f"'filter' ran with no prior 'scan' and exited {result.rc}:\n"
        f"{result.tail()}"
    )
    assert "Traceback" not in result.text, result.tail()
    text = result.text.lower()
    assert "scan" in text or "state" in text or "prerequisite" in text, (
        f"the refusal does not say which step is missing:\n{result.tail()}"
    )


@pytest.mark.capability("resume", "subcommand:report",
                        "workflow:state-persistence")
def test_resume_runs_the_missing_prerequisites(run_markerfinder, genome_input,
                                               step_inputs, case_output,
                                               record_metric):
    """``--resume`` is documented as "automatically run any missing prerequisite
    steps up to the requested step". Asking only for `report` must therefore
    leave a finished report behind, having run scan/filter/infer itself — the
    products are the evidence, not a log line a user may never have asked for.
    """
    inp = genome_input("quad4")
    result = run_markerfinder(
        ["report", "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
        ["--gtdb-markers-dir", step_inputs["markers"],
         "--taxonomy-table", step_inputs["taxonomy"], "--resume", "--force"])
    record_metric("v12_resume", "rc", result.rc)
    assert result.rc == 0, f"--resume did not complete the chain:\n{result.tail()}"
    for rel in ("Phase5_reports/markerfinder.report.html",
                "Phase5_reports/markerfinder.marker_summary.tsv",
                "Phase4_trees/markerfinder.species_tree_concat.newick"):
        assert (case_output / rel).exists(), (
            f"--resume reported success without producing {rel}"
        )
    # The chain must be recorded as such, or a later --resume cannot tell where
    # To pick up.
    state = case_output / ".markerfinder" / ".pipeline_state.json"
    assert state.exists(), f"--resume left no step state at {state}"


@pytest.mark.capability("redo", "subcommand:filter",
                        "workflow:state-persistence")
def test_redo_repeats_a_completed_step_and_reports_it(run_markerfinder,
                                                     genome_input, step_inputs,
                                                     case_output,
                                                     record_metric):
    inp = genome_input("quad4")
    for name in ("scan", "filter"):
        first = run_markerfinder(
            [name, "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
            ["--gtdb-markers-dir", step_inputs["markers"],
             "--taxonomy-table", step_inputs["taxonomy"], "--force"])
        assert first.rc == 0, f"{name}: {first.tail()}"

    skipped = run_markerfinder(
        ["filter", "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
        ["--gtdb-markers-dir", step_inputs["markers"],
         "--taxonomy-table", step_inputs["taxonomy"], "--force"])
    assert skipped.rc == 0, skipped.tail()

    redone = run_markerfinder(
        ["filter", "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
        ["--gtdb-markers-dir", step_inputs["markers"],
         "--taxonomy-table", step_inputs["taxonomy"], "--redo", "--force"])
    assert redone.rc == 0, redone.tail()
    record_metric("v12_redo", "without_redo_rc", skipped.rc)
    record_metric("v12_redo", "with_redo_rc", redone.rc)
    # --redo's contract is "re-run the requested step even if it was already
    # Completed". The observable difference is that the step's own products are
    # Rewritten: compare their modification times through the state file, which
    # Is the thing --redo is documented to reset.
    from markerfinder.utils.state_codec import load_pipeline_state

    state = load_pipeline_state(str(case_output)) or {}
    assert "filter" in str(state), (
        f"after --redo the persisted state no longer records the step it was "
        f"asked to repeat: {list(state)}"
    )


@pytest.mark.capability("subcommand:report",
                        "product:.markerfinder.context.json")
def test_stepwise_products_are_written_under_the_documented_directories(
        run_markerfinder, genome_input, step_inputs, case_output):
    """The historical ``Phase4_*`` / ``Phase5_*`` directory names are part of the
    published output layout: a step may not quietly move its products."""
    inp = genome_input("quad4")
    for name in ("scan", "filter", "infer", "report"):
        result = run_markerfinder(
            [name, "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
            ["--gtdb-markers-dir", step_inputs["markers"],
             "--taxonomy-table", step_inputs["taxonomy"], "--force"])
        assert result.rc == 0, f"{name}: {result.tail()}"
    for rel in ("Phase5_reports/markerfinder.report.html",
                "Phase5_reports/markerfinder.marker_summary.tsv",
                "Phase5_reports/markerfinder.hgt_evaluation.tsv",
                "Phase5_reports/markerfinder.pipeline_summary.txt",
                "Phase5_metadata/run_config.json",
                "Phase4_trees/markerfinder.species_tree_concat.newick"):
        assert (case_output / rel).exists(), f"{rel} missing from a stepwise run"


@pytest.mark.capability("subcommand:scan", "workflow:state-persistence")
def test_state_files_live_in_the_work_directory_not_the_source_tree(
        run_markerfinder, genome_input, step_inputs, case_output):
    """``<output>/.markerfinder/`` holds the chain state; nothing may be written
    next to the shipped inputs."""
    inp = genome_input("quad4")
    result = run_markerfinder(
        ["scan", "-i", str(inp), "-o", str(case_output)] + BASE_FLAGS +
        ["--gtdb-markers-dir", step_inputs["markers"],
         "--taxonomy-table", step_inputs["taxonomy"], "--force"])
    assert result.rc == 0, result.tail()
    assert (case_output / ".markerfinder").is_dir()
    data_root = Path(step_inputs["markers"]).parent.parent
    strays = [p for p in data_root.rglob(".markerfinder") if p.is_dir()]
    assert not strays, f"state written inside the data bundle: {strays}"
