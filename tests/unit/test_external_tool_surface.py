"""The external-tool pre-check table must match what the code invokes.

The finding this gate retires: ``EXTERNAL_TOOLS`` listed six binaries while the
package actually shells out to more. ``mag_optimization._run_checkm`` runs
``checkm lineage_wf`` (mag_optimization.py:121-133) and, when CheckM is missing,
degrades to default quality estimates -- quietly, from the user's point of view,
because the only place that could have warned them (the startup pre-check and
``--check``) had never heard of CheckM. A run then reported completeness numbers
that were estimates, labelled honestly in the report (``quality_source``) but
disclosed nowhere beforehand.

So the rule here is not "CheckM is in the list" (that would only pin this one
instance). The rule is: **derive** the invoked-binary set from the package's own
AST and require the declared table to agree with it. A new ``subprocess`` call
site for a tool nobody registered will now fail this test instead of failing a
user's run.

Two binaries are deliberately NOT in the table and stay that way:

  * ``git`` -- used to be stamped into the version by ``_version.py``; the
    version is a constant now, so nothing in the package may invoke it (the
    control test below asserts that it is gone rather than that it is found);
  * ``diamond`` -- belongs to the Phase 1.5 ortholog-resolution path, which the
    pipeline does not execute (``phases.PHASES["1.5"]``; disclosed as INFO in
    ``--check``). Telling users to install a tool that cannot change their
    results would be a new defect, not a fix.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PACKAGE = REPO / "markerfinder"

SUBPROCESS_FUNCS = {"run", "Popen", "call", "check_call", "check_output"}

# Binaries that are invoked by the package but are not pipeline capabilities.
NOT_A_PIPELINE_TOOL = {"git", "diamond"}

EXTERNAL_TOOLS, TOOL_STEPS, PIPELINE_STEPS = (
    __import__(
        "markerfinder.utils.dependency_check", fromlist=["EXTERNAL_TOOLS"]
    ).EXTERNAL_TOOLS,
    __import__(
        "markerfinder.utils.dependency_check", fromlist=["TOOL_STEPS"]
    ).TOOL_STEPS,
    __import__(
        "markerfinder.utils.dependency_check", fromlist=["PIPELINE_STEPS"]
    ).PIPELINE_STEPS,
)
TOOL_ALIASES = __import__(
    "markerfinder.utils.dependency_check", fromlist=["_TOOL_ALIASES"]
)._TOOL_ALIASES

# Canonical phases reachable from each CLI step (``scan`` covers Phases 0-1.5).
STEP_PHASES = {
    "scan": {"0", "0.1", "0.2", "1", "1.5"},
    "filter": {"2"},
    "infer": {"3"},
    "report": {"4"},
}

# Which pipeline step runs each module that shells out. ``utils.gtdb_tk_markers``
# Builds gene trees with mafft/trimal/fasttree and is reached from
# Marker_selection.py:449, i.e. during ``scan``.
MODULE_STEP = {
    "modules/mag_optimization.py": ["scan"],
    "modules/marker_selection.py": ["scan"],
    "modules/hgt_filter.py": ["filter"],
    "modules/phylogenetic_inference.py": ["infer"],
    "modules/report_generator.py": ["report"],
    "modules/ortholog_resolver.py": [],  # Reserved interface: not wired in
    "utils/gtdb_tk_markers.py": ["scan"],
}

IS_ALIASED = re.compile(r"^[A-Za-z0-9_.+-]+$")


def _first_word(node: ast.AST):
    """The binary a subprocess argument names, if it is a literal command."""
    if isinstance(node, ast.List) and node.elts:
        head = node.elts[0]
        if isinstance(head, ast.Constant) and isinstance(head.value, str):
            return head.value
    return None


def invoked_binaries() -> dict[str, list[str]]:
    """Binary -> modules that invoke it, read off the package's own AST."""
    found: dict[str, set[str]] = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        # ``cmd = ["checkm",...]`` -- remember the local, then follow its use.
        assigned: dict[str, str] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.List):
                word = _first_word(node.value)
                if word and node.targets and isinstance(node.targets[0], ast.Name):
                    assigned[node.targets[0].id] = word
            if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)):
                continue
            if node.func.attr not in SUBPROCESS_FUNCS:
                continue
            base = node.func.value
            if not (isinstance(base, ast.Name) and base.id == "subprocess"):
                continue
            if not node.args:
                continue
            arg = node.args[0]
            word = _first_word(arg)
            if word is None and isinstance(arg, ast.Name):
                word = assigned.get(arg.id)
            if word and IS_ALIASED.match(word):
                found.setdefault(word, set()).add(str(path.relative_to(PACKAGE)))
    return {k: sorted(v) for k, v in found.items()}


def test_control_the_binary_discovery_actually_finds_something():
    """A scanner that finds nothing would make every check below pass vacuously."""
    binaries = invoked_binaries()
    assert "hmmsearch" in binaries, binaries
    assert "checkm" in binaries, "the F-4 call site must be discoverable"
    assert "git" not in binaries, (
        "the version is a constant in _version.py: nothing may shell out to git "
        "to stamp a checkout into it")
    assert "mafft" in binaries and "astral" in binaries


def _explained_names() -> set[str]:
    """Declared names plus every alias they stand in for (fasttree ~ FastTree)."""
    declared = {name for name, _purpose, _critical in EXTERNAL_TOOLS}
    return declared | {a for aliases in TOOL_ALIASES.values() for a in aliases}


def test_every_invoked_binary_is_declared_or_explicitly_excused():
    binaries = set(invoked_binaries())
    unexplained = sorted(binaries - _explained_names() - NOT_A_PIPELINE_TOOL)
    assert not unexplained, (
        f"subprocess call sites for tools missing from EXTERNAL_TOOLS: {unexplained}; "
        f"declared: {sorted(_explained_names())}; excused: {sorted(NOT_A_PIPELINE_TOOL)}"
    )


def test_the_table_names_no_tool_the_package_never_calls():
    binaries = set(invoked_binaries())
    unexplained = sorted(
        name for name in {n for n, _p, _c in EXTERNAL_TOOLS} - binaries
        if all(alias not in binaries for alias in TOOL_ALIASES.get(name, [name]))
    )
    assert not unexplained, (
        f"EXTERNAL_TOOLS advertises tools with no call site: {unexplained}; "
        f"invoked: {sorted(binaries)}"
    )


def test_control_alias_names_count_as_one_tool():
    """``FastTree`` on macOS and ``fasttree`` on Linux must not read as a
    missing or phantom tool -- that would make both gates above lie."""
    assert "FastTree" in invoked_binaries()
    assert "FastTree" in _explained_names()
    assert "FastTree" not in (set(invoked_binaries()) - _explained_names())


def test_tool_purpose_phase_citations_match_the_real_steps():
    """The purpose string is dependency documentation: a phase in it is a claim.

    This is what caught ``iqtree3`` being described as a Phase 3 tool when its
    only call site (gtdb_tk_markers.py:287) runs during ``scan``.
    """
    offenders = {}
    for name, purpose, _critical in EXTERNAL_TOOLS:
        cited = set(re.findall(r"Phase ([0-9.]+)", purpose))
        cited |= {p for chunk in re.findall(r"Phase ([0-9./]+)", purpose)
                  for p in re.split(r"[/.]", chunk) if p}
        allowed = set().union(*(STEP_PHASES[s] for s in TOOL_STEPS[name]))
        if cited - allowed:
            offenders[name] = (sorted(cited - allowed), sorted(TOOL_STEPS[name]))
    assert not offenders, (
        f"purpose text cites phases the tool is not invoked in (cited, steps): {offenders}"
    )


def test_checkm_is_part_of_the_disclosed_dependency_surface():
    """The regression itself: Phase 0 can silently fall back without CheckM."""
    entry = {name: (purpose, critical) for name, purpose, critical in EXTERNAL_TOOLS}
    assert "checkm" in entry, "CheckM is invoked by mag_optimization and must be listed"
    purpose, critical = entry["checkm"]
    assert critical is False, "CheckM absence degrades, it must not hard-fail the run"
    assert "Phase 0" in purpose and "default" in purpose.lower(), purpose


def test_diamond_stays_out_of_the_startup_pre_check():
    declared = {name for name, _purpose, _critical in EXTERNAL_TOOLS}
    assert "diamond" not in declared, (
        "diamond belongs to the unwired Phase 1.5 path; if Phase 1.5 is ever "
        "wired in, move this assertion and disclose it in --check"
    )


def test_tool_step_mapping_covers_the_table_and_matches_the_call_sites():
    declared = {name for name, _purpose, _critical in EXTERNAL_TOOLS}
    assert set(TOOL_STEPS) == declared, (
        f"TOOL_STEPS and EXTERNAL_TOOLS disagree: {set(TOOL_STEPS) ^ declared}"
    )
    binaries = invoked_binaries()
    offenders = {}
    for tool, steps in TOOL_STEPS.items():
        real_steps = sorted(
            {
                step
                for module in binaries.get(tool, [])
                for step in MODULE_STEP.get(module, ["<unknown module>"])
            }
        )
        if set(steps) != set(real_steps):
            offenders[tool] = (sorted(steps), real_steps)
    assert not offenders, (
        f"TOOL_STEPS must name the steps that really invoke each binary "
        f"(declared, derived): {offenders}"
    )


def test_control_step_filtering_rejects_an_invented_step():
    from markerfinder.utils.dependency_check import check_external_tools

    # Without this gate the argument is accepted and ignored, so
    # ["repot"] silently meant "check everything" and reported a clean bill.
    with pytest.raises(ValueError):
        check_external_tools(required_phases=["repot"])


def test_step_filtering_narrows_the_checked_set():
    from markerfinder.utils.dependency_check import check_external_tools

    scan_only = check_external_tools(required_phases=["scan"])
    assert "hmmsearch" in scan_only and "checkm" in scan_only, scan_only
    assert "astral" not in scan_only, "ASTRAL is an infer-step tool"
    assert set(check_external_tools(required_phases=["report"])) == set()
    assert set(check_external_tools()) == {
        name for name, _purpose, _critical in EXTERNAL_TOOLS
    }


def test_check_discloses_every_declared_tool_without_scoring_it():
    from markerfinder.cli.self_test import _test_external_tools

    rows = _test_external_tools()
    declared = [name for name, _purpose, _critical in EXTERNAL_TOOLS]
    assert len(rows) == len(declared), rows
    for name, status, detail in rows:
        assert name.startswith("External tool: "), name
        assert status == "INFO", (
            f"{name} is scored; PATH is per-machine and would move the "
            "documented --check item count"
        )
        assert detail
    checkm_row = [r for r in rows if r[0].endswith("checkm (optional)")][0]
    assert "default_no_checkm" in checkm_row[2] or "on PATH" == checkm_row[2], checkm_row
