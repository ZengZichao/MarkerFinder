"""Phase numbers in logs and comments must resolve to the one table.

 created ``markerfinder/phases.py`` as the single source of truth and kept
the legacy directory names with an explicit alias table. asks for the
consequence: the log labels and comments must not point anywhere else. The
baseline symptom was comments calling the coalescent step "Phase 4" and the
reports "Phase 5" while the canonical numbers are 3 and 4 -- exactly the kind of
drift that makes a reader distrust every other number in the product.

Nothing enforced that relationship: no test named. The extractor below is
shared by both directions, and each negative lock carries a positive control so
"grep found nothing" cannot come from a broken pattern.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

from markerfinder.phases import HISTORICAL_ALIASES, PHASES

PACKAGE = Path(__file__).resolve().parents[2] / "markerfinder"

TAG = re.compile(r"\[Phase ([0-9.]+)\]")
LEGACY_DIR = re.compile(r"\bPhase(\d)_([A-Za-z_]+)")


def _python_files():
    return sorted(PACKAGE.rglob("*.py"))


def _undeclared_tags(text: str) -> list:
    declared = {phase.number for phase in PHASES}
    return sorted({t for t in TAG.findall(text) if t not in declared})


def test_control_the_tag_extractor_fires_on_a_planted_stray():
    assert _undeclared_tags("logger.info('[Phase 9] nope')") == ["9"]
    assert _undeclared_tags("logger.info('[Phase 3] ok')") == []


def test_no_log_line_uses_an_undeclared_phase_tag():
    offenders = {}
    for path in _python_files():
        stray = _undeclared_tags(path.read_text(encoding="utf-8"))
        if stray:
            offenders[str(path.relative_to(PACKAGE))] = stray
    assert not offenders, (
        "log labels must come from phases.PHASES; declared numbers are "
        f"{sorted(p.number for p in PHASES)}; offenders: {offenders}"
    )


def test_every_tag_matches_the_canonical_form_its_number_would_produce():
    """`[Phase 3 ]`, `[phase 3]`, `[Phase3]` all name phase 3 but bypass the table."""
    offenders = []
    for path in _python_files():
        for line_no, line in enumerate(
            path.read_text(encoding="utf-8").split("\n"), 1
        ):
            for match in re.finditer(r"\[\s*[Pp]hase\s*([0-9.]+)\s*\]", line):
                number = match.group(1)
                try:
                    canonical = f"[Phase {number}]"
                    from markerfinder.phases import canonical_log_tag

                    declared = canonical_log_tag(number)
                except KeyError:
                    offenders.append((path.name, line_no, number, "undeclared"))
                    continue
                if match.group(0) != declared or canonical != declared:
                    offenders.append(
                        (path.name, line_no, match.group(0), declared)
                    )
    assert not offenders, offenders


def test_legacy_phase_directory_names_stay_registered_as_aliases():
    """Any `PhaseN_something` still in the tree must be a registered alias, so a
    reader (and tooling) can map the artefact back to its canonical phase."""
    registered = {
        key for key in HISTORICAL_ALIASES
        if LEGACY_DIR.fullmatch(key)
    }
    found = set()
    for path in _python_files():
        for match in LEGACY_DIR.finditer(path.read_text(encoding="utf-8")):
            found.add(match.group(0))
    unknown = sorted(found - registered)
    assert not unknown, (
        f"unregistered legacy phase directories {unknown}; registered aliases: "
        f"{sorted(registered)}"
    )


def test_control_the_alias_scanner_is_not_vacuous():
    text = "written into Phase5_reports and Phase7_mystery"
    assert {m.group(0) for m in LEGACY_DIR.finditer(text)} == {
        "Phase5_reports",
        "Phase7_mystery",
    }


def test_phase_table_is_unique_and_ordered_by_declaration():
    numbers = [p.number for p in PHASES]
    assert len(set(numbers)) == len(numbers), numbers
    assert all(p.log_tag == f"[Phase {p.number}]" for p in PHASES), numbers
    assert numbers == ["0", "0.1", "0.2", "1", "1.5", "2", "3", "4"]


# --- CLI help text ('s user-visible form) ------------------------------
#
# The bracketed-log ratchet above cannot see ``--help`` prose, and an
# Undeclared-number scan cannot catch the actual drift this test found: the
# Subcommand help cited DECLARED numbers on the WRONG steps (filter "Phase 3",
# Infer "Phase 4", report "Phase 5" -- the historical DIRECTORY numbering),
# So the mapping itself is locked, per subcommand.

EXPECTED_HELP_PHASES = {
    "scan": {"0", "1", "1.5"},
    "filter": {"2"},
    "infer": {"3"},
    "report": {"4"},
}


def _help_phase_numbers(spec) -> set:
    from markerfinder.cli.parser import _SUBCOMMAND_SPECS

    body = " ".join(piece for piece in spec[1:] if isinstance(piece, str))
    numbers: set = set()
    # "Phase 0+1+1.5" is one compound citation of three phases; split it.
    for compound in re.findall(r"\bPhase ([0-9.+]+)", body):
        numbers.update(part for part in compound.split("+") if part)
    return numbers


def test_cli_help_names_exactly_the_canonical_phase_of_its_step():
    from markerfinder.cli.parser import _SUBCOMMAND_SPECS

    by_name = {spec[0]: spec for spec in _SUBCOMMAND_SPECS}
    assert set(by_name) == set(EXPECTED_HELP_PHASES), (
        "a subcommand was added without declaring its canonical phase here"
    )
    offenders = {
        name: sorted(_help_phase_numbers(spec))
        for name, spec in by_name.items()
        if _help_phase_numbers(spec) != EXPECTED_HELP_PHASES[name]
    }
    assert not offenders, (
        "subcommand --help must cite the canonical phase from phases.PHASES; "
        f"expected {EXPECTED_HELP_PHASES}, offenders: {offenders}"
    )


def test_control_the_help_mapping_lock_fires_on_historical_numbering():
    """Positive control: the wrong report help ("Phase 5: generates...")
    must fail the mapping check, not slip through as a declared number."""
    from markerfinder.cli.parser import _SUBCOMMAND_SPECS

    report_spec = next(s for s in _SUBCOMMAND_SPECS if s[0] == "report")
    drifted = ("report", report_spec[1], report_spec[2].replace("Phase 4:", "Phase 5:", 1))
    assert _help_phase_numbers(drifted) == {"4", "5"}
    assert _help_phase_numbers(drifted) != EXPECTED_HELP_PHASES["report"]


# --- step prose and emitted log labels --------------------------------
#
# Nothing above could see the drift that was actually found, because it lived
# In two forms the bracketed scanner cannot match:
#
# * ``# Step 2: filter (Phase 3)`` in ``pipeline.py`` -- contradicting the
# ``Phase 2`` docstring a few lines below it, in the same file;
# * ``logger.info("Phase 5: Generating interactive HTML report")`` in
# ``report_generator.py`` -- an un-bracketed label naming a phase that does
# Not exist, printed at runtime on every report run.
#
# 's other half is locked too: the table used to have zero runtime callers,
# Which is how it drifted from being "the single source of truth" to being a
# Parallel fact. Its labels must now flow through ``canonical_log_tag``.

STEP_PROSE = re.compile(r"#\s*Step\s+\d:\s+(\w+)\s+\(([^)]*)\)")
DOC_SUBCOMMAND = re.compile(r'"""Execute the ``(\w+)`` subcommand \(([^)]*)\)')
LITERAL_LABEL = re.compile(r"(\[Phase [0-9.]+\]|Phase [0-9.]+:)")


def _prose_phase_numbers(text: str) -> set:
    """Numbers cited as ``Phase 2``, ``Phase 0+1+1.5`` or ``Phase 0 + 1 + 1.5``."""
    numbers: set = set()
    pattern = re.compile(
        r"Phase\s+([0-9](?:\.[0-9])?(?:\s*\+\s*[0-9](?:\.[0-9])?)*)"
    )
    for match in pattern.finditer(text):
        numbers.update(p for p in re.split(r"\s*\+\s*", match.group(1)) if p)
    return numbers


def test_step_comments_name_the_canonical_phase_of_their_step():
    text = (PACKAGE / "pipeline.py").read_text(encoding="utf-8")
    found = STEP_PROSE.findall(text)
    assert len(found) >= 4, f"the step-comment scanner found nothing to check: {found}"
    offenders = {
        name: sorted(_prose_phase_numbers(cited))
        for name, cited in found
        if name in EXPECTED_HELP_PHASES
        and _prose_phase_numbers(cited) != EXPECTED_HELP_PHASES[name]
    }
    assert not offenders, (
        "# Step k comments must cite the canonical phase of that step "
        f"({EXPECTED_HELP_PHASES}); offenders: {offenders}"
    )


def test_subcommand_docstrings_name_the_canonical_phase_of_their_step():
    commands = PACKAGE / "cli" / "commands"
    checked = 0
    offenders = {}
    for path in sorted(commands.glob("*.py")):
        for name, cited in DOC_SUBCOMMAND.findall(path.read_text(encoding="utf-8")):
            checked += 1
            got = _prose_phase_numbers(cited)
            if got != EXPECTED_HELP_PHASES[name]:
                offenders[path.name] = (sorted(got), sorted(EXPECTED_HELP_PHASES[name]))
    assert checked == len(EXPECTED_HELP_PHASES), (
        f"one docstring per subcommand expected, matched {checked}"
    )
    assert not offenders, f"subcommand docstring phase numbers drifted: {offenders}"


def test_control_the_step_prose_scanner_fires_on_the_retired_numbering():
    """The wrong comment is the planted positive: it must be catchable."""
    assert _prose_phase_numbers("Step 2: filter (Phase 3)") == {"3"}
    assert _prose_phase_numbers("Step 1: scan (Phase 0 + 1 + 1.5)") == {"0", "1", "1.5"}
    assert _prose_phase_numbers("(canonical phase 4)") == set(), "lowercase is not a citation"


def _label_source(arg):
    """The hand-written phase label an emitted log line opens with, if any."""
    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
        text = arg.value
    elif isinstance(arg, ast.JoinedStr) and arg.values:
        first = arg.values[0]
        if not (isinstance(first, ast.Constant) and isinstance(first.value, str)):
            return None  # Opens with an expression: canonical_log_tag -- correct
        text = first.value
    else:
        return None
    match = LITERAL_LABEL.match(text.lstrip())
    return match.group(0) if match else None


def _log_call_label_args(tree: ast.AST):
    """First argument of every ``logger.<level>`` call in a parsed module."""
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        fn = node.func
        if (
            isinstance(fn, ast.Attribute)
            and isinstance(fn.value, ast.Name)
            and fn.value.id == "logger"
            and fn.attr in ("info", "warning", "error", "debug", "critical")
        ):
            yield node


def test_phase_log_labels_are_sourced_from_the_table():
    offenders = {}
    table_sourced = 0
    for path in _python_files():
        if path.name == "phases.py":
            continue  # The table itself is where the literals legitimately live
        text = path.read_text(encoding="utf-8")
        for call in _log_call_label_args(ast.parse(text)):
            hand_written = _label_source(call.args[0])
            if hand_written:
                offenders.setdefault(str(path.relative_to(PACKAGE)), []).append(
                    (call.lineno, hand_written)
                )
            elif "canonical_log_tag" in ast.dump(call.args[0]):
                table_sourced += 1
    assert table_sourced >= 20, (
        f"only {table_sourced} labels come from the table; the scanner is not "
        "measuring what it claims to"
    )
    assert not offenders, (
        f"log labels must be built by canonical_log_tag(), not typed out: {offenders}"
    )


def test_control_the_log_label_scanner_fires_on_both_label_forms():
    hand_written = ast.parse(
        "logger.info('[Phase 9] nope')\n"
        'logger.info(f"[Phase 9] {n} tips")\n'
        "logger.info('Phase 9: nope')\n"
    )
    labels = [_label_source(call.args[0]) for call in _log_call_label_args(hand_written)]
    assert labels == ["[Phase 9]", "[Phase 9]", "Phase 9:"], labels

    from_table = ast.parse('logger.info(f"{canonical_log_tag(chr(57))} nope")')
    assert [_label_source(c.args[0]) for c in _log_call_label_args(from_table)] == [
        None
    ], "a table-sourced label must not be flagged"


def test_the_table_has_runtime_callers():
    """``phases.py`` declared itself the single source of truth while no
    module imported it. A table nobody calls cannot stay in sync by itself."""
    callers = [
        path.name
        for path in _python_files()
        if path.name != "phases.py" and "canonical_log_tag" in path.read_text(encoding="utf-8")
    ]
    assert len(callers) >= 6, f"expected the table to be used across the pipeline, got {callers}"
    assert {"pipeline.py", "main.py"} <= set(callers), callers
