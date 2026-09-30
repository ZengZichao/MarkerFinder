"""User-facing help text may not name a path the program never writes.

The symptom this gate was written for: ``--help`` told the user that state lives
in ``<output>/.markerfinder/.pipeline_state.pkl`` while ``state_codec`` has been
writing ``.pipeline_state.json`` (schema v2) the whole time. Three copies of the
wrong string were live in the shipped help text and the ``main`` dispatch
comment. The number of tests in this repository is irrelevant to that bug: a
user who cannot find the file the manual promised will reasonably suspect the
install, and the author is the one who pays for the debugging.

 already locks phase numbers in help text; nothing locked *paths*, so this
file adds the missing rule in the repo's existing shape -- a single-source
comparison rather than a copy-editing one:

    every ``<something>/<file>.<ext>`` literal in the CLI help surface must name
    a basename that some module in the package actually writes or reads.

Every check is two-sided, per the convention in this suite: the scanner is
proved able to fire on a planted path before it is allowed to report clean, so
"grep found nothing" cannot come from a broken pattern.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PACKAGE = REPO / "markerfinder"

# The user-facing surface: argparse help strings and the dispatch comments that
# Are mirrored into them.
HELP_SOURCES = [
    PACKAGE / "cli" / "parser.py",
    PACKAGE / "cli" / "main.py",
]

# ``<path/to/name.ext>`` -- only names that carry an extension are claims about a
# File the program produces; bare directory mentions are checked by.
HELP_PATH = re.compile(r"[\w.<>/-]*?/?([\w.-]+\.(?:json|pkl|tsv|csv|txt|html|htm|domtblout|fasta|faa|fna|nwk|newick))\b")

# A name counts as "real" when the package itself spells it out in a string
# Literal (write side, read side, or a constant such as STATE_FILE_NAME). Names
# Built with a placeholder -- ``f"{prefix}.threshold_scan.tsv"`` -- register both
# Their full form and their trailing component, which is what help text quotes.
ANY_LITERAL = re.compile(r'"([^"\n]*)"|\'([^\'\n]*)\'')
FILE_TAIL = re.compile(
    r"([\w{}.\-/]*\.(?:json|pkl|tsv|csv|txt|html|htm|domtblout|fasta|faa|fna|nwk|newick))\b"
)


def real_output_names() -> set[str]:
    """Every file name the package itself writes or reads, in claim-friendly form."""
    names: set[str] = set()
    for path in PACKAGE.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        for literal in ANY_LITERAL.findall(text):
            for token in FILE_TAIL.findall(literal[0] or literal[1]):
                names.add(token)
                tail = token.rsplit(".", 2)
                if len(tail) == 3 and tail[1]:
                    names.add(f"{tail[1]}.{tail[2]}")
        # Values reached through a constant (STATE_FILE_NAME = ".pipeline_state.json").
        names.update(re.findall(r'_NAME\s*[:=]\s*["\']([\w.-]+\.\w+)["\']', text))
    return names


def help_paths(path: Path) -> list[str]:
    """File-name claims made by the help surface.

    Usage-example lines (``markerfinder -i... --species-tree tree.nwk``) are
    skipped on purpose: they name *inputs the user supplies*, whose names are
    arbitrary by definition. Everything else in this file that mentions a file is
    a promise about what the run produces.
    """
    claims: list[str] = []
    for line in path.read_text(encoding="utf-8").split("\n"):
        stripped = line.strip()
        if stripped.startswith(("import ", "from ")) or "markerfinder " in stripped:
            continue
        claims.extend(HELP_PATH.findall(line))
    return claims


def test_control_the_path_scanner_fires_on_a_retired_filename():
    """Without this, a null scanner would report the real check as passing."""
    planted = 'state to <output>/.markerfinder/.pipeline_state.pkl for the next step'
    assert HELP_PATH.findall(planted) == [".pipeline_state.pkl"], "the F-1 bug must be catchable"
    assert ".pipeline_state.pkl" not in real_output_names(), (
        "if the package ever writes a .pkl again, this gate's premise needs "
        "re-reading rather than deleting"
    )


def test_control_the_scanner_ignores_prose_without_a_file_claim():
    text = "Requires a prior 'filter' step. Reuses cached gene trees."
    assert HELP_PATH.findall(text) == []


def test_control_placeholder_built_names_are_recognised():
    """``{prefix}.threshold_scan.tsv`` is written; help quoting the tail is fine,
    help quoting a name nobody writes is not."""
    real = real_output_names()
    assert "threshold_scan.tsv" in real, "the suffix of a placeholder-built name must count"
    assert "run_config.json" in real
    assert "no_such_artifact.tsv" not in real


def test_no_help_text_names_a_file_the_program_never_writes():
    real = real_output_names()
    offenders: dict[str, list[str]] = {}
    for source in HELP_SOURCES:
        missing = sorted(
            name for name in set(help_paths(source))
            if name not in real and not any(name in candidate for candidate in real)
        )
        if missing:
            offenders[str(source.relative_to(REPO))] = missing
    assert not offenders, (
        f"help text advertises files that appear nowhere in the package: {offenders}; "
        f"{len(real)} output names are actually written or read"
    )


def test_the_state_file_help_text_matches_the_codec_constant():
    """The specific regression, stated as a positive assertion."""
    from markerfinder.utils.state_codec import STATE_DIR_NAME, STATE_FILE_NAME

    combined = "\n".join(p.read_text(encoding="utf-8") for p in HELP_SOURCES)
    assert STATE_FILE_NAME in combined, "help never mentions the real state file"
    assert f"{STATE_DIR_NAME}/{STATE_FILE_NAME}" in combined, (
        "help must name the real directory + file pair the codec writes"
    )
    assert ".pipeline_state.pkl" not in combined, (
        "a retired pickle path is back in the user-facing help text"
    )
    assert STATE_FILE_NAME.endswith(".json"), (
        "schema v2 state is JSON; if that ever changes, update the help AND this"
    )
