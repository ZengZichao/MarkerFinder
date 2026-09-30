"""The capability surface the validation suite must cover, read from the code.

Why this module exists
----------------------
"Test everything the software claims" is only checkable if the claim is
enumerated rather than remembered. So the inventory is derived:

* ``cli_options`` — every option argparse accepts (``markerfinder --help``)
* ``subcommands`` — the four resumable steps
* ``exit_codes`` — the documented terminal states
* ``products`` — the files a run is supposed to leave behind

and ``test_v90_capability_matrix.py`` fails when any of them carries no case.
Cases declare what they cover with ``@pytest.mark.capability("<dest>")``; the
markers are collected at import/collection time into ``.work/capabilities.json``
by ``conftest.py``, so the matrix is built from the tests themselves and cannot
drift from them.

Anything listed in ``OPTION_ROLES`` with a role other than ``"positive"`` is an
option whose behaviour is deliberately exercised as a refusal, a deprecation
notice, or a documented no-op — the role is recorded in the exported matrix so a
reader can see *how* a capability was covered, not only that it was.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(__file__).parent.resolve()
SNAPSHOT_DIR = ROOT / ".work" / "capability_snapshot"

# Options that cannot be exercised by adding them to a normal run, grouped by
# Why. Keys are argparse dests; values are one-line rationales that end up in
# The exported matrix.
OPTION_ROLES: Dict[str, str] = {
    # ── plain behavioural knobs ────────────────────────────────────────────
    "input": "positive", "output": "positive", "threads": "positive",
    "mode": "positive", "marker_mode": "positive",
    "marker_hmm_dir": "positive", "marker_db_source": "positive",
    "min_hmm_score": "positive", "gtdb_markers_dir": "positive",
    "species_tree": "positive", "sequences": "positive",
    "taxonomy_table": "positive", "taxonomy_format": "positive",
    "taxonomy_source_priority": "positive",
    "taxonomy_delimiter_mode": "positive", "table_sep": "positive",
    "auto_embed_taxonomy": "positive", "multi_tree_mode": "positive",
    "strip_annotations": "positive", "mol_type": "positive",
    "skip_length_check": "positive", "no_cross_check": "positive",
    "ignore_malformed": "positive", "taxonomy_levels": "positive",
    "hgt_steps": "positive", "hgt_threshold_l1l2": "positive",
    "hgt_threshold_l2l3": "positive", "max_markers": "positive",
    "marker_preset": "positive", "hgt_scan": "positive",
    "scan_stability_min": "positive", "require_evidence_coverage": "positive",
    "hgt_mode": "positive", "consistency_stringency": "positive",
    "min_informative_sites": "positive", "cog_category_map": "positive",
    "strict_assertions": "positive", "allow_assertion_failure": "positive",
    "taxonomy_mustpass": "positive", "monophyly_rank": "positive",
    "monophyly_threshold": "positive", "ufboot": "positive",
    "gene_tree_builder": "positive", "coalescent_mode": "positive",
    "min_gene_tree_support": "positive", "checkm_results": "positive",
    "skip_checkm": "positive", "min_occupancy": "positive",
    "report_format": "positive", "force": "positive", "no_clobber": "positive",
    "save_intermediates": "positive", "redo": "positive", "resume": "positive",
    "config": "positive", "log_file": "positive", "tmp_dir": "positive",
    "keep_tmp": "positive", "db_dir": "positive", "verbose": "positive",
    "hgt_adaptive_thresholds": "positive",
    # ── deprecated aliases: the claim is the notice, not a new behaviour ────
    "tree": "deprecated-alias",
    "hgt_threshold": "deprecated-alias",
    "min_marker_coverage": "deprecated-alias",
    "fast_tree": "deprecated-alias",
    "no_hgt_adaptive_thresholds": "documented-no-op",
    # ── terminal-state switches ────────────────────────────────────────────
    "version": "terminal-flag", "self_test": "terminal-flag",
    "command": "subcommand-dispatch",
}

# The four resumable steps, plus the aggregate run with no subcommand.
SUBCOMMAND_CAPABILITIES = [
    "subcommand:scan", "subcommand:filter", "subcommand:infer",
    "subcommand:report", "subcommand:aggregate-run",
]

# Output products, by the directory a run writes them into. Names are the
# Defaults (prefix ``markerfinder``); a case that runs with another prefix
# Asserts the same suffix.
PRODUCTS = [
    "product:Phase5_reports.report.html",
    "product:Phase5_reports.marker_summary.tsv",
    "product:Phase5_reports.hgt_evaluation.tsv",
    "product:Phase5_reports.pipeline_summary.txt",
    "product:Phase5_reports.assertions.tsv",
    "product:Phase5_reports.threshold_scan.tsv",
    "product:Phase5_reports.composition.tsv",
    "product:Phase5_reports.excluded_profile.tsv",
    "product:Phase5_reports.mustpass.tsv",
    "product:Phase5_evidence.decision.json",
    "product:Phase5_metadata.run_config.json",
    "product:Phase4_trees.species_tree_concat.newick",
    "product:Phase4_trees.species_tree_astral.newick",
    "product:Phase4_trees.gene_trees.newick",
    "product:Phase4_trees.gene_trees.cache",
    "product:Phase4_alignments.partition.nex",
    "product:Phase4_intermediate.saved",
    "product:.markerfinder.pipeline_state.json",
    "product:.markerfinder.context.json",
]

# Behavioural capabilities that are not a single option or file.
WORKFLOW_CAPABILITIES = [
    "workflow:occupancy-selection",
    "workflow:hgt-grading",
    "workflow:state-persistence",
    "workflow:determinism",
    "workflow:provenance-recorded",
    "workflow:unknown-not-placeholder",
    "workflow:failure-loudness",
    "data:provenance-integrity",
    "data:manifest-checksums",
]


def cli_options() -> List[Tuple[str, str]]:
    """``(option strings, capability id)`` for every option the parser accepts.

    A deprecated alias that stores into another option's dest (``--tree`` ->
    ``species_tree``) is still advertised in ``--help`` and still carries its own
    behavioural claim (the notice *is* the behaviour), so its recording
    attribute is reported as an id beside the canonical dest — otherwise the
    alias would be testable but unnameable in the matrix.
    """
    from markerfinder.cli.parser import _build_parser

    parser = _build_parser()
    out = []
    for action in parser._actions:              # Noqa: SLF001 - introspection
        if not action.option_strings or action.dest == "help":
            continue
        joined = ",".join(action.option_strings)
        out.append((joined, action.dest))
        if getattr(action, "alias_attr", None):
            # Named after the option the user types (--tree -> 'tree'), not
            # After the attribute that records its use: the capability being
            # Advertised is the option, and cases declare it by that name.
            alias_id = action.option_strings[0].lstrip("-").replace("-", "_")
            out.append((joined, alias_id))
    return out


def subcommands() -> List[str]:
    from markerfinder.cli.parser import _SUBCOMMAND_SPECS

    return [f"subcommand:{name}" for name, _h, _d in _SUBCOMMAND_SPECS]


def exit_codes() -> List[str]:
    from markerfinder.cli import constants

    return sorted(
        name for name in vars(constants)
        if name.startswith("EXIT_") and isinstance(getattr(constants, name), int)
    )


def required_capabilities() -> Dict[str, str]:
    """Capability id -> where it came from (used in the exported matrix)."""
    required: Dict[str, str] = {}
    for _opts, dest in cli_options():
        required[dest] = OPTION_ROLES.get(dest, "positive")
    for cap in SUBCOMMAND_CAPABILITIES + PRODUCTS + WORKFLOW_CAPABILITIES:
        required[cap] = "positive"
    for code in exit_codes():
        required[f"exit:{code}"] = "terminal-state"
    return required


def collected_case_count() -> int:
    """How many items the collection snapshot was taken over.

    A guard against a partial run reporting full coverage: ``-k somefilter``
    shrinks the snapshot, and a matrix built from it would silently claim that
    the missing cases were never needed.
    """
    counts = sorted(SNAPSHOT_DIR.glob("*.count"))
    if not counts:
        return 0
    return max(int(p.read_text(encoding="utf-8")) for p in counts)


def covered_capabilities() -> Dict[str, List[str]]:
    """Capability id -> test node ids that declare it, from the collection dump."""
    files = sorted(SNAPSHOT_DIR.glob("*.json"))
    if not files:
        raise FileNotFoundError(
            f"no collection snapshot under {SNAPSHOT_DIR}: it is written during "
            "pytest collection (conftest.pytest_collection_modifyitems). Run "
            "`python -m pytest validation/cases --collect-only -q` first."
        )
    merged: Dict[str, List[str]] = {}
    for path in files:
        mapping = json.loads(path.read_text(encoding="utf-8"))
        for capability, node_ids in mapping.items():
            merged.setdefault(capability, [])
            merged[capability].extend(
                n for n in node_ids if n not in merged[capability]
            )
    return {k: sorted(v) for k, v in merged.items()}


def export_matrix(path: Path) -> List[str]:
    """Write the coverage matrix as TSV; returns the uncovered capability ids."""
    required = required_capabilities()
    covered = covered_capabilities()
    lines = ["capability\torigin\tcovered_by_cases"]
    missing: List[str] = []
    for capability in sorted(required):
        cases = covered.get(capability, [])
        if not cases:
            missing.append(capability)
        lines.append(
            "\t".join([capability, required[capability],
                       ";".join(sorted(cases)) or "NOT COVERED"])
        )
    # Anything a case declared but that does not exist in the surface is a
    # Typo or a capability that was removed from the code while its marker
    # Stayed behind: both must surface.
    for capability in sorted(set(covered) - set(required)):
        lines.append(f"{capability}\tUNDECLARED\t"
                     f"{';'.join(sorted(covered[capability]))}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return missing
