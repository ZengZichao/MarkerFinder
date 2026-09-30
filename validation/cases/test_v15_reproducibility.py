"""V-15 — reproducibility: the same input must give the same answer.

Two claims are tested, and they are different claims:

*Determinism* — running the pipeline twice on the same inputs, in the same
environment, produces the same scientific content. Timestamps, run durations,
absolute paths and temp directories legitimately differ, so they are excluded
by name rather than by ignoring the whole file.

*Replay* — a recorded ``run_config.json`` reproduces the recorded
configuration (the CLI case V-02 covers the parameter layer; this one compares
the products).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

# Lines/fields that legitimately differ between two runs of the same input.
# Not anchored to line starts: the same values appear inline in the HTML report
# (``<p>Generated: …</p>``), and a per-case exception for "the report has a
# Timestamp" would defeat the point of comparing it at all.
VOLATILE_PATTERNS = (
    re.compile(r"Generated:\s*\d{4}-\d{2}-\d{2}T[\d:.]+"),
    re.compile(r"Runtime:\s*[\d.]+s"),
    re.compile(r"(?i)(input|output) directory:\s*[^\n<]*"),
    re.compile(r"(?i)timestamp[^\n<]*"),
    # The automatic temp directory carries a per-run random suffix, and it is
    # Written into run_config.json both absolutely and relative to the output
    # Directory; only the suffix differs, so normalise the token itself.
    re.compile(r"markerfinder-[0-9a-f]{8}"),
)
VOLATILE_KEYS = {
    "timestamp", "run_duration_seconds", "markerfinder_version",
}


def _scrub_text(text: str) -> str:
    """Apply the line-level volatile patterns (non-JSON products).

    Only for text files: ``timestamp[^\n<]*`` swallows the rest of a JSON line
    including the closing quote and comma, which makes the document
    unparseable — and an unparseable snapshot falls back to raw text, where the
    per-run ``run_duration_seconds`` survives and every run "differs".
    """
    for pattern in VOLATILE_PATTERNS:
        text = pattern.sub("<VOLATILE>", text)
    return text


def _normalise(path: Path, root: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    # Absolute paths differ because each run gets its own output directory.
    text = text.replace(str(root), "<RUN>")
    if path.suffix == ".json":
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            return _scrub_text(text)
        _strip_volatile(payload)
        return json.dumps(payload, sort_keys=True, indent=1)
    return _scrub_text(text)


def _scrub(value):
    """Normalise a leaf string: drop per-run path/temp tokens."""
    if not isinstance(value, str):
        return value
    for pattern in VOLATILE_PATTERNS:
        value = pattern.sub("<VOLATILE>", value)
    return value


def _strip_volatile(node):
    if isinstance(node, dict):
        for key in list(node):
            if key in VOLATILE_KEYS or "_dir" in key or "tmp" in key:
                node[key] = "<VOLATILE>"
            elif isinstance(node[key], str):
                node[key] = _scrub(node[key])
            else:
                _strip_volatile(node[key])
    elif isinstance(node, list):
        for i, item in enumerate(node):
            if isinstance(item, str):
                node[i] = _scrub(item)
            else:
                _strip_volatile(item)


@pytest.mark.capability("workflow:determinism",
                        "product:Phase5_reports.hgt_evaluation.tsv",
                        "product:Phase4_trees.species_tree_concat.newick")
def test_two_runs_produce_identical_scientific_content(mf, case_output,
                                                      record_metric):
    a = mf(out_dir=case_output / "a", extra=["--force"])
    a.assert_ok("first run")
    b = mf(out_dir=case_output / "b", extra=["--force"])
    b.assert_ok("second run")

    compared = []
    differing = []
    for path in sorted(a.out_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = str(path.relative_to(a.out_dir))
        if rel.startswith((".markerfinder/", "Phase4_intermediate/")):
            continue          # Internal caches, not reported content
        other = b.out_dir / rel
        if not other.exists():
            differing.append(f"{rel} (missing in the second run)")
            continue
        compared.append(rel)
        if _normalise(path, a.out_dir) != _normalise(other, b.out_dir):
            differing.append(rel)

    missing_in_b = [
        str(p.relative_to(b.out_dir)) for p in b.out_dir.rglob("*")
        if p.is_file() and not str(p.relative_to(b.out_dir)).startswith(
            (".markerfinder/", "Phase4_intermediate/"))
        and not (a.out_dir / p.relative_to(b.out_dir)).exists()
    ]
    record_metric("v15_determinism", "files_compared", len(compared))
    record_metric("v15_determinism", "files_differing", differing)
    assert not differing, (
        f"two runs on the same input differ in {differing} "
        f"(of {len(compared)} comparable files)"
    )
    assert not missing_in_b, f"only the first run wrote: {missing_in_b}"
    assert len(compared) >= 6, (
        f"only {len(compared)} files were comparable — the check is vacuous"
    )


@pytest.mark.capability("workflow:determinism",
                        "product:Phase5_evidence.decision.json")
def test_marker_ranking_is_identical_across_two_runs(mf, read_tsv, case_output,
                                                    record_metric):
    """The tree files can legitimately differ in float formatting while the
    decision stays the same; the ranking is the decision, so it is compared
    directly."""
    a = mf(out_dir=case_output / "rank_a", extra=["--force"])
    b = mf(out_dir=case_output / "rank_b", extra=["--force"])
    a.assert_ok()
    b.assert_ok()
    ra = read_tsv(a.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    rb = read_tsv(b.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    ranking_a = [(r["marker_id"], r["hgt_risk_level"], r["overall_risk"]) for r in ra]
    ranking_b = [(r["marker_id"], r["hgt_risk_level"], r["overall_risk"]) for r in rb]
    record_metric("v15_ranking", "markers", len(ranking_a))
    assert ranking_a == ranking_b, (
        f"marker decisions differ between runs:\n{ranking_a}\n{ranking_b}"
    )


@pytest.mark.capability("workflow:determinism",
                        "product:Phase5_metadata.run_config.json")
def test_seed_and_support_values_are_pinned_not_left_to_the_clock(mf,
                                                                 record_metric):
    """Deterministic tree building requires the external tools to be seeded.

    The recorded configuration is where a reader checks that: UFBOOT replicates
    and the builder are both written, and the pipeline passes ``-seed`` to
    IQ-TREE and FastTree (asserted here from the products, and in the unit
    suite from the call sites).
    """
    run = mf(extra=["--coalescent-mode", "always", "--force"])
    run.assert_ok()
    assert run.recorded("phylo_config", "ufboot_replicates"), \
        run.recorded("phylo_config", "ufboot_replicates")
    assert run.recorded("phylo_config", "gene_tree_builder") == "fasttree"
    trees = run.out_dir / "Phase4_trees" / "gene_trees"
    assert trees.is_dir() and list(trees.glob("*.nwk")), (
        "no per-marker gene trees to be reproduced"
    )
    record_metric("v15_seeds", "gene_tree_cache", len(list(trees.glob("*.nwk"))))


@pytest.mark.capability("workflow:determinism", "config",
                        "product:Phase5_metadata.run_config.json")
def test_a_replayed_run_reproduces_the_products(mf, case_output, read_tsv,
                                               record_metric):
    original = mf(out_dir=case_output / "orig", extra=["--force"])
    original.assert_ok()
    recorded = original.product("Phase5_metadata/run_config.json")
    replay = mf(omit=("--marker-mode", "--gtdb-markers-dir", "--taxonomy-table",
                      "--max-markers", "--monophyly-rank",
                      "--gene-tree-builder", "--coalescent-mode",
                      "--skip-checkm"),
                out_dir=case_output / "replay",
                extra=["--config", str(recorded), "--force"])
    replay.assert_ok("replayed run")

    before = read_tsv(original.product("Phase5_reports/markerfinder.marker_summary.tsv"))
    after = read_tsv(replay.product("Phase5_reports/markerfinder.marker_summary.tsv"))
    ids_before = {r["marker_id"] for r in before}
    ids_after = {r["marker_id"] for r in after}
    record_metric("v15_replay", "markers_before", sorted(ids_before))
    record_metric("v15_replay", "markers_after", sorted(ids_after))
    assert ids_before == ids_after, (
        f"replaying the recorded configuration selected a different marker set: "
        f"{sorted(ids_before ^ ids_after)}"
    )
