"""V-17 — measured detection performance on the 12-genome benchmark.

The claim under test is MarkerFinder's second headline: "phylogenetic-only HGT
screening removes phylogenetically misleading markers". A claim about detection
needs a positive control with a known answer, and this bundle ships one.

Design (paired, so the background rate is measured and not assumed)
------------------------------------------------------------------
``markers/bench12`` the clean marker set: 20 markers × 12 genomes.
``markers/chimera_bench12`` byte-identical except that in three markers the
                              *Bacillus anthracis* Ames (Caryophanales)
                              sequence is replaced by the *Lactococcus lactis*
                              Il1403 (Lactobacillales) ortholog — a cross-order
                              transfer at a known position.

Both runs see the same genomes, the same taxonomy table and the same
parameters, so a marker whose risk rises between them rose because of the
planted sequence, not because of the dataset. The clean run is therefore used
as the negative control: markers it excludes are the background false positives.

Everything measured is written to ``results/chimera_metrics.json`` and quoted by
the test report, so the report's numbers come from the run rather than from
whoever wrote the report.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.slow

# The whole spanning directory must reach the screen: with a budget below the
# Number of candidates, some planted positives are never evaluated at all, and
# "every planted marker raised its risk" would then be a claim about the subset
# The budget happened to keep (measured: 12 of the 20 candidates scored, and 2 of
# The 3 planted markers absent from both runs).
BUDGET = ["--max-markers", "20", "--monophyly-rank", "family",
          "--gene-tree-builder", "fasttree", "--coalescent-mode", "off"]


@pytest.fixture(scope="module")
def benchmark_runs(request, mf, data_dir, genome_input):
    """The two runs this file compares. Module-scoped: they are the expensive
    part of the whole suite (12 real genomes)."""
    out_root = Path(request.config.rootdir) / "validation" / ".work" / "outputs"
    clean = mf(genome_input("bench12"),
               markers=data_dir / "markers" / "bench12",
               taxonomy=data_dir / "taxonomy" / "taxonomy_bench12.tsv",
               extra=BUDGET + ["--force"],
               out_dir=out_root / "v17_clean")
    clean.assert_ok("clean benchmark run")
    chimera = mf(genome_input("bench12"),
                 markers=data_dir / "markers" / "chimera_bench12",
                 taxonomy=data_dir / "taxonomy" / "taxonomy_bench12.tsv",
                 extra=BUDGET + ["--force"],
                 out_dir=out_root / "v17_chimera")
    chimera.assert_ok("planted-chimera benchmark run")
    return clean, chimera


def _evaluation(run, read_tsv):
    rows = read_tsv(run.product("Phase5_reports/markerfinder.hgt_evaluation.tsv"))
    return {r["marker_id"]: r for r in rows}


@pytest.mark.capability("workflow:hgt-grading", "workflow:failure-loudness")
def test_planted_chimeras_raise_the_risk_of_exactly_their_own_markers(
        benchmark_runs, read_tsv, data_dir, record_metric):
    """The directional core of the detection claim.

    For each planted marker the risk in the chimera run is compared with the
    risk of the SAME marker in the clean run; the unplanted markers are the
    internal negative control and must not move.
    """
    clean, chimera = benchmark_runs
    truth = json.loads((data_dir / "markers" /
                        "chimera_ground_truth_bench12.json")
                       .read_text(encoding="utf-8"))
    planted = {c["marker_id"] for c in truth["planted"]}
    assert planted, "the ground truth lists no planted chimeras"

    before, after = _evaluation(clean, read_tsv), _evaluation(chimera, read_tsv)
    shared = set(before) & set(after)
    assert shared, "the two runs share no scored marker"

    def risk(row):
        try:
            return float(row["overall_risk"])
        except (TypeError, ValueError):
            return None

    moved_up, unmoved, bystanders = [], [], []
    for marker in sorted(shared):
        a, b = risk(before[marker]), risk(after[marker])
        if a is None or b is None:
            continue
        if marker in planted:
            (moved_up if b > a else unmoved).append((marker, a, b))
        elif b != a:
            bystanders.append((marker, a, b))

    metrics = {
        "planted": sorted(planted),
        "planted_risk_in_clean_run": {m: a for m, a, _b in moved_up + unmoved},
        "planted_risk_in_chimera_run": {m: b for m, _a, b in moved_up + unmoved},
        "planted_with_risk_increase": sorted(m for m, _a, _b in moved_up),
        "planted_unchanged_or_lower": sorted(m for m, _a, _b in unmoved),
        "unplanted_markers_that_moved": sorted(bystanders),
        "n_markers_scored_in_both": len(shared),
    }
    record_metric("v17_detection", "planted_risk_raised",
                  len(metrics["planted_with_risk_increase"]))
    record_metric("v17_detection", "planted_total", len(planted))
    # Merged, not written: the two cases in this file publish different views of
    # The same paired runs, and whichever finished last used to erase the other.
    _merge_metrics(data_dir.parent / "results" / "chimera_metrics.json",
                   {"directionality": metrics})

    assert metrics["planted_with_risk_increase"], (
        "not one planted cross-order transfer raised its marker's HGT risk; the "
        f"phylogenetic screen is blind on this data: {metrics}"
    )
    # Every planted marker must reach the maximum risk in the chimeric run: a
    # Positive that is not excluded is a miss, whatever its baseline was.
    not_maxed = {m: r for m, r in
                 metrics["planted_risk_in_chimera_run"].items() if r < 1.0}
    assert not not_maxed, (
        f"planted chimeras left unexcluded in the chimeric run: {not_maxed}"
    )
    # A rise is only defined against the baseline, and the baseline is measured,
    # Not assumed: a marker already flagged in the CLEAN run cannot rise, and
    # Counting that as a detection failure would hide the real number — the
    # Screen's false-positive rate on this data. It is recorded and bounded here
    # Instead of being wished away.
    already_flagged = sorted(
        m for m in planted
        if metrics["planted_risk_in_clean_run"].get(m, 0.0) >= 1.0
    )
    record_metric("v17_detection", "planted_already_flagged_clean", already_flagged)
    assert len(already_flagged) < len(planted), (
        f"the clean run already maxed out every planted marker "
        f"({already_flagged}): the paired design has no negative direction left "
        "to measure, so this comparison cannot say anything about detection"
    )
    assert not bystanders, (
        f"unplanted markers changed between two runs whose only difference is "
        f"three planted sequences: {bystanders}"
    )


@pytest.mark.capability("workflow:hgt-grading")
def test_detection_quality_at_the_exclusion_decision(benchmark_runs, read_tsv,
                                                    data_dir, record_metric):
    """Precision / recall / F1 where a "detection" means the marker is excluded
    (Level 3), reported against both runs.

    Level 2 is deliberately NOT counted as a detection: it is the suspicious
    band that the pipeline keeps, so scoring it as a hit would inflate the
    apparent sensitivity of the screen. The clean run gives the background
    exclusion rate, i.e. the false-positive side of the same decision.
    """
    clean, chimera = benchmark_runs
    truth = json.loads((data_dir / "markers" /
                        "chimera_ground_truth_bench12.json")
                       .read_text(encoding="utf-8"))
    planted = {c["marker_id"] for c in truth["planted"]}

    before, after = _evaluation(clean, read_tsv), _evaluation(chimera, read_tsv)
    excluded_clean = {m for m, r in before.items() if r["hgt_risk_level"] == "level_3"}
    excluded_chimera = {m for m, r in after.items() if r["hgt_risk_level"] == "level_3"}

    # Differential call: excluded because of the planted sequence.
    called = excluded_chimera - excluded_clean
    tp = called & planted
    fp = called - planted
    fn = planted - called
    precision = len(tp) / len(called) if called else 0.0
    recall = len(tp) / len(planted) if planted else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) else 0.0)

    metrics = {
        "decision": "hgt_risk_level == level_3",
        "n_planted": len(planted),
        "n_excluded_in_clean_run": len(excluded_clean),
        "n_excluded_in_chimera_run": len(excluded_chimera),
        "background_exclusion_rate": round(
            len(excluded_clean) / len(before), 4) if before else None,
        "differentially_excluded": sorted(called),
        "true_positives": sorted(tp),
        "false_positives": sorted(fp),
        "missed": sorted(fn),
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
    }
    record_metric("v17_quality", "precision", precision)
    record_metric("v17_quality", "recall", recall)
    record_metric("v17_quality", "background_exclusion_rate",
                  metrics["background_exclusion_rate"])
    _merge_metrics(data_dir.parent / "results" / "chimera_metrics.json",
                   {"quality_at_exclusion": metrics})
    assert tp, f"no planted chimera was differentially excluded: {metrics}"
    assert precision > 0.0, metrics


@pytest.mark.capability("data:provenance-integrity",
                        "workflow:provenance-recorded")
def test_chimera_ground_truth_matches_the_shipped_marker_files(data_dir):
    """The positive control is only a control if the swap is actually there.

    Checked against the files rather than against the JSON alone: a ground truth
    that describes a directory nobody built would make every number in the next
    two cases meaningless.
    """
    truth_path = data_dir / "markers" / "chimera_ground_truth_bench12.json"
    assert truth_path.exists(), f"{truth_path} missing — re-run 02_build_marker_sets.py"
    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    clean_dir = data_dir / "markers" / "bench12"
    chimera_dir = data_dir / "markers" / "chimera_bench12"
    recipient = truth["recipient"]["accession"]
    donor = truth["donor"]["accession"]

    for entry in truth["planted"]:
        marker = entry["marker_id"]
        clean_lines = (clean_dir / f"{marker}.faa").read_text(
            encoding="utf-8").split("\n")
        chimeric = (chimera_dir / f"{marker}.faa").read_text(encoding="utf-8")
        assert f">{recipient}" in chimeric, f"{marker}: recipient absent"

        def seq_of(lines, name):
            for i, line in enumerate(lines):
                if line == f">{name}":
                    return lines[i + 1]
            return None

        clean_recipient = seq_of(clean_lines, recipient)
        clean_donor = seq_of(clean_lines, donor)
        chimeric_recipient = seq_of(chimeric.split("\n"), recipient)
        assert chimeric_recipient == clean_donor, (
            f"{marker}: the recipient sequence is not the donor's, so this is "
            "not the planted transfer the ground truth describes")
        assert chimeric_recipient != clean_recipient, (
            f"{marker}: recipient unchanged")
        assert entry["donor_sequence_length"] == len(clean_donor)

    # And the untouched markers must really be untouched, or the paired design
    # Has no negative control.
    planted = {e["marker_id"] for e in truth["planted"]}
    changed = []
    for faa in sorted(clean_dir.glob("*.faa")):
        if faa.stem in planted:
            continue
        other = chimera_dir / faa.name
        if not other.exists() or other.read_bytes() != faa.read_bytes():
            changed.append(faa.name)
    assert not changed, f"non-planted marker files differ between runs: {changed}"


def _merge_metrics(path: Path, extra: dict) -> None:
    payload = {}
    if path.exists():
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}
    payload.update(extra)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
