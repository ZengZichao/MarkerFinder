"""Product-level acceptance evidence from a real end-to-end run.

The same mocked external tools the existing integration suite relies on, so what
this proves is that the *pipeline writes the promised artifacts with the promised
columns and states* -- not that the biology is right (that needs real
mafft/FastTree/ASTRAL output on real genomes).

Why it matters: a channel can be implemented, unit-tested in isolation, and then
never actually exercised by the run path -- so the only honest check is to run
the pipeline end to end and read the files it leaves behind.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from markerfinder.config import PipelineConfig, SelectionConfig
from markerfinder.pipeline import MarkerFinderPipeline
from markerfinder.utils.io import load_genomes_from_directory

from .test_pipeline import _mock_subprocess_run, _tree_newick

pytestmark = pytest.mark.integration

GENOME_IDS = ["g1", "g2", "g3", "g4"]
MARKERS = ("M1", "M2", "M3", "M4")

CARD_KEYS = (
    "schema_version", "marker_id", "detector", "gene_tree_source",
    "trimming_regime", "rf", "quartet", "monophyly", "rank_used", "n_total",
    "n_mono", "weights", "n_signals_used", "thresholds_in_effect",
    "far_active", "risk_basis", "overall_risk", "level", "confidence",
    "assertion_ids_fired", "notes",
)

FOUR_STATES = {"measured", "not_measurable", "not_applicable", "rejected"}


def _write_inputs(tmp_path: Path) -> Path:
    indir = tmp_path / "input"
    indir.mkdir(parents=True, exist_ok=True)
    for gid in GENOME_IDS:
        (indir / f"{gid}.faa").write_text(
            ">p1\nACDEFGHIKLMNPQRSTVWY\n>p2\nACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8"
        )
    gtdb = tmp_path / "gtdb_markers"
    gtdb.mkdir(exist_ok=True)
    for marker in MARKERS:
        body = "\n".join(f">{gid}\nACDEFGHIKLMNPQRSTVWY" for gid in GENOME_IDS)
        (gtdb / f"{marker}.faa").write_text(body + "\n", encoding="utf-8")
    tree = tmp_path / "species.nwk"
    tree.write_text(_tree_newick(GENOME_IDS), encoding="utf-8")
    return indir


def _taxa():
    return {
        gid: {"domain": "Bacteria", "phylum": "P1", "class": "C1",
              "order": order, "family": family, "genus": genus}
        for gid, (order, family, genus) in zip(
            GENOME_IDS,
            [("O1", "F1", "G1"), ("O1", "F1", "G1"), ("O1", "F2", "G2"),
             ("O2", "F3", "G3")],
        )
    }


def _run(tmp_path: Path, **tweaks) -> Path:
    """One full pipeline run with the external tools mocked, modes configurable."""
    indir = _write_inputs(tmp_path)
    out = tmp_path / "output"
    cfg = PipelineConfig(
        input_dir=str(indir),
        output_dir=str(out),
        tmp_dir=str(tmp_path / "tmp"),
        output_prefix="markerfinder",
        selection_config=SelectionConfig(
            marker_mode="gtdb_tk",
            gtdb_markers_dir=str(tmp_path / "gtdb_markers"),
            species_tree=str(tmp_path / "species.nwk"),
            min_occupancy=0.0,
            max_markers=10,
        ),
    )
    cfg.phylo_config.fast_mode = True
    cfg.phylo_config.use_fasttree = True
    cfg.phylo_config.coalescent_mode = "post-filter"
    cfg.report_config.output_dir = cfg.output_dir
    cfg.report_config.output_prefix = cfg.output_prefix

    cfg.hgt_config.hgt_mode = tweaks.get("hgt_mode", "risk")
    cfg.hgt_config.min_informative_sites = int(tweaks.get("min_informative_sites", 0))
    if tweaks.get("composition"):
        cfg.hgt_config.hgt_steps = "phylogenetic,composition"
        # Config_build derives this pair from --hgt-steps; a programmatic run has
        # To set both, which is exactly the coupling the composition column keys on.
        cfg.report_config.composition_screen = True
    # 's switch lives on ReportConfig, not HGTConfig. Setting the wrong one
    # Would fail silently, so the scan flag is written where the scanner reads it.
    cfg.report_config.hgt_scan = bool(tweaks.get("hgt_scan", False))
    if tweaks.get("mustpass"):
        # Reads the requirements file off ReportConfig; the CLI fills it
        # From --taxonomy-mustpass, so a programmatic run sets the same field.
        cfg.report_config.taxonomy_mustpass = str(tweaks["mustpass"])

    genomes = load_genomes_from_directory(str(indir))
    with patch("subprocess.run",
               side_effect=_mock_subprocess_run(tmp_path, GENOME_IDS)), \
            patch("markerfinder.utils.gtdb_tk_markers.shutil.which",
                  return_value="/usr/bin/tool"):
        MarkerFinderPipeline(cfg).run(genomes, table_taxa=_taxa())
    return out


@pytest.fixture
def risk_run(tmp_path):
    return _run(tmp_path)


def _report(out: Path, name: str) -> list:
    return (out / "Phase5_reports" / name).read_text(
        encoding="utf-8"
    ).strip().split("\n")


# ──: the card and the tables, as written to disk ────────────

def test_decision_cards_are_written_with_the_required_fields(risk_run):
    cards = sorted((risk_run / "Phase5_evidence").glob("decision_*.json"))
    assert cards, "every graded marker must have one machine-readable card"
    payload = json.loads(cards[0].read_text(encoding="utf-8"))
    missing = [key for key in CARD_KEYS if key not in payload]
    assert not missing, (missing, sorted(payload))
    assert payload["schema_version"] == 2, payload
    # Every evidence slot declares one of the four states. The state is
    # What makes a null interpretable, and it has to be in the product.
    for slot in ("rf", "quartet", "monophyly"):
        assert payload[slot]["state"] in FOUR_STATES, (slot, payload[slot])


def test_ac05_every_non_unknown_marker_explains_itself(risk_run):
    """'s own criterion, read off the written cards: a marker that received a
    grade must say why. An empty note on a graded marker means the verdict cannot
    be audited from the product, which is what is for."""
    cards = [json.loads(c.read_text(encoding="utf-8"))
             for c in sorted((risk_run / "Phase5_evidence").glob("decision_*.json"))]
    assert cards
    graded = [c for c in cards if c["level"] != "UNKNOWN"]
    for card in graded:
        assert card["notes"], card["marker_id"]
        assert card["thresholds_in_effect"], card["marker_id"]
        assert card["weights"] is not None, card["marker_id"]


def test_hgt_table_carries_the_verdict_and_its_inputs(risk_run):
    header = _report(risk_run, "markerfinder.hgt_evaluation.tsv")[0].split("\t")
    for column in ("marker_id", "overall_risk", "hgt_risk_level", "notes",
                   "detector", "rank_used", "rf", "rf_state", "quartet",
                   "quartet_state", "weights", "thresholds_in_effect"):
        assert column in header, (column, header)


def test_marker_summary_columns_and_na_rendering(risk_run):
    lines = _report(risk_run, "markerfinder.marker_summary.tsv")
    header = lines[0].split("\t")
    assert header == ["marker_id", "occupancy_score", "marker_quality_score",
                      "marker_quality_level", "pis", "effective_columns",
                      "consistency_grade"], header
    assert len(lines) > 1, lines
    grade = header.index("consistency_grade")
    # The criterion did not run in the default mode, so the column
    # Reads NA -- never a borrowed risk verdict, never a 0.
    assert all(row.split("\t")[grade] == "NA" for row in lines[1:]), lines


# ── the two paths that crashed before any test reached them ───────────────

def test_pis_floor_runs_end_to_end(tmp_path):
    """``--min-informative-sites`` used to raise NameError in Phase 3.

    The floor reads the HGT evaluations, which inside ``run_infer`` live in
    ``filter_data`` rather than in a local called ``hgt_report``. No test had
    driven Phase 3 with a floor above zero, so the documented flag crashed the
    run on the way to writing nothing.
    """
    out = _run(tmp_path, min_informative_sites=1)
    assert _report(out, "markerfinder.marker_summary.tsv")[0].startswith("marker_id")
    assert (out / "Phase5_reports" / "markerfinder.pipeline_summary.txt").exists()


def test_threshold_scan_product_only_when_requested(tmp_path, risk_run):
    assert not list((risk_run / "Phase5_reports").glob("*threshold_scan*")), (
        "a default run must not grow a product nobody asked for"
    )
    scanned = _run(tmp_path / "scan", hgt_scan=True)
    produced = list((scanned / "Phase5_reports").glob("*threshold_scan*"))
    assert produced, "--hgt-scan must produce threshold_scan.tsv"
    assert len(produced[0].read_text(encoding="utf-8").strip().split("\n")) > 1


# ──: the two opt-in evidence channels ───────────────────────

def test_composition_product_is_parallel_and_na_where_not_applicable(tmp_path):
    lines = _report(_run(tmp_path, composition=True),
                    "markerfinder.composition.tsv")
    header = lines[0].split("\t")
    assert header[:5] == ["marker_id", "rcv", "gc_bias", "n_sequences",
                          "outlier_flag"], header
    # Read off the artifact instead of asserted about the code: the
    # Composition table carries no verdict column of any kind.
    assert not any(c in header for c in ("overall_risk", "overall_score", "level"))
    # G1: GC does not apply to protein input, so NA -- never 0.
    assert all(row.split("\t")[2] == "NA" for row in lines[1:]), lines


def test_consistency_product_carries_both_legs(tmp_path):
    out = _run(tmp_path, hgt_mode="consistency")
    lines = _report(out, "markerfinder.excluded_profile.tsv")
    header = lines[0].split("\t")
    assert header[-4:] == ["concat_strength", "concat_quartets",
                           "coalescent_strength", "coalescent_quartets"], header
    summary = (out / "Phase5_reports" / "markerfinder.pipeline_summary.txt").read_text(
        encoding="utf-8"
    )
    assert "RANKING ONLY" in summary, (
        "the run graded on consistency, so the summary must say which criterion "
        f"adjudicated; summary was:\n{summary[:1500]}"
    )
    assert "fr59:" in "\n".join(
        _report(out, "markerfinder.hgt_evaluation.tsv")
    ), "the demotion must also be per-marker evidence, not one summary sentence"


def test_hybrid_mode_writes_the_conjunction_into_the_products(tmp_path):
    """End to end: hybrid = risk AND consistency, both must pass.

    Read off the artifacts: the summary names the conjunction, the per-marker
    notes name which leg decided each marker, the grade column is populated,
    and the profile of the consistency-removed markers exists. (The mock-tool
    caveat from the module docstring applies here too: this proves the run
    path *records* the conjunction, not that the biology is right.)
    """
    out = _run(tmp_path, hgt_mode="hybrid")
    summary = (out / "Phase5_reports" / "markerfinder.pipeline_summary.txt").read_text(
        encoding="utf-8"
    )
    assert "HYBRID conjunction" in summary, summary[:1500]
    assert "both must pass" in summary, summary[:1500]

    rows = _report(out, "markerfinder.hgt_evaluation.tsv")
    header = rows[0].split("\t")
    notes_col = header.index("notes")
    graded = [row.split("\t")[notes_col] for row in rows[1:]]
    assert any("hybrid:" in note for note in graded), graded[:5]
    # Every graded marker names its deciding leg; a plain risk verdict under a
    # Hybrid label is exactly the mislabelling this mode must not do.
    assert all(
        ("hybrid:" in note) or ("UNKNOWN" in note) or (note == "")
        for note in graded
    ), [n for n in graded if n][:8]

    summary_lines = _report(out, "markerfinder.marker_summary.tsv")
    grade_col = summary_lines[0].split("\t").index("consistency_grade")
    grades = {row.split("\t")[grade_col] for row in summary_lines[1:]}
    assert grades <= {"consistent", "inconsistent", "inconclusive"}, grades
    assert grades != {"NA"}, "hybrid mode ran the screen, grades cannot all be NA"

    assert (out / "Phase5_reports" / "markerfinder.excluded_profile.tsv").exists(), (
        "the markers the consistency leg removed must be profiled"
    )


# ──: independence and informativeness, read off the page ────

def test_independence_measurements_and_the_verdict_they_gate_are_on_the_page(risk_run):
    """ At product level. The first execution found a new instance of the
     disease: the summary asserted a property ("recommendation
    confidence is capped at medium") of a dataset-level verdict whose object
    had zero consumers -- ``recommend_tree`` was never called in the run, so
    the page made a claim about thin air. Now the verdict is computed inside
    the run and rendered, so the cap note and the verdict can be checked
    against each other by a reader."""
    summary = (risk_run / "Phase5_reports" / "markerfinder.pipeline_summary.txt").read_text(
        encoding="utf-8"
    )
    assert "Leg independence:" in summary, summary[:1500]
    assert "Tree recommendation: confidence=" in summary, summary[:1500]
    capped = "capped at medium" in summary
    verdict_high = "Tree recommendation: confidence=high" in summary
    # The invariant that makes both lines trustworthy: they may not disagree.
    assert not (capped and verdict_high), summary[:2000]
    # The mock run shares every gene tree across the two legs (100% cache
    # Reuse), so the cap must actually be in force on the verdict.
    assert capped and "confidence=medium" in summary, summary[:2000]


def test_marker_summary_pis_columns_carry_measured_values(risk_run):
    """ At product level: the informativeness columns hold measured
    integers, not NA and not placeholders. The mock alignment is invariant, so
    ``pis=0`` is the TRUE measured answer here -- a value, exactly the case a
    placeholder-to-NA conversion must not destroy."""
    lines = _report(risk_run, "markerfinder.marker_summary.tsv")
    header = lines[0].split("\t")
    pis_col = header.index("pis")
    eff_col = header.index("effective_columns")
    assert len(lines) > 1, lines
    for row in lines[1:]:
        fields = row.split("\t")
        assert fields[pis_col].isdigit(), fields
        assert fields[eff_col].isdigit(), fields
        assert int(fields[eff_col]) > 0, fields


# ── at product level: an unusable interpreter declares its gaps ─────

def test_unmeasured_declaration_reaches_the_summary(tmp_path, monkeypatch):
    import markerfinder.utils.etree as etree_module
    from markerfinder.exceptions import PhyloToolUnavailable

    def _broken():
        raise PhyloToolUnavailable("simulated: ete3 not importable here")

    monkeypatch.setattr(etree_module, "require_ete3", _broken)
    out = _run(tmp_path)
    summary = (out / "Phase5_reports" / "markerfinder.pipeline_summary.txt").read_text(
        encoding="utf-8"
    )
    assert "METRICS NOT MEASURED THIS RUN" in summary, summary[:1500]
    # Named with a cause, not just a coverage fraction.
    assert "ete3 unusable" in summary or "no readable internal support label" in (
        summary
    ), summary[:1500]


# ──: two identical runs must agree byte for byte ───────────────────

def test_control_two_runs_produce_identical_evidence(tmp_path):
    first = _run(tmp_path / "a")
    second = _run(tmp_path / "b")
    for name in ("markerfinder.marker_summary.tsv", "markerfinder.assertions.tsv",
                 "markerfinder.hgt_evaluation.tsv"):
        left = (first / "Phase5_reports" / name).read_bytes()
        right = (second / "Phase5_reports" / name).read_bytes()
        assert left == right, name


# ──: the two gates whose verdict must be observable ──────────

BASELINE = (Path(__file__).resolve().parents[2]
            / "tests" / "benchmark" / "expected" / "taxonomy_mustpass.yaml")


def test_mustpass_gate_writes_a_verdict_for_every_state(tmp_path):
    """Read off the product: the gate must state what it checked, what
    failed, and what it could not check -- and the summary's verdict must agree
    with the rows underneath it.

    The abort path itself (a planted non-monophyletic relation exiting with 4) is
    covered by ``tests/unit/test_taxonomy_mustpass.py``; the shipped baseline's
    marker ids are still the author's to supply, so this test pins the shape and
    the internal consistency of what a real run writes.
    """
    if not BASELINE.exists():  # Pragma: no cover - guard against a moved fixture
        pytest.skip(f"NOT EXECUTED: shipped baseline missing at {BASELINE}")
    out = _run(tmp_path, mustpass=BASELINE)
    path = out / "Phase5_reports" / "markerfinder.mustpass.tsv"
    assert path.exists(), "--taxonomy-mustpass must write mustpass.tsv"
    lines = path.read_text(encoding="utf-8").strip().split("\n")
    header = lines[0].split("\t")
    assert header == ["kind", "marker_id", "relation", "level", "taxon",
                      "members", "result"], header
    rows = [line.split("\t") for line in lines[1:]]
    kinds = {row[0] for row in rows}
    assert kinds <= {"summary", "violation", "not_checked", "detail"}, kinds
    summary = next(row for row in rows if row[0] == "summary")
    violations = [row for row in rows if row[0] == "violation"]
    # The invariant that makes the file trustworthy: one FAIL row means the
    # Summary cannot say PASS, and a PASS summary cannot hide a FAIL row.
    expected = "FAIL" if violations else "PASS"
    assert summary[-1] == expected, (summary, violations)
    if not violations:
        # A clean verdict still has to say WHY it is clean: either every required
        # Relation was checked and held (detail), or the gate names what it could
        # Not check. Silence is not a verdict.
        assert {"detail", "not_checked"} & kinds, rows


def _declared_db(tmp_path: Path) -> Path:
    """A temp --db-dir whose single entry is pinned by the production hash code."""
    import json

    from markerfinder.utils.db_versioning import DatabaseVersionManager

    db = tmp_path / "db"
    (db / "ar53").mkdir(parents=True)
    (db / "ar53" / "TIGR00001.hmm").write_text("HMMER3/x\n", encoding="utf-8")
    manager = DatabaseVersionManager(str(db))
    correct = manager.verify_database("ar53", str(db / "ar53")).get("hash")
    assert correct, "hash computation returned nothing -- the fixture is wrong"
    (db / "expected_hashes.json").write_text(
        json.dumps({"ar53": correct}), encoding="utf-8"
    )
    return db


def _check(db_dir: Path):
    import subprocess
    import sys

    return subprocess.run(
        [sys.executable, "-m", "markerfinder", "--check", "--db-dir", str(db_dir)],
        cwd=str(Path(__file__).resolve().parents[2]),
        capture_output=True, text=True, timeout=600,
    )


def test_ac12_hash_mismatch_fails_the_process_and_a_match_does_not(tmp_path):
    """ At process level, both directions.

    The unit-level trio (mismatch FAIL / unpopulated INFO / absent NOT CHECKED)
    lives in ``test_db_hash_selfcheck.py``; what is pinned here is the consequence
    a user actually sees -- the exit code and the named line -- plus the control
    that a matching database does not fail.
    """
    good = _declared_db(tmp_path / "ok")
    passing = _check(good)
    assert "Database hash [ar53]" in passing.stdout, passing.stdout[-1200:]
    assert "SHA-256 matches" in passing.stdout, passing.stdout[-1200:]

    tampered = _declared_db(tmp_path / "bad")
    (tampered / "ar53" / "TIGR00001.hmm").write_text("HMMER3/tampered\n", encoding="utf-8")
    failing = _check(tampered)
    assert failing.returncode != 0, failing.stdout[-1200:]
    assert "[FAIL]" in failing.stdout, failing.stdout[-1200:]
    assert "Database hash [ar53]" in failing.stdout, failing.stdout[-1200:]
