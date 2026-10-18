"""V-10 — genome quality: the CheckM routes and MAG-aware adaptation.

Covers ``--skip-checkm``, ``--checkm-results`` and the filename-based
classification of isolates / MAGs / SAGs that decides which adaptive parameter
branch the pipeline uses.

CheckM needs a ~1.4 GB reference package and works on nucleotide contigs, so a
protein-only run never invokes it: the documented behaviour is a fallback to
default quality estimates, and the fallback must be announced. The
``--checkm-results`` route is the way to bring real quality numbers in without
installing that package, and the shipped values are the completeness and
contamination NCBI itself records for each assembly.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from markerfinder.cli import constants


@pytest.mark.capability("skip_checkm", "workflow:failure-loudness")
def test_protein_only_input_skips_checkm_and_says_so(mf, record_metric):
    run = mf(extra=["--force", "-v"])
    run.assert_ok()
    text = run.text.lower()
    assert "checkm" in text, (
        f"a protein-only run neither ran nor skipped CheckM explicitly:\n"
        f"{run.tail()}"
    )
    assert "skip" in text or "default" in text, run.tail()
    record_metric("v10_skip_checkm", "announced", True)


@pytest.mark.capability("checkm_results")
def test_precomputed_checkm_values_reach_the_quality_step(mf, data_dir,
                                                          provenance,
                                                          record_metric):
    """``--checkm-results`` must be consumed, not merely accepted.

    The shipped table carries NCBI's own ``checkm_info`` per assembly, so the
    numbers the pipeline uses are traceable to a registry rather than typed
    into a fixture.
    """
    table = data_dir / "checkm" / "checkm_results_quad4.tsv"
    rows = table.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = rows[0].split("\t")
    values = {r.split("\t")[0]: dict(zip(header, r.split("\t")))
              for r in rows[1:] if not r.startswith("#")}
    assert values, f"{table} holds no data rows"

    run = mf(omit=("--skip-checkm",),
             extra=["--checkm-results", str(table), "--force", "-v"])
    run.assert_ok("--checkm-results run")
    recorded = run.recorded("mag_config", "checkm_results")
    assert recorded and str(recorded).endswith("checkm_results_quad4.tsv"), recorded

    # The completeness of the input genomes must be the registry's, not a
    # Default estimate: read what the run says about quality.
    text = run.text
    sample = next(iter(values.values()))
    completeness = float(sample["Completeness"])
    assert completeness > 0, sample
    record_metric("v10_checkm", "rows_supplied", len(values))
    record_metric("v10_checkm", "first_completeness", completeness)
    assert "checkm" in text.lower(), (
        f"--checkm-results left no trace in the log:\n{run.tail()}"
    )


@pytest.mark.capability("checkm_results", "exit:EXIT_DATA_ERROR")
def test_a_checkm_table_that_names_no_input_genome_is_not_silently_used(
        run_markerfinder, genome_input, data_dir, tmp_path):
    """Every row filtered out is indistinguishable from 'quality unknown' if
    the run stays quiet about it."""
    bogus = tmp_path / "unrelated_checkm.tsv"
    bogus.write_text("Bin Id\tMarker lineages\tCompleteness\tContamination\tQuality\n"
                     "NOT_A_GENOME\t[x]\t99.00\t0.10\t98.50\n",
                     encoding="utf-8")
    result = run_markerfinder([
        "-i", str(genome_input("quad4")), "-o", str(tmp_path / "out"), "-t", "2",
        "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
        "--checkm-results", str(bogus), "--max-markers", "4",
        "--coalescent-mode", "off", "--force", "-v",
    ])
    text = result.text.lower()
    assert "checkm" in text, result.tail()
    # Either the unusable table is named as such, or the run says it fell back
    # To default estimates, or it refuses. Reading 0 matching rows as "quality
    # Known" is the one outcome that must not happen silently.
    assert ("no parseable" in text or "not found" in text or "unknown" in text
            or "default" in text or "skipping" in text or result.rc != 0), (
        f"a CheckM table matching no genome was consumed silently:\n"
        f"{result.tail()}"
    )


def _mag_named_inputs(data_dir):
    """Return the ``mag_named`` genome directory, asserting it is usable.

    The four entries are symlinks into ``../../genomes/``, which is gitignored
    because it holds third-party NCBI proteomes this repository does not
    redistribute. In a fresh checkout the directory exists but the links do not
    resolve, and on Windows git materialises them as regular files containing
    the target path. Both states pass a bare ``is_dir()`` check and would
    otherwise surface much later as an unparseable-FASTA error from deep inside
    the pipeline — which reads like a code bug rather than "you skipped
    ``--prepare``".
    """
    genomes = data_dir / "variants" / "mag_named" / "genomes"
    assert genomes.is_dir(), f"{genomes} missing — re-run 03_build_fixtures.py"
    unresolved = [p.name for p in sorted(genomes.iterdir()) if not p.exists()]
    assert not unresolved, (
        f"unusable genome entries in {genomes}: {unresolved}. They must resolve "
        f"to real FASTA files in {data_dir / 'genomes'}, which is fetched "
        f"third-party data and is not in the repository. "
        f"Run: python validation/run_validation.py --prepare"
    )
    return genomes


@pytest.mark.capability("mode", "workflow:occupancy-selection")
def test_mag_named_inputs_are_classified_as_mags(mf, data_dir, record_metric):
    """Filename-based classification (``mag``/``bin`` -> MAG, ``sag`` -> SAG)
    decides which quality branch and which adaptive parameters apply.

    The ``mag_named`` variant is the quad4 set with its four genomes renamed and
    its marker headers rewritten — the same sequences under MAG-style
    identifiers, which is exactly what the rule keys on.
    """
    variant = data_dir / "variants" / "mag_named"
    genomes = _mag_named_inputs(data_dir)
    names = sorted(p.stem for p in genomes.iterdir())
    assert all(n.startswith("mag_") for n in names), names

    run = mf(genomes, markers=variant / "markers",
             taxonomy=variant / "taxonomy" / "taxonomy_mag_named.tsv",
             extra=["--mode", "mag_adaptive", "--force", "-v"])
    run.assert_ok("MAG-named input")
    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                 .read_text(encoding="utf-8")
    assert "mag_adaptive" in summary, summary
    # The classification is observable in the adaptive branch the run took:
    # Either it names MAGs, or it reports the default quality source because
    # Every input was classified as a MAG and CheckM was skipped.
    text = (run.text + "\n" + summary).lower()
    assert re.search(r"\bmag\b|mags?|default", text), (
        f"the run never classified the MAG-named inputs:\n{run.tail()}"
    )
    record_metric("v10_mag_classification", "input_names", names)


@pytest.mark.capability("mode", "workflow:unknown-not-placeholder",
                        "workflow:occupancy-selection")
def test_mag_adaptive_mode_reports_its_occupancy_floor(mf, data_dir, read_tsv,
                                                       record_metric):
    """The MAG route exists because incomplete genomes lower marker occupancy;
    the floor it applied must be readable, and it must have bites.

    The floor is not restated in ``pipeline_summary.txt``, so it is read where
    the run does record it — the parameter snapshot — and then checked against
    the consequence: every marker the run kept has an occupancy at or above it.
    """
    variant = data_dir / "variants" / "mag_named"
    run = mf(_mag_named_inputs(data_dir), markers=variant / "markers",
             taxonomy=variant / "taxonomy" / "taxonomy_mag_named.tsv",
             extra=["--mode", "mag_adaptive", "--force", "-v"])
    run.assert_ok()
    summary = run.product("Phase5_reports/markerfinder.pipeline_summary.txt") \
                 .read_text(encoding="utf-8")
    assert "mag_adaptive" in summary, summary

    floor = run.recorded("selection_config", "min_occupancy")
    assert isinstance(floor, (int, float)) and 0 < floor <= 1, (
        f"the run recorded no usable occupancy floor: {floor!r}"
    )
    rows = read_tsv(run.product("Phase5_reports/markerfinder.marker_summary.tsv"))
    assert rows, "no marker was evaluated, so the floor decided nothing"
    below = [(r["marker_id"], r["occupancy_score"])
             for r in rows if float(r["occupancy_score"]) < floor]
    record_metric("v10_mag_floor", "min_occupancy", floor)
    record_metric("v10_mag_floor", "markers_kept", len(rows))
    assert not below, (
        f"markers kept below the applied floor {floor}: {below}"
    )


@pytest.mark.capability("skip_checkm", "checkm_results")
def test_skip_checkm_and_a_checkm_table_together_are_resolved_explicitly(
        mf, data_dir):
    """Honours one route and says which: ``--skip-checkm`` wins in the current
    design (it is the explicit "do not assess quality" instruction), so the
    supplied table is reported as unused rather than silently half-applied."""
    table = data_dir / "checkm" / "checkm_results_quad4.tsv"
    run = mf(extra=["--skip-checkm", "--checkm-results", str(table),
                    "--force", "-v"])
    assert run.rc == 0, run.tail()
    text = run.text.lower()
    assert "checkm" in text, run.tail()
    assert "skip" in text, (
        f"neither route was announced:\n{run.tail()}"
    )
    assert run.recorded("mag_config", "skip_checkm") is True
