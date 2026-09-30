"""V-16 — robustness: degenerate and malformed inputs must fail in the right way.

Three properties are checked over the failure paths, none of which a
happy-path test can see:

1. nothing crashes with a bare traceback (failures must be
   distinguishable and actionable);
2. nothing "succeeds" by producing an empty or placeholder result
   (the silent-zero class of defect);
3. the supported range is respected — a dataset below the four-tip minimum is
   refused, not forced through.

Also exercised here: the shipped marker files that hold fewer tips than the
pipeline requires, and non-FASTA junk inside the input directory.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from markerfinder.cli import constants

EXTRA_ARGS = ["-t", "2", "--max-markers", "4", "--coalescent-mode", "off",
              "--skip-checkm", "--force"]


def _run(run_markerfinder, data_dir, genome_input, name, out_dir, extra=()):
    return run_markerfinder([
        "-i", str(genome_input(name)), "-o", str(out_dir),
        "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / f"{name}_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / f"taxonomy_{name}.tsv"),
    ] + EXTRA_ARGS + list(extra))


def _assert_loud_failure(result, out_dir: Path):
    assert result.rc != 0, (
        f"expected a refusal, got exit 0:\n{result.tail()}"
    )
    assert "Traceback" not in result.text, (
        f"the failure was an unhandled exception:\n{result.tail()}"
    )
    assert result.text.strip(), "the failing run produced no message at all"


@pytest.mark.capability("input", "exit:EXIT_DATA_ERROR",
                        "workflow:failure-loudness")
def test_empty_input_directory_is_refused(run_markerfinder, empty_input_dir,
                                         tmp_path, data_dir):
    empty = empty_input_dir("v16_empty")
    result = run_markerfinder([
        "-i", str(empty), "-o", str(tmp_path / "out"), "-t", "2",
        "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
        "--skip-checkm", "--force",
    ])
    _assert_loud_failure(result, tmp_path / "out")
    assert "no " in result.text.lower() or "empty" in result.text.lower(), (
        result.tail())


@pytest.mark.capability("input", "workflow:failure-loudness")
def test_a_single_genome_cannot_be_phyligenetically_placed(run_markerfinder,
                                                          data_dir, tmp_path,
                                                          genome_input):
    result = _run(run_markerfinder, data_dir, genome_input, "single1",
                  tmp_path / "single")
    _assert_loud_failure(result, tmp_path / "single")


@pytest.mark.capability("input", "workflow:failure-loudness")
def test_three_tips_are_below_the_supported_minimum(run_markerfinder, data_dir,
                                                   tmp_path, genome_input,
                                                   record_metric):
    """The pipeline refuses to build a gene tree under four tips; the refusal
    must be a data error, not a zero-marker "success"."""
    result = _run(run_markerfinder, data_dir, genome_input, "pair3",
                  tmp_path / "three")
    record_metric("v16_three_tips", "exit_code", result.rc)
    _assert_loud_failure(result, tmp_path / "three")
    text = result.text.lower()
    assert "tree" in text or "marker" in text, result.tail()
    trees = tmp_path / "three" / "Phase4_trees"
    if trees.exists():
        empties = [p.name for p in trees.glob("*.newick") if p.stat().st_size == 0]
        assert not empties, f"empty tree files published as products: {empties}"


@pytest.mark.capability("gtdb_markers_dir", "workflow:unknown-not-placeholder")
def test_marker_files_below_the_informative_minimum_are_reported(
        run_markerfinder, data_dir, genome_input, tmp_path, record_metric):
    """``markers/degenerate`` holds a 2-tip and a 1-tip marker carved out of a
    real one (see its README). The pipeline must not build a tree from them and
    must not pretend the dataset was usable."""
    degenerate = data_dir / "markers" / "degenerate"
    assert degenerate.is_dir(), f"{degenerate} missing — re-run 02_build_marker_sets.py"
    result = run_markerfinder([
        "-i", str(genome_input("quad4")), "-o", str(tmp_path / "degen"),
        "-t", "2", "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(degenerate),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
        "--max-markers", "4", "--coalescent-mode", "always", "--skip-checkm",
        "--force",
    ])
    record_metric("v16_degenerate", "exit_code", result.rc)
    if result.rc == 0:
        # Allowed only if the run said the markers were unusable rather than
        # Publishing a tree built from one or two tips.
        text = result.text.lower()
        assert "skip" in text or "not" in text or "warn" in text, (
            f"sub-minimal marker files produced a quiet success:\n"
            f"{result.tail()}"
        )
        tree = tmp_path / "degen" / "Phase4_trees" / "markerfinder.species_tree_concat.newick"
        if tree.exists():
            from markerfinder.utils.tree_utils import parse_newick_tips
            tips = set(parse_newick_tips(tree.read_text(encoding="utf-8")))
            assert len(tips) >= 4, (
                f"a species tree with {len(tips)} tips was published from "
                "sub-minimal markers"
            )
    else:
        _assert_loud_failure(result, tmp_path / "degen")


@pytest.mark.capability("input", "workflow:failure-loudness")
def test_non_fasta_files_in_the_input_directory_are_not_silently_counted(
        run_markerfinder, data_dir, genome_input, tmp_path, case_workdir):
    junk = case_workdir / "input"
    junk.mkdir()
    for acc in (data_dir / "sets" / "quad4.txt").read_text().split():
        src = data_dir / "genomes" / f"{acc.strip()}.faa"
        if acc.strip() and src.exists():
            shutil.copy2(src, junk / src.name)
    (junk / "README.md").write_text("# not a genome\n", encoding="utf-8")
    (junk / "notes.txt").write_text("hello\n", encoding="utf-8")
    (junk / "empty.faa").write_text("", encoding="utf-8")

    result = run_markerfinder([
        "-i", str(junk), "-o", str(tmp_path / "junk"),
        "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
    ] + EXTRA_ARGS + ["--force"])
    text = result.text.lower()
    # Junk must be *named*, and an empty FASTA must not be counted as a genome:
    # A silently skipped input is how "the run used 4 genomes" and "the
    # Directory held 5.faa files" come to disagree.
    assert "not read" in text or "skipping" in text or "skipped" in text, (
        f"non-FASTA and empty inputs left no trace:\n{result.tail()}"
    )
    for junk in ("readme.md", "notes.txt", "empty.faa"):
        assert junk in text, f"{junk} was not named in the report:\n{result.tail()}"
    if result.rc == 0:
        summary = (tmp_path / "junk" / "Phase5_reports" /
                   "markerfinder.pipeline_summary.txt")
        if summary.exists():
            body = summary.read_text(encoding="utf-8")
            genomes = [ln for ln in body.splitlines() if "Genomes:" in ln]
            assert genomes, body[:400]
            count = int(genomes[0].split(":")[1].strip())
            assert count == 4, (
                f"the run counted {count} genomes; the directory holds 4 non-"
                "empty FASTA files (the empty one must be dropped and named)"
            )
    else:
        _assert_loud_failure(result, tmp_path / "junk")


@pytest.mark.capability("input", "exit:EXIT_ARG_ERROR")
def test_a_file_passed_as_the_input_directory_is_refused(run_markerfinder,
                                                        data_dir, tmp_path):
    one_faa = data_dir / "genomes" / "GCF_000009045.1.faa"
    result = run_markerfinder(["-i", str(one_faa), "-o", str(tmp_path / "o"),
                               "-t", "2"])
    assert result.rc == constants.EXIT_ARG_ERROR, result.tail()
    assert "directory" in result.text.lower(), result.tail()


@pytest.mark.capability("config", "exit:EXIT_DATA_ERROR")
def test_an_unparseable_config_file_is_a_named_error(run_markerfinder,
                                                    data_dir, tmp_path,
                                                    case_workdir):
    broken = case_workdir / "broken.yaml"
    broken.write_text("input: [\n  unclosed\n", encoding="utf-8")
    result = run_markerfinder([
        "-i", str(data_dir / "genomes"), "-o", str(tmp_path / "o"),
        "--config", str(broken),
    ])
    assert result.rc != 0, result.tail()
    assert "Traceback" not in result.text, result.tail()
    assert "config" in result.text.lower(), result.tail()


@pytest.mark.capability("marker_hmm_dir", "workflow:failure-loudness")
def test_a_marker_directory_of_unusable_profiles_is_reported(run_markerfinder,
                                                            genome_input,
                                                            tmp_path,
                                                            case_workdir,
                                                            data_dir):
    fake = case_workdir / "bad_hmms"
    fake.mkdir()
    (fake / "NOTAPROFILE.HMM").write_text("this is not an HMM profile\n",
                                          encoding="utf-8")
    result = run_markerfinder([
        "-i", str(genome_input("quad4")), "-o", str(tmp_path / "badhmm"),
        "-t", "2", "--marker-mode", "hmm", "--marker-hmm-dir", str(fake),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
        "--skip-checkm", "--force",
    ])
    if result.rc == 0:
        assert "warn" in result.text.lower() or "no marker" in result.text.lower(), (
            f"an unreadable HMM directory produced a silent success:\n"
            f"{result.tail()}"
        )
    else:
        _assert_loud_failure(result, tmp_path / "badhmm")


@pytest.mark.capability("output", "workflow:failure-loudness")
def test_an_unwritable_output_path_is_reported_not_traced(run_markerfinder,
                                                         genome_input,
                                                         data_dir, tmp_path):
    locked = tmp_path / "locked"
    locked.mkdir()
    (locked / "blocks").write_text("x", encoding="utf-8")
    # Ask for the output to be a path occupied by a FILE: mkdir must fail and
    # The pipeline has to turn that into its own error.
    result = run_markerfinder([
        "-i", str(genome_input("quad4")), "-o", str(locked / "blocks"),
        "-t", "2", "--marker-mode", "gtdb_tk",
        "--gtdb-markers-dir", str(data_dir / "markers" / "quad4_core"),
        "--taxonomy-table", str(data_dir / "taxonomy" / "taxonomy_quad4.tsv"),
        "--skip-checkm", "--force",
    ])
    assert result.rc != 0, result.tail()
    assert "Traceback" not in result.text, result.tail()
    assert "exist" in result.text.lower() or "director" in result.text.lower() \
        or "output" in result.text.lower(), result.tail()
