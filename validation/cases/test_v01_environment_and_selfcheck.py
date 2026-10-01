"""V-01 — environment, self-check, terminal states and data integrity.

Covers: ``--check`` / ``--self-test``, ``--version``, the external-tool surface
the pipeline really shells out to, the exit-code alphabet the documentation
promises (failure classes must be distinguishable), and the two
integrity gates over the shipped test data itself.
"""

from __future__ import annotations

import hashlib
import os
import re

import pytest

from markerfinder.cli import constants
from markerfinder.utils import dependency_check


@pytest.mark.capability("self_test", "version", "db_dir", "exit:EXIT_SUCCESS")
def test_check_reports_every_item_and_exits_zero(run_markerfinder):
    result = run_markerfinder(["--check"])
    text = result.text
    assert result.rc == 0, (
        f"--check must exit 0 on a complete install, got exit code {result.rc}"
        f"\n{text}"
    )
    # A green wall of text is not a self-test: the item classes must appear.
    for probe in ("python", "ete3", "external", "database", "assertion"):
        assert probe in text.lower(), f"--check never mentions {probe!r}:\n{text}"
    fails = [line for line in text.splitlines() if re.search(r"\bFAIL\b", line)]
    assert not fails, f"--check reported failures:\n" + "\n".join(fails)


@pytest.mark.capability("self_test", "db_dir")
def test_check_with_an_explicit_db_dir_checks_that_directory(run_markerfinder,
                                                             tmp_path):
    """'s second clause: ``--check --db-dir X`` must examine X, not the
    repository default. A check pointed at the wrong artifact still prints a
    verdict — about something the user never named."""
    other = tmp_path / "emptydb"
    other.mkdir()
    result = run_markerfinder(["--check", "--db-dir", str(other)])
    text = result.text.lower()
    assert other.name in result.text or "not present" in text or "no populated" in text, (
        f"--check did not look at the db dir it was handed:\n{result.text}"
    )


@pytest.mark.capability("version", "exit:EXIT_SUCCESS")
def test_version_string_is_the_installed_version(run_markerfinder):
    import markerfinder

    result = run_markerfinder(["--version"])
    assert result.rc == 0
    assert markerfinder.__version__.split("+")[0] in result.text, (
        f"--version says {result.text!r} but the package imports as "
        f"{markerfinder.__version__!r}"
    )


@pytest.mark.capability("workflow:provenance-recorded")
def test_banner_claims_only_declared_runtime_dependencies(run_markerfinder):
    """The startup screen attributes licences; it may only name libraries the
    distribution actually requires.

    The banner used to print pandas / numpy / scipy / rich on every run. None
    of them is imported anywhere in the package and none is a declared
    dependency since so the attribution was false on its face.
    """
    from markerfinder import banner

    declared = set(banner.runtime_dependencies())
    assert declared, "runtime dependency discovery returned nothing"
    for forbidden in ("pandas", "numpy", "scipy", "rich"):
        assert forbidden not in {d.replace("-", "") for d in declared}, forbidden

    result = run_markerfinder(["--version"])  # Does not print the banner
    run = run_markerfinder(["--help"])
    printed = {name.lower() for name, _lic in banner.third_party_licenses()}
    assert printed == {d.lower() for d in declared}, (
        f"banner list {sorted(printed)} != declared {sorted(declared)}"
    )
    assert run.rc == 0


@pytest.mark.capability("workflow:failure-loudness")
def test_every_external_tool_the_code_calls_is_found_or_reported():
    """``EXTERNAL_TOOLS`` is derived from the subprocess call sites, so this is
    a statement about the installation, not about the code."""
    status = dependency_check.check_external_tools()
    assert status, "no tools were checked at all"
    missing = sorted(t for t, ok in status.items() if not ok)
    critical = {name for name, _p, is_critical in dependency_check.EXTERNAL_TOOLS
                if is_critical}
    assert not (set(missing) & critical), (
        f"critical tools missing: {sorted(set(missing) & critical)}"
    )
    if missing:
        print(f"optional tools not installed on this machine: {missing}")


@pytest.mark.capability("workflow:failure-loudness")
def test_tool_to_step_mapping_covers_every_tool():
    status = dependency_check.check_external_tools()
    unmapped = sorted(set(status) - set(dependency_check.TOOL_STEPS))
    assert not unmapped, f"tools with no recorded call site: {unmapped}"
    bad_steps = sorted(
        {step for steps in dependency_check.TOOL_STEPS.values() for step in steps}
        - set(dependency_check.PIPELINE_STEPS)
    )
    assert not bad_steps, f"TOOL_STEPS names unknown steps: {bad_steps}"


@pytest.mark.capability("exit:EXIT_SUCCESS", "exit:EXIT_ARG_ERROR",
                        "exit:EXIT_DATA_ERROR", "exit:EXIT_ASSERTION_FAILED",
                        "exit:EXIT_INCONCLUSIVE", "exit:EXIT_INTERRUPT",
                        "exit:EXIT_RUNTIME_ERROR", "workflow:failure-loudness")
def test_exit_codes_are_pairwise_distinct():
    codes = {n: v for n, v in vars(constants).items()
             if n.startswith("EXIT_") and isinstance(v, int)}
    assert len(set(codes.values())) == len(codes), codes
    assert codes["EXIT_SUCCESS"] == 0
    assert codes["EXIT_ARG_ERROR"] == 2
    assert codes["EXIT_DATA_ERROR"] == 3
    assert codes["EXIT_ASSERTION_FAILED"] == 4
    assert codes["EXIT_INCONCLUSIVE"] == 5


@pytest.mark.capability("exit:EXIT_ARG_ERROR", "mode")
def test_bad_argument_value_exits_with_the_argument_code(run_markerfinder,
                                                         tmp_path):
    """An unusable option value is an argument error (2), not a runtime error
    (1) and not a traceback."""
    result = run_markerfinder(["-i", str(tmp_path), "-o", str(tmp_path / "o"),
                               "--mode", "not_a_mode"])
    assert result.rc == constants.EXIT_ARG_ERROR, result.tail()
    assert "invalid choice" in result.text.lower(), result.tail()


@pytest.mark.capability("input", "exit:EXIT_ARG_ERROR")
def test_missing_input_directory_is_refused_before_any_work(run_markerfinder,
                                                           tmp_path):
    result = run_markerfinder(["-i", str(tmp_path / "nope"),
                               "-o", str(tmp_path / "out")])
    assert result.rc == constants.EXIT_ARG_ERROR, result.tail()
    assert "does not exist" in result.text.lower(), result.tail()
    assert not (tmp_path / "out").exists(), (
        "an output directory was created by a run that never started")


@pytest.mark.capability("data:provenance-integrity")
def test_shipped_genomes_are_the_organisms_the_provenance_claims(data_dir,
                                                                provenance):
    """Data-integrity gate over the bundle itself.

    A hand-typed lineage once labelled *Salmonella enterica* GCF_000008105.1 as
    ``o__Bacillales;f__Bacillaceae;g__Bacillus``. The taxonomy table drives the
    monophyly screen, so that is an invalid experiment, not a typo. NCBI writes
    the organism into every protein header, so the claim is checked against the
    sequence file instead of against whoever typed it.

    NCBI labels a MULTISPECIES protein only down to the rank the name is shared
    at (``[Bacteria]``, ``[Bacillus]``), so the first header alone proves
    nothing: the genus must be the majority organism among a sample of headers.
    """
    genomes = data_dir / "genomes"
    sample = 500
    for row in provenance:
        faa = genomes / f"{row['accession']}.faa"
        assert faa.exists(), f"{faa} missing — re-run 01_fetch_genomes.py"
        genus = row["genus"].split()[0]
        headers = []
        with open(faa, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    headers.append(line)
                    if len(headers) >= sample:
                        break
        assert headers, f"{faa} has no sequence records"
        hits = sum(1 for h in headers if genus in h)
        fraction = hits / len(headers)
        assert fraction >= 0.5, (
            f"{row['accession']}: only {hits}/{len(headers)} of the first "
            f"sampled protein headers mention genus {genus!r}, so the shipped "
            f"proteome does not look like the provenanced organism "
            f"({row['organism_name']!r})"
        )


@pytest.mark.capability("data:provenance-integrity")
def test_provenance_rows_are_complete(provenance):
    for row in provenance:
        for rank in ("domain", "phylum", "class", "order", "family", "genus"):
            assert row[rank], f"{row['accession']}: empty {rank} in PROVENANCE.tsv"
        assert row["domain"] in {"Bacteria", "Archaea"}, row["domain"]
        assert row["tax_id"].isdigit(), row["tax_id"]
        assert row["sha256"] and len(row["sha256"]) == 64, row["accession"]
        assert row["source_registry"] == "NCBI RefSeq", row["source_registry"]


@pytest.mark.capability("data:manifest-checksums")
def test_manifest_checksums_match_the_shipped_files(data_dir):
    """The bundle ships with the release, so every byte of it must verify.

    One mismatch means the published data is not the data the cases ran
    against, which would make every number in the report unverifiable.

    ``genomes/`` and ``hmms/`` are the parts of the manifest that are not
    committed: the NCBI RefSeq proteomes are re-downloaded by
    ``01_fetch_genomes.py`` and the working copy of the HMM profiles is
    materialised by ``02_build_marker_sets.py`` out of ``db/gtdb_markers``. The
    manifest is what tells the builder (and this check) that it got the same
    bytes the report was measured on. Any other absent file is a real defect.

    The four ``variants/mag_named/genomes/mag_0*.faa`` entries are the same
    category of uncommitted data reached through a symlink, so the prefix test
    has to look at where a link *points*, not only at where it sits. They are
    committed as links precisely because the proteomes are not redistributed,
    and ``03_build_fixtures.py`` recreates them. Judging them by their own path
    made this check pass only on a machine where the fetch had already run, and
    fail everywhere else — the manifest is a record of what a *prepared*
    environment contains, so an un-prepared one must not report corruption.
    """
    manifest = data_dir / "MANIFEST.sha256"
    assert manifest.exists(), f"{manifest} missing — re-run 03_build_fixtures.py"

    def _is_uncommitted(rel: str) -> bool:
        """True when this entry's bytes come from fetched, not committed, data."""
        if rel.startswith(("genomes/", "hmms/")):
            return True
        # A committed symlink into one of those trees: resolve it and classify
        # by the target. Path.is_symlink() survives the target being absent,
        # which is exactly the state that used to be reported as corruption.
        path = data_dir / rel
        if path.is_symlink():
            target = (path.parent / os.readlink(path)).resolve()
            data_root = data_dir.resolve()
            try:
                inside = target.relative_to(data_root).as_posix()
            except ValueError:
                return False
            return inside.startswith(("genomes/", "hmms/"))
        return False

    checked = 0
    fetched = 0
    mismatches = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, _, rel = line.partition("  ")
        rel = rel.strip()
        if rel == "MANIFEST.sha256":
            continue
        path = data_dir / rel
        if not path.exists():
            if _is_uncommitted(rel):
                fetched += 1
                continue
            mismatches.append(f"{rel} (absent)")
            continue
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            mismatches.append(rel)
        checked += 1
    assert not mismatches, f"checksum problems: {mismatches[:10]}"
    if fetched:
        print(f"{fetched} fetched data file(s) not present — run 01_fetch_genomes.py "
              f"and 02_build_marker_sets.py before the slow cases")
    assert checked > 50, f"only {checked} files covered by MANIFEST.sha256"
