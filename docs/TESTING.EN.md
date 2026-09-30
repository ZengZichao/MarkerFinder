# MarkerFinder Test Plan and Validation Methodology

Applies to: MarkerFinder 0.1.0
[中文版本](TESTING.CN.md) · Companion report: [TEST-REPORT.EN.md](TEST-REPORT.EN.md)

---

## 1. What is being tested, and why in this shape

MarkerFinder's published claims are of three kinds, and each kind needs a
different instrument:

| Claim kind | Example | Instrument |
|---|---|---|
| Value-level | `rf_distance` of two topologies is in 0–1 | unit tests, `tests/unit/` |
| Structure-level | `run_config.json` records the database version | integration tests with mocked tools, `tests/integration/` |
| Consequence-level | "HGT screening removes phylogenetically misleading markers" | **real end-to-end runs against a positive control**, `validation/` |

A consequence-level claim cannot be established by mocking: if MAFFT is mocked,
the alignment is whatever the mock returns. The retracted *Science* 2009 case
that motivated this project's self-check layer failed on exactly this gap —
every internal check passed, and
the topology was biologically absurd. So the suite is split:

```
tests/       988 fast cases, no external tools required, no data download
validation/  153 end-to-end cases on real genomes with real tools
```

Both layers ship with the release, and `validation/` additionally ships its
**data**, its **provenance**, and the **archived results** of the published run.

## 2. The capability surface: how "everything supported" is defined

"Cover every advertised feature and the whole workflow" is only checkable if the
claim is enumerated from a source of truth. That source is the code:

* `markerfinder.cli.parser._build_parser` → 68 CLI options (including four
  documented deprecated aliases and one documented no-op);
* `markerfinder.cli.parser._SUBCOMMAND_SPECS` → the four resumable steps;
* `markerfinder.cli.constants` → seven distinct exit codes;
* the product inventory in `validation/capabilities.py` → 19 output artefacts
  (`Phase5_reports/*`, `Phase5_evidence/*`, `Phase4_trees/*`,
  `Phase4_alignments/*`, `Phase4_intermediate/*`, `.markerfinder/*`);
* nine workflow/data properties: occupancy-based selection, HGT grading, state
  persistence, determinism, provenance recorded, unknown-not-placeholder,
  failure loudness, data provenance integrity, manifest checksums. (Stepwise vs
  aggregate equivalence is asserted by cases marked under the options it is
  driven by, not as a property of its own.)

Coverage is not a hand-maintained list. Each case declares what it covers:

```python
@pytest.mark.capability("hgt_scan", "product:Phase5_reports.threshold_scan.tsv")
def test_threshold_scan_writes_its_own_artifact(...):
```

`conftest.pytest_collection_modifyitems` dumps those markers to
`validation/.work/capability_snapshot/*.json` at collection time, and
`test_v90_capability_coverage.py` enforces two directions:

1. a surface entry with no case → **fail**;
2. a case marker with no surface entry → **fail** (an option was removed while
   its test marker stayed behind).

The matrix is exported to `validation/results/capability_matrix.tsv` with the
covering case ids per capability, so the report's coverage claim is the same
data the test enforces. A partial collection (e.g. `pytest -k something`)
triggers an explicit skip rather than a weakened pass.

Option *roles* are recorded too: a deprecated alias is covered by asserting that
it still works **and** announces itself; `--no-hgt-adaptive-thresholds` is
covered as the documented no-op, i.e. by asserting it behaves like the default.

## 3. Test data

### 3.1 Provenance

All sequence data is public and re-derivable; nothing was typed by hand.

| Item | Source | Recorded where |
|---|---|---|
| 12 proteomes | NCBI RefSeq, `datasets download genome accession … --include protein` — fetched by `01_fetch_genomes.py`; `genomes/` is git-ignored and verified against `MANIFEST.sha256` | `validation/data/genomes/` + `PROVENANCE.tsv` |
| Organism + 7 lineage ranks | NCBI `datasets summary genome accession` → `datasets summary taxonomy taxon <tax_id>` | `PROVENANCE.tsv` columns `domain…species`, `tax_id` |
| Assembly stats (contigs, length, GC, protein count) | same NCBI report | `PROVENANCE.tsv` |
| Completeness / contamination | NCBI's own `checkm_info` per assembly | `PROVENANCE.tsv`, emitted as `data/checkm/*.tsv` |
| Checksums | SHA-256 per file | `PROVENANCE.tsv` + `data/MANIFEST.sha256` |
| Marker HMM profiles | `db/gtdb_markers/bac120`, assembled locally by `scripts/fetch_marker_db.py` (third-party models, not committed) | copied, not modified, into `data/hmms/` — also a local copy, git-ignored |

Three build scripts reproduce the whole bundle from NCBI (`01_fetch_genomes.py`,
`02_build_marker_sets.py`, `03_build_fixtures.py`); they are idempotent and
`--force`-able, and `validation/data/.hit_cache/` holds the HMM scan cache they
derive from.

### 3.2 An error the data layer caught

The earlier revision of this dataset hand-typed its taxonomy table and labelled
`GCF_000008105.1` as `o__Bacillales; f__Bacillaceae; g__Bacillus`. NCBI reports
that assembly as **Salmonella enterica** serovar Choleraesuis str. SC-B67 —
phylum Pseudomonadota, order Enterobacterales. Because the taxonomy table drives
the MAD-rooting + monophyly screen, a wrong rank is not a typo but an invalid
experiment: it decides which clade is expected to be monophyletic.

Two counters exist now:

* lineages are fetched from NCBI by `01_fetch_genomes.py`, and a table without a
  `tax_id` cannot be produced;
* `test_shipped_genomes_are_the_organisms_the_provenance_claims` reads each
  shipped proteome's own FASTA headers (NCBI writes `[organism]` per protein)
  and requires the recorded genus to be the majority organism.

### 3.3 Marker inputs are a layout reconstruction, and say so

`--marker-mode gtdb_tk` consumes GTDB-Tk's per-marker unaligned FASTA. GTDB-Tk
needs a ~100 GB reference-tree package that cannot ship in a test bundle, so the
shipped marker files are **a reconstruction of that layout**: for each marker
profile in `db/gtdb_markers/bac120`, each genome's best hmmsearch hit at
bitscore ≥ 20, one file per marker, filename stem = marker id, tips named by
assembly accession. What is preserved is exactly what the pipeline reads (file
layout, marker ids, unaligned protein sequences, per-marker occupancy). What is
**not** reproduced is GTDB-Tk's taxonomy-aware best-hit choice, so these files
must not be read as a GTDB-Tk release output. `OCCUPANCY_ALL.tsv` records the
occupancy of all 121 scanned profiles; each marker directory carries its own
`OCCUPANCY.tsv`.

Two derived sets exist for deliberate reasons:

* `markers/<set>` (20 markers spanning occupancy 0.875–1.0) versus
  `markers/<set>_core` (8 complete markers) — a threshold cannot be shown to
  reject on data where everything is complete, and the core set keeps runs fast;
* `markers/degenerate/` — a 2-tip and a 1-tip marker carved out of a real one by
  subsetting, with a README saying so, for the "below the informative minimum"
  robustness cases.

### 3.4 The positive control

`markers/chimera_bench12/` differs from the clean set in exactly three files. In
each, the *Bacillus anthracis* Ames (Caryophanales) sequence is replaced by the
*Lactococcus lactis* Il1403 (Lactobacillales) ortholog — a cross-order transfer
at a known position, drawn from a recorded candidate pool with
`random.Random(20260929)`. Everything else is byte-identical, so the clean run is
the paired negative control for the same 12 genomes.

Detection is scored at the **exclusion** decision (`hgt_risk_level == level_3`),
because Level 2 is the suspicious band the pipeline keeps: counting it as a hit
would inflate apparent sensitivity. Two quantities are reported:

* **differential call** — markers excluded in the chimera run but *not* in the
  clean run, so the background exclusion rate is measured, not assumed;
* **direction** — for each planted marker, risk in the chimera run must exceed
  its risk in the clean run, and no unplanted marker may move at all.

## 4. Methodology of the end-to-end cases

1. **Observable consequences only.** A case asserts on a product, a column, an
   exit code, or a measured difference. "The function was called" is not
   accepted as evidence for a behavioural claim.
2. **Differential over absolute.** Most cases run the pipeline twice with one
   factor changed and compare (`--max-markers 2` vs `8`, threshold 0.0 vs 1.01,
   floor 0 vs 100000). An absolute expected value would silently encode the
   dataset; a difference isolates the knob.
3. **A knob that changes nothing is a failure.** `test_l2l3_boundary_decides_exclusion`,
   `test_occupancy_floor_rejects_markers_below_it`,
   `test_informative_site_floor_marks_markers_inconclusive_and_says_so` and
   friends exist to catch "read, recorded, then ignored" options — the defect
   class this codebase has been systematically audited for.
4. **Failure paths are first-class.** `test_v16_robustness_edges.py` requires:
   a nonzero exit, no traceback, a message naming the reason, and no empty
   placeholder products.
5. **The supported range is asserted, not assumed.** A dataset below four tips
   per marker cannot yield an inferrable gene tree; the case pins the refusal
   (exit ≠ 0, explained) and that no zero-byte tree file is published.
6. **Numbers used by the report are written by the tests** into
   `validation/results/metrics.json`, so prose cannot quietly drift from the run.
7. **Skips are reasoned and never a pass.** A gated criterion that refuses to run
   skips with the refusal quoted; the coverage matrix still requires some case to
   cover the option.

## 5. Environment

Required: Python 3.10–3.12 (ete3 cannot be imported on 3.13+, where the standard
library `cgi` module is gone), the declared runtime dependencies
(biopython, pyyaml, ete3, tomli), and the external tools `hmmsearch`, `mafft`,
`trimal`, `FastTree`/`fasttree`, `astral`, plus `iqtree3` for ML gene trees and
`checkm` if quality assessment is not skipped.

The suite pins its own environment (`validation/conftest.py::_tool_env`):
the invoking interpreter's `bin` is prepended to `PATH`, and `TMPDIR` /
`MAFFT_TMPDIR` point inside `validation/.work`. The second part is not
tidiness: MAFFT creates its scratch directory with `mktemp -dt` and, if that
fails, continues with an empty path, produces no alignment, and the pipeline
reports "no markers could be aligned" — an unwritable `/tmp` turning into a
wrong scientific verdict.

The measured environment of the published run (interpreter, package versions,
tool versions as each tool reports itself, host, CPU count) is archived in
`validation/results/environment.txt`.

## 6. How to run

```bash
# fast layer
pytest tests                      # ~35 s, no external tools needed
pytest tests/benchmark            # contract tests for the benchmark driver

# acceptance layer (real pipeline, real data)
python validation/run_validation.py -n 8            # excludes the slow case
python validation/run_validation.py --all -n 8      # includes detection metrics
python validation/run_validation.py --prepare       # rebuild data from NCBI first

# one feature area
pytest validation/cases/test_v06_hgt_screening.py -k composition -n 4
```

Interpreters other than the default are exercised by running the same suites
under each: 3.10, 3.11 and 3.12 are the measured, supported range, and the upper
bound is itself under test (`tests/unit/test_supported_range_is_earned.py`).

## 7. Results, findings and current limits

The measured outcome of the published run — pass/fail counts, timings, coverage
figures, detection metrics, and every product/schema defect the suite exposed
during this rework — is in
**[docs/TEST-REPORT.EN.md](TEST-REPORT.EN.md)** (中文:
[docs/TEST-REPORT.CN.md](TEST-REPORT.CN.md)).

Known limits, stated rather than left implicit:

* Marker inputs are a GTDB-Tk **layout** reconstruction (§3.3), so end-to-end
  results here do not measure GTDB-Tk's own marker choice.
* The functional-category map's text is each profile's `# DESCRIPTION` line: it
  exercises the `--cog-category-map` plumbing, not COG biology.
* Detection metrics rest on three planted chimeras over 12 genomes. That is a
  positive control with a known answer and a measured background rate, not an
  epidemiologically representative HGT sample; the wider external benchmark in
  `tests/benchmark/` (30 genomes, GTDB r53 truth) is the route to a stronger
  estimate and is user-triggered because of its download size.
* The 20-marker occupancy span reaches 0.875 at the low end: the bundled
  library is a conserved-marker library, so very sparse markers do not occur on
  this (closely related) genome set and are simulated by subsetting instead.
