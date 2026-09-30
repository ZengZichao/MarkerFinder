# MarkerFinder Validation Suite

**Full-feature, full-workflow validation on real data.**
[中文文档](README.CN.md)

This directory is the acceptance evidence that ships with MarkerFinder. It runs
the real pipeline against real genomes and real HMM profiles — nothing here is
mocked, and nothing here is a unit test. The fast, mock-based layer lives in
[`../tests/`](../tests).

---

## What it establishes

| | Claim being established | Where |
|---|---|---|
| 1 | Every capability the command line advertises has a case | `cases/test_v90_capability_coverage.py`, `results/capability_matrix.tsv` |
| 2 | The whole workflow runs, step by step and in one shot, with the same result | `cases/test_v12_stepwise_workflow.py` |
| 3 | Every documented output file exists, with the documented columns | `cases/test_v14_product_audit.py` |
| 4 | Detection performance is measured against a positive control with a known answer | `cases/test_v17_benchmark_detection.py` |
| 5 | Two runs on the same input give the same answer | `cases/test_v15_reproducibility.py` |
| 6 | Failure paths fail loudly, in the right exit code, without a traceback | `cases/test_v16_robustness_edges.py` |
| 7 | The test data itself is what it says it is | `cases/test_v01_environment_and_selfcheck.py` |

---

## Layout

```
validation/
├── README.EN.md / README.CN.md   this document (English is primary)
├── capabilities.py               the capability surface, read from the code
├── conftest.py                   fixtures: runners, inputs, per-case outputs
├── run_validation.py             one-command entry point
├── cases/                        test_vNN_*.py, numbered by feature area
├── data/                         fixtures that ship + data the suite fetches
│   ├── PROVENANCE.tsv            per-assembly organism, tax_id, lineage, checksum
│   ├── MANIFEST.sha256           checksum of every file in this directory
│   ├── genomes/                  12 NCBI RefSeq proteomes (.faa) — FETCHED, not committed
│   ├── sets/                     which accessions form each named input set
│   ├── markers/                  per-marker FASTA, occupancy tables, chimeras
│   ├── hmms/                     TIGRFAM/Pfam profiles for --marker-mode hmm
│   ├── taxonomy/                 Format A / Format B tables, malformed controls
│   ├── trees/                    reference species trees and illegal references
│   ├── checkm/                   quality tables taken from NCBI's own records
│   ├── cog/                      marker → functional-category map
│   ├── mustpass/                 taxonomy must-pass baselines (pass and fail)
│   ├── configs/                  YAML / TOML / JSON config files
│   ├── sequences/                protein and nucleotide side files
│   └── variants/mag_named/       identifier-scheme variant for the MAG pathway
├── scripts/                      how the data was built, from public sources
│   ├── 01_fetch_genomes.py       NCBI download + authoritative lineage
│   ├── 02_build_marker_sets.py   per-marker FASTA, occupancy, planted chimeras
│   └── 03_build_fixtures.py      tables, trees, quality files, configs, manifest
└── results/                      archived outputs of the published run
    ├── metrics.json              every number quoted by the test report
    ├── capability_matrix.tsv     capability → covering cases
    ├── chimera_metrics.json      detection precision / recall / F1
    ├── junit.xml                 per-case pass/fail record of the run
    └── environment.txt           interpreter, packages and tool versions measured
```

`data/.work/` (actually `.work/`) is scratch: materialised input directories,
per-case outputs and temp directories. It is disposable and git-ignored.

---

## Running it

Requirements: the MarkerFinder environment (Python 3.10–3.12 plus
`hmmsearch`, `mafft`, `trimal`, `FastTree`/`fasttree`, `astral`, `iqtree3`), and
`pip install -e ".[dev]"`.

```bash
# Everything except the slow 12-genome detection case, 8 parallel workers
python validation/run_validation.py -n 8

# Including the slow case (real detection metrics on 12 genomes)
python validation/run_validation.py --all -n 8

# Rebuild the test data from NCBI first, then run
python validation/run_validation.py --prepare
```

Equivalent direct pytest invocations:

```bash
pytest validation/cases -m "not slow" -n 8
pytest validation/cases -m slow            # detection metrics; ~10 min
pytest validation/cases/test_v06_hgt_screening.py -k composition
```

The suite pins its own environment: the interpreter's `bin` is prepended to
`PATH`, and `TMPDIR`/`MAFFT_TMPDIR` point inside `validation/.work`, so a
read-only or full `/tmp` cannot turn into a misleading "no markers could be
aligned".

### Runtime

Measured on 96 cores with `-n 14`: the non-slow suite completes in about 15
minutes. A single pipeline run on the default four-genome, four-marker
configuration costs roughly 40 seconds; the same dataset over 30 markers costs
7 m 19 s, which is why feature-level cases run on small marker budgets and only
the occupancy and detection cases use the larger sets.

---

## The test data, honestly described

**Genomes.** 12 complete NCBI RefSeq assemblies, protein FASTA (`--include
protein`), spanning two orders of the phylum Bacillota and one Pseudomonadota
outgroup. `data/PROVENANCE.tsv` records for each: organism name, tax_id, the
seven NCBI lineage ranks, assembly level, contig count, length, GC%, protein
count, NCBI's own CheckM completeness/contamination, the download URL, the
retrieval date and the SHA-256 of the shipped file.

Lineages are **read from NCBI**, never typed. An earlier revision of this
dataset hand-typed its taxonomy table and labelled *Salmonella enterica*
`GCF_000008105.1` as `o__Bacillales;f__Bacillaceae;g__Bacillus`. Because the
taxonomy table drives the monophyly screen, that was an invalid experiment, not
a typo. `test_shipped_genomes_are_the_organisms_the_provenance_claims` now
checks each shipped proteome against its recorded genus, and
`01_fetch_genomes.py` cannot produce a table without a `tax_id`.

**Per-marker FASTA (`markers/`).** GTDB-Tk writes one unaligned FASTA per marker
into `align/marker_genes/`, and `--marker-mode gtdb_tk` consumes that layout.
GTDB-Tk itself needs a ~100 GB reference-tree package, which cannot ship in a
test bundle, so these files are **a reconstruction of the layout**: each genome's
own best HMM hit (bitscore ≥ 20) against the TIGRFAM/Pfam profiles in
`db/gtdb_markers/bac120` (assembled by `scripts/fetch_marker_db.py`), one file
per marker, filename stem = marker id, tips
named by assembly accession. What is preserved is everything the pipeline reads;
what is not reproduced is GTDB-Tk's taxonomy-aware best-hit selection.
`markers/OCCUPANCY_ALL.tsv` records the occupancy of every profile scanned, and
each directory carries its own `OCCUPANCY.tsv`.

**Planted chimeras (`markers/chimera_bench12`).** In three markers the
*Bacillus anthracis* Ames sequence is replaced by the *Lactococcus lactis*
Il1403 ortholog — a cross-order transfer at a known position, drawn from the
candidate pool with `random.Random(20260929)` (the pool and the rule are recorded
in `markers/chimera_ground_truth_bench12.json`). Every other byte is identical to
the clean directory, which makes the clean run the negative control for the same
12 genomes.

**Quality values (`checkm/`).** Completeness and contamination are NCBI's own
`checkm_info` figures per assembly, not invented numbers. Assemblies for which
NCBI records none are omitted, with the omission written into the file.

**Functional categories (`cog/`).** The category text is the `# DESCRIPTION`
line of the marker's HMM profile. This exercises the plumbing of
`--cog-category-map`; it is not an assertion about COG biology, and the file
says so.

**MAG-named variant (`variants/mag_named/`).** The four quad4 genomes renamed to
`mag_01..mag_04`, with marker headers and the taxonomy table rewritten to match.
Sequences are unchanged — the classification rule keys on the file name, so this
is a labelling experiment and is recorded as one.

---

## Conventions the cases follow

* Every case gets its own output directory; the suite is safe under `-n`.
* Cases assert **observable** consequences: a product, a column, an exit code,
  or a difference between two runs that differ in exactly one knob. A knob that
  changes nothing is reported as a dead control rather than a passing test.
* Differential assertions (run A vs run B, one factor changed) are preferred
  over absolute thresholds, because an absolute threshold silently encodes the
  dataset.
* Skips are loud and reasoned; a skip is never presented as a pass.
* Numbers that the report quotes are written to `results/metrics.json` by the
  cases themselves, so a stale number in prose becomes a diff in the record.
