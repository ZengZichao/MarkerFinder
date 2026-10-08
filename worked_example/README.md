# Worked example — reproduction recipe

**Status of this directory: recipe only.** The repository contains **no output of the
six-genome archaeal run** quoted in the manuscript's worked example. Nothing here
records a runtime, an alignment length, a marker count or a tree statistic, because
none exists in the tree to record it from; see
`审阅报告_MarkerFinder_2026-10-08.md` batch note X-015 / FACT-020, which searched the
whole repository for `176 s`, `9,771`, `Halobacteria` and the marker counts and found
zero hits. `validation/data/` holds bacterial proteomes only (`PROVENANCE.tsv` lists
twelve, all Bacteria), and the committed run sets are `quad4` and `small8`.

Filling this gap requires **one real run by the author**, after which the artefacts
listed under *Artefacts to archive* belong in this directory (or in the data deposit)
and the manuscript's numbers get re-checked against them. Until that happens the
worked-example figures in the paper are unverifiable from the software release, and no
amount of documentation here changes that.

## What the example is

Six GTDB r220 domain-Archaea representative genomes — three class *Halobacteria*
(order Halobacteriales) and three class *Methanosarcinia* — taken as `.faa` protein
predictions, run end to end in HMM marker mode. The accessions are those listed in
Supplementary Information **Table S6c**; they are repeated here so the recipe is
self-contained.

| GTDB accession | NCBI assembly | Class; order | Species |
|---|---|---|---|
| `RS_GCF_000517625.1` | `GCF_000517625.1` | Halobacteria; Halobacteriales | *Halostagnicola larsenii* |
| `RS_GCF_010119195.1` | `GCF_010119195.1` | Halobacteria; Halobacteriales | *Haloarcula salina* |
| `RS_GCF_002844195.1` | `GCF_002844195.1` | Halobacteria; Halobacteriales | *Halegenticoccus soli* |
| `RS_GCF_000970325.1` | `GCF_000970325.1` | Methanosarcinia; Methanosarcinales | *Methanococcoides methylutens_A* |
| `RS_GCF_000711905.1` | `GCF_000711905.1` | Methanosarcinia; Methanosarcinales_A | *Methermicoccus shengliensis* |
| `RS_GCF_000204415.1` | `GCF_000204415.1` | Methanosarcinia; Methanotrichales | *Methanothrix soehngenii* |

Protein counts, total amino acids and `.faa` SHA-256 prefixes per genome are in Table
S6c; recompute them locally rather than trusting the table when you re-run.

## Prerequisites

External tools on `PATH`: `hmmsearch` (HMMER 3.4), `mafft`, `trimal`, `iqtree3` (or
`fasttree` with `--fast-tree`), `astral`, and optionally `checkm`. Python 3.10–3.12
(3.13 is excluded because `ete3` needs the removed `cgi` module).

## Step 1 — obtain the marker HMM library

`db/gtdb_markers/` is **git-ignored**: the TIGRFAM/Pfam profiles are third-party work
and are assembled locally, so a package built from git does not contain them and
`--marker-mode hmm` cannot run from a bare checkout. Build them first:

```bash
python scripts/fetch_marker_db.py --auto          # split a local GTDB-Tk install
python scripts/fetch_marker_db.py --auto --record # and stamp the reference hash
```

The ar53 set is 53 per-marker profiles (12 `.hmm` + 41 `.HMM`); `db/expected_hashes.json`
holds the SHA-256 the release was validated against. A mismatch is a **warning during a
run** (the run continues) but a **FAIL for `markerfinder --check`**, which exits non-zero —
see `_note_gtdb_markers` in that file.

## Step 2 — obtain the six genomes

Download each assembly's predicted protein set (`.faa`) into `genomes/`, named after
the accession, e.g. `genomes/RS_GCF_000517625.1.faa`. `provenance_manifest.json` lists,
per genome, the accessions, the GTDB r220 lineage, the download URL, the retrieval date,
the SHA-256 and the sequence counts the manuscript quotes, so a reviewer can fetch the
identical files and check them by digest; see *Provenance record* below.

## Step 3 — supply the taxonomy table

`--taxonomy-table` takes a TSV of lineage assignments. The example uses the GTDB r220
archaeal marker taxonomy file (`ar53_taxonomy_r220.tsv`) filtered to the six
accessions, in the Format-A (`table` / `reverse`) layout `markerfinder` parses.

## Step 4 — run

```bash
markerfinder -i genomes/ -o output/ \
  --marker-mode hmm \
  --taxonomy-table taxonomy.tsv \
  --mode standard \
  -t 4
```

What those flags do — stated precisely, because two of them are commonly
over-read (see `README.EN.md` → *Analysis Modes* / *Marker Presets*):

* `--marker-mode hmm` selects the profile-HMM marker route and is what makes the
  ar53 library load-bearing. It is also the **only** mode in which the non-default
  selection strategies are reachable for a real dataset.
* `--mode standard` selects the **selection strategy only** (Greedy). It does *not*
  set the occupancy floor or the marker budget; those keep their configuration
  defaults (0.75 / 60) or the Phase-0 adaptive values.
* `-t 4` is the thread count for the parallel stages. ASTRAL is invoked `-t 1`
  regardless (pinned in `markerfinder/modules/phylogenetic_inference.py`).

Because the input is protein-only, CheckM cannot run and the pipeline records
`Quality source: default_no_nucleotide` in the summary — the quality tier is then an
assumption from file naming, not a measurement. That is the expected state for this
example, not an error.

## Artefacts to archive with the paper

Copy these out of the run directory (names are the pipeline's own products, written
with the default `markerfinder.` prefix). `output/` is git-ignored, so committing here
needs `git add -f`, or archive the run in the data deposit instead.

| Artefact | Path as written | Why it is the evidence |
|---|---|---|
| Parameter snapshot | `Phase5_metadata/run_config.json` | the exact invocation, versions and database hashes; replayable with `--config` |
| Headline summary | `Phase5_reports/markerfinder.pipeline_summary.txt` | marker counts, alignment length, runtime, quality source, species-tree source |
| Marker inventory | `Phase5_reports/markerfinder.marker_summary.tsv` | per-marker occupancy and selection record |
| HGT audit | `Phase5_reports/markerfinder.hgt_evaluation.tsv` | per-marker risk components and exclusion decision |
| Assertion log | `Phase5_reports/markerfinder.assertions.tsv` | which output assertions fired |
| Decision cards | `Phase5_evidence/decision_*.json` | per-marker provenance behind the HGT calls |
| Trees | `Phase4_trees/markerfinder.species_tree_concat.newick`, `...species_tree_astral.newick` | the two topologies the RF/quartet statistics compare |
| Provenance record | `provenance_manifest.json` (ships with the release) | registry, licence, retrieval date, download URL and SHA-256 per genome, plus the counts Table S6c quotes |

## Provenance record

`tests/benchmark/provenance.py` defines what "provenance" has to mean here; its
`REQUIRED_FIELDS` are the contract, and placeholder values (`tbd`, `n/a`, empty) are
rejected on purpose so a blank record cannot pass as a filled one.

`provenance_manifest.json` is that record for the six worked-example genomes. It follows
`validation/data/PROVENANCE.tsv`, the shipped provenance file for the bacterial
validation set, and carries per genome: `source_registry`, `license`, `retrieved_utc`,
`checksum_sha256`, `gtdb_release`, the GTDB and NCBI accessions, the GTDB r220 lineage,
both download URLs, and the protein count and total residue count that Table S6c of the
Supporting Information quotes.

One field of the contract has no honest value to write here. Its `doi` is defined in
`tests/benchmark/expected/accessions.yaml` as the "DOI of the dataset/release you actually
downloaded" — a release-level identifier, which the bacterial benchmark file leaves as
`doi: ""` until whoever runs it chooses the concrete order and fills it. A RefSeq
assembly has no data DOI of its own, and `validation/data/PROVENANCE.tsv` likewise has no
`doi` column. So this record identifies the source by `gtdb_release`, `source_registry`,
the assembly accessions and the download URLs, and leaves a release DOI to whoever mints
or cites one; nothing here is guessed. `tests/benchmark/provenance.py` enforces its
`REQUIRED_FIELDS` on `accessions.yaml` and on the benchmark runners, not on this file.

The manifest's counts and digests are measured from the files themselves; it was written
by a generator that refuses to emit a record whose counts disagree with Table S6c.

## What still needs the author

1. **One real run** at the release the paper cites, and the artefacts above archived
   with it. The numbers in the manuscript's worked example cannot be verified,
   corrected, or reproduced from this repository until that is done.
2. **Confirmation of the accession list** against Table S6c after any GTDB update.
3. A decision on whether `worked_example/` ships inside the software release or in the
   data deposit (dryad/Zenodo) referenced by the Data Availability statement.
