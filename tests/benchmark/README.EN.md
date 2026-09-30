# Benchmark: external reproduction & chimera positive controls

Purpose: give MarkerFinder the ability to discover its own biases using
EXTERNAL data (G6). Nothing in this directory is auto-downloaded — GB-level
downloads and bioinformatics tools are the user's trigger.

## Contents

| Path | Role |
|---|---|
| `expected/accessions.yaml` | GTDB ar53 order-level set: 30 genomes (25 in-order + 5 out-group), per |
| `expected/chimera_positives.yaml` | Planted-chimera ground truth template (position-known HGT positives) |
| `expected/taxonomy_mustpass.yaml` | "must-pass" taxonomic relations skeleton |
| `fetch_datasets.sh` | Downloads the accession set from NCBI into `downloads/` (gitignored) |
| `make_chimeras.py` | Builds chimera FASTA: replaces marker X of genome A with the homologue from distant genome B |
| `run_benchmark.py` | Runs the pipeline, reads decision cards, computes metrics |
| `metrics.py` | precision / recall / F1 + block bootstrap CI (block = genome) |
| `provenance.py` |: DOI / licence / retrieval date / checksum record; `run_benchmark` refuses when incomplete |
| `ablation.py` | ****: reverse ablation — compares two rankings that differ in exactly ONE factor, refuses otherwise |
| `perf_probe.py` | runtime/scale probe against the shipped skeleton |
| `test_metrics_selfcheck.py` | pytest: metrics self-check incl. must-fail control |
| `test_provenance_contract.py` | pytest: gate is two-sided (shipped file must FAIL, a filled one must PASS) |
| `test_ablation_contract.py` | pytest: refusal rules + tau controls + driver wiring |

## How to run (user-triggered)

```bash
# 1. fetch genomes (GB-level; requires network + ~2 GB disk)
bash tests/benchmark/fetch_datasets.sh

# 2. build planted chimeras
python tests/benchmark/make_chimeras.py --downloads tests/benchmark/downloads

# 3. run the benchmark (requires mafft/trimal/FastTree/ASTRAL on PATH)
python tests/benchmark/run_benchmark.py --downloads tests/benchmark/downloads

# 4. reverse ablation — two completed runs, ONE declared factor change.
#    This stage needs no aligners: it reads two ranking TSVs and compares them.
python tests/benchmark/ablation.py \
    --run-a run_with_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --run-b run_without_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --factor hgt_filter=on:off --out-dir ablation_out

# or through the driver, which prints the same verdict before probing tools:
python tests/benchmark/run_benchmark.py --ablation-a ... --ablation-b ... \
    --ablation-factor hgt_filter=on:off
```

 rules worth knowing: omitting `--factor` is refused (a ranking difference
with no named cause is not an ablation), and passing two `--factor` flags is
refused too — the difference would be unattributable. Identical declared
settings are also refused, because in this codebase that historically meant the
switch was dead. Disagreement is reported as tie-corrected Kendall tau-b plus
per-marker rank displacement in `<prefix>.ablation.tsv`; markers whose score is
`NA` are excluded and named on stderr rather than ranked as zero.

If mafft/trimal/FastTree/ASTRAL are missing, run_benchmark prints
`NOT EXECUTED` per stage and exits non-zero — a benchmark that silently
"passes" without tools is exactly the failure we must not reproduce.

## Sources / licence / citation

* Genomes: NCBI RefSeq accessions listed in `expected/accessions.yaml`
  (selection rule per: one GTDB ar53 order, 25 representative genomes
  within the order + 5 representatives from outside as out-group).
* Taxonomy truth: GTDB release r53 (Parks et al. 2021, doi:10.1038/s41587-021-00960-6).
* Downloaded data is NOT committed; `downloads/` is gitignored and
  `fetch_datasets.sh` records SHA-256 per file into `downloads/manifest.sha256`.
