# Benchmark: external reproduction & chimera positive controls

[中文文档](README.CN.md)

Purpose: give MarkerFinder the ability to discover its own biases using
EXTERNAL data. Nothing in this directory is auto-downloaded — GB-level
downloads and bioinformatics tools must be triggered explicitly by the user.

## Contents

| Path | Role |
|---|---|
| `expected/accessions.yaml` | GTDB ar53 order-level set: 30 genomes (25 in-order + 5 out-group), per the rule under [Sources / licence / citation](#sources--licence--citation) |
| `expected/chimera_positives.yaml` | Planted-chimera ground truth template (position-known HGT positives) |
| `expected/taxonomy_mustpass.yaml` | "must-pass" taxonomic relations skeleton |
| `fetch_datasets.sh` | Downloads the accession set from NCBI into `downloads/` (gitignored) |
| `make_chimeras.py` | Chimera ground-truth gate: validates `expected/chimera_positives.yaml` and `downloads/manifest.sha256` so the swap cannot silently no-op; the per-marker swap itself (genome A's marker X → the ortholog from distant genome B) is NOT implemented |
| `run_benchmark.py` | Probes the tools, the downloads manifest and the provenance record; the pipeline → decision card → metrics stages it gates are NOT implemented, so even a run that passes every gate ends in `NOT EXECUTED` |
| `metrics.py` | precision / recall / F1 + block bootstrap CI (block = genome) |
| `provenance.py` | The `provenance:` block of [`expected/accessions.yaml`](expected/accessions.yaml): DOI / licence / retrieval date / checksum record; `run_benchmark` refuses when incomplete |
| `ablation.py` | **step 4 of [How to run](#how-to-run-user-triggered), reverse ablation**: compares two rankings that differ in exactly ONE factor, refuses otherwise |
| `perf_probe.py` | performance-budget probe: times the added pure-Python paths on synthetic data (60 markers x 60 taxa, seed `20260920`) against a 60 s ceiling; it reads no shipped data, calls no external tool, and leaves the "5% of total run time" denominator to the user's own log |
| `test_metrics_selfcheck.py` | pytest: metrics self-check incl. must-fail control |
| `test_provenance_contract.py` | pytest: the [`provenance.py`](provenance.py) gate is two-sided (shipped file must FAIL, a filled one must PASS) |
| `test_ablation_contract.py` | pytest: [`ablation.py`](ablation.py) refusal rules + tau controls + driver wiring |

## How to run (user-triggered)

```bash
# 1. fetch genomes (GB-level; requires network + at least 5 GB of free disk — below that the script REFUSES and exits 2)
bash tests/benchmark/fetch_datasets.sh

# 2. build planted chimeras
python tests/benchmark/make_chimeras.py --downloads tests/benchmark/downloads

# 3. run the benchmark (requires mafft/trimal/FastTree/fasttree/astral on PATH; iqtree3 is NOT in this required list)
python tests/benchmark/run_benchmark.py --downloads tests/benchmark/downloads

# 4. reverse ablation — two completed runs, ONE declared factor change.
#    This stage needs no aligners: it reads two ranking TSVs and compares them.
python tests/benchmark/ablation.py \
    --run-a run_with_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --run-b run_without_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --score-column marker_quality_score \
    --factor hgt_filter=on:off --out-dir ablation_out

# or through the driver, which prints the same verdict before probing tools:
python tests/benchmark/run_benchmark.py --ablation-a ... --ablation-b ... \
    --ablation-factor hgt_filter=on:off
```

[`ablation.py`](ablation.py) rules worth knowing: omitting `--factor` is refused (a ranking difference
with no named cause is not an ablation), and a second `--factor` naming a
different factor is refused too — the difference would be unattributable.
Identical declared settings are also refused, because in this codebase that
historically meant the switch was dead — `ablation.py` names three that were:
`--scan-stability-min`, `--min-informative-sites` and the `composition` step.
Disagreement is reported as tie-corrected Kendall tau-b plus
per-marker rank displacement in `<prefix>.ablation.tsv`; markers whose score is
`NA` are excluded and named on stderr rather than ranked as zero. `--score-column`
defaults to `overall_score`, a column `markerfinder.marker_summary.tsv` does not
have, so against a real product it must be named; the driver form above does not
forward it and therefore refuses on a real `marker_summary.tsv`.

If any of mafft/trimal/FastTree/fasttree/astral is missing, `run_benchmark` prints
`NOT EXECUTED` per stage and exits non-zero — a benchmark that silently
"passes" without tools is exactly the failure we must not reproduce.

## Sources / licence / citation

* Genomes: NCBI RefSeq accessions listed in `expected/accessions.yaml`
  (selection rule per [`expected/accessions.yaml`](expected/accessions.yaml): one GTDB ar53 order, 25 representative genomes
  within the order + 5 representatives from outside as out-group).
* Prerequisite: the shipped `expected/accessions.yaml` is still a placeholder
  template — `order` reads `REPLACE_WITH_ORDER_NAME` and all five
  `provenance:` fields (`source_registry`, `doi`, `license`, `retrieved_utc`,
  `checksum_sha256`) are empty. Until you choose the concrete GTDB ar53 order,
  fill in the accessions and those five fields, `fetch_datasets.sh` REFUSES to
  download and exits 2, and `run_benchmark` likewise refuses to run on an
  incomplete provenance record. Passing every gate is still not enough for a
  number: the swap and the scoring stages do not exist, so `run_benchmark` ends
  by printing `NOT EXECUTED` and exiting 2 (see the two rows above).
* Taxonomy truth: GTDB release r53 (Parks et al. 2021, doi:10.1038/s41587-021-00960-6).
* Downloaded data is NOT committed; `downloads/` is gitignored and
  `fetch_datasets.sh` records SHA-256 per file into `downloads/manifest.sha256`
  (the script calls `shasum -a 256` with no `sha256sum` fallback; on Linux
  `shasum` is the Perl wrapper — present wherever `perl` is installed, absent
  from minimal images — so where it is missing the script's `set -euo pipefail`
  aborts the run at the manifest step).
