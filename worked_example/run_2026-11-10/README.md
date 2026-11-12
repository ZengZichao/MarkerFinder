# Worked example — archived run (2026-11-10)

This directory contains the complete output of the six-genome worked example
quoted in the manuscript (Section 3; Table S6). It was produced with the
released code (`markerfinder 0.1.1`) and the inputs listed below.

## Command

```bash
markerfinder -i genomes/ -o output/ \
  --marker-mode hmm \
  --taxonomy-table taxonomy.tsv \
  --mode standard \
  -t 4
```

The `db/` directory used by the run was the repository's `db/gtdb_markers/ar53`
set (53 profile HMMs, 6.3 MB), assembled locally from a GTDB-Tk 2.4.1
installation. The run logs a database-hash warning because that build's SHA-256
(`84cf08b9…`) differs from the reference recorded for the earlier validation
build (`67e931f3…`); this is the documented behaviour for a different upstream
build, not a failed run (see `output/Phase5_metadata/run_config.json`).

## Environment

- Host: Linux x86_64 (AMD EPYC 9654), 4 threads requested
- Python 3.11.0; Biopython 1.88; PyYAML 6.0.3; ete3 3.1.3
- HMMER 3.4 (Aug 2023)
- MAFFT v7.526 (2024/Apr/26)
- trimAl v1.5.rev0 build[2024-05-27]
- FastTree 2.1.11 (double precision)
- IQ-TREE 3.0.1
- ASTRAL 5.7.8

Exact versions are in `tool_versions.txt` and `run_config.json`.

## Result

- 6 genomes, 53 ar53 candidate markers, 44 selected at occupancy ≥ 0.75
- HGT screening: 43 Level 1, 0 Level 2, 1 Level 3 (`PF07541.13`)
- Mean HGT risk: 0.0227 (43 × 0, 1 × 1.0); HGT evidence coverage 1.00
- Supermatrix: 9,771 amino-acid sites, 43 partitions
- Concatenation vs coalescent: normalised RF 0.0000; quartet agreement 1.0000
- Runtime: **238 s** on four threads (this Linux host; the earlier Apple
  Silicon run quoted in the manuscript took 176 s)
- Species tree source: ASTRAL

The 44 risk values take only two values (43 at 0 and one at 1.0), so the
Level 2 band is empty. This is expected with six genomes in two classes
(two informative taxa for the monophyly proportion); the example is a smoke
test, not a calibration of the 0.25/0.60 thresholds.

## Files

- `genomes/` — the six `.faa` files, SHA-256-verified against
  `../provenance_manifest.json` and Table S6c
- `taxonomy.tsv` — the six-line GTDB r220 lineage table in Format A
- `tool_versions.txt` — versions of the external tools used
- `run_config.json` — the parameter snapshot written by the pipeline
- `output/Phase4_alignments/` — the partition file (9,771 sites)
- `output/Phase4_trees/` — concatenation tree, ASTRAL tree, all 53 gene trees
  with per-tree input SHA-256 stamps
- `output/Phase5_reports/` — HTML report, marker summary, HGT evaluation,
  assertion log and pipeline summary
- `output/Phase5_metadata/` — the run configuration
- `output/Phase5_evidence/` — decision cards for the screened markers
- `run_output.txt` — the successful run's log
- `run_config.json` — a copy of the parameter snapshot

`SHA256SUMS` records the digest of every file in this directory.
