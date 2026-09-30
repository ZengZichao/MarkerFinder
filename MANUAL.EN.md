# MarkerFinder Manual

**Adaptive HGT-Aware Phylogenomic Pipeline — Technical Reference Manual**

Version 0.1.0

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Installation Guide](#2-installation-guide)
3. [Configuration Reference](#3-configuration-reference)
4. [Pipeline Phases — Detailed Reference](#4-pipeline-phases--detailed-reference)
5. [External Tool Integration](#5-external-tool-integration)
6. [Database Setup](#6-database-setup)
7. [Output Files — Complete Reference](#7-output-files--complete-reference)
8. [Advanced Usage](#8-advanced-usage)
9. [Troubleshooting](#9-troubleshooting)
10. [Algorithm Details](#10-algorithm-details)
11. [Self-Check Artifacts Quick Reference](#self-check-artifacts-quick-reference)
12. [Glossary](#glossary)

---

## 1. Introduction

### 1.1 What is MarkerFinder?

MarkerFinder is a phylogenomic pipeline that automates marker gene selection, horizontal gene transfer (HGT) screening, and species tree inference for prokaryotic genomes. It provides a dynamic, data-driven alternative framework to traditional fixed-marker approaches (e.g., 16/27/37/38 CSCG sets), while preserving the biological basis of established marker collections.

### 1.2 Design Principles

| Principle | Description |
|-----------|-------------|
| **Modular** | Five innovation paths operate as independent modules within a unified pipeline |
| **Adaptive** | Marker set and parameter thresholds adjust dynamically based on input data |
| **HGT-Aware** | Standalone Phase 2 phylogenetic HGT screening (MAD rooting + monophyly proportion or RF/quartet) |
| **Dual-Strategy** | Parallel concatenation (Supermatrix) and coalescent (ASTRAL) tree inference |
| **MAG-First** | Native support for metagenome-assembled genomes with quality-aware analysis |
| **Engineered** | Modular test suite (case counts live in the README badge; stated in one place only). The integration layer runs the pipeline under mocked external tools and reads the artifacts back against the acceptance clauses; the benchmark layer ships fetch/chimera/metrics/provenance/ablation scripts. Neither is a placeholder any more |

### 1.3 Five Innovation Paths

| Path | Module | Innovation |
|------|--------|------------|
| 1 | `marker_selection` | Adaptive marker gene selection with 4 strategies |
| 2 | `hgt_filter` | Phylogenetic HGT screening with marker grading (MAD rooting + monophyly proportion or RF/quartet) |
| 3 | `phylogenetic_inference` | Dual-strategy concatenation + coalescent with conflict detection |
| 4 | `report_generator` | Static HTML report + plain text output formats |
| 5 | `mag_optimization` | Quality-aware adaptive parameters for MAG datasets |
| — | `ortholog_resolver` | Multi-copy gene ortholog resolution (BBH, graph clustering) |

---

## 2. Installation Guide

### 2.1 System Requirements

- **Operating System:** Linux (recommended), macOS
- **Python:** >= 3.10
- **Disk Space:** the self-contained mode needs only ~1–12 GB, plus a marker HMM of a few MB; the tool itself is ~500 MB. HGT screening uses gene-tree incongruence, so it requires no large external protein database
- **Memory:** Minimum 4 GB, recommended 16 GB for datasets > 50 genomes

### 2.2 Conda Environment Setup

```bash
# Create and activate environment
conda create -n markerfinder python=3.10 -y
conda activate markerfinder

# Install bioinformatics tools (blast is not needed because the main flow skips BBH/graph clustering; diamond is optional as a reserved interface)
conda install -c bioconda -c conda-forge \
    hmmer>=3.3.2 \
    mafft>=7.520 \
    trimal>=1.4 \
    iqtree>=3.0.0 \
    astral-tree>=5.7.0 \
    fasttree>=2.1.11 \
    checkm-genome>=1.2 \
    ete3

# Install diamond only if you plan to use the OrthologResolver interface in future versions
# conda install -c bioconda diamond>=2.1.0

# Install Python dependencies
pip install biopython>=1.81 pandas>=2.0 numpy>=1.24 \
    scipy>=1.10 rich>=13.0

# Install MarkerFinder
pip install -e ".[dev]"
```

### 2.3 Verification

```bash
# Verify MarkerFinder
markerfinder --version
# Expected: markerfinder 0.1.0
# (the version is a constant in markerfinder/_version.py, so every install path
#  and every unpacked source tree report the same string)

# Verify external tools
hmmsearch -h | head -1         # HMMER 3.3.2+
mafft --version                 # MAFFT 7.520+
trimal --version                # trimAl 1.4+
iqtree3 --version               # IQ-TREE3 3.0.0+
diamond version                 # DIAMOND 2.1.0+

# Run unit tests (case counts are stated in the README "Testing" section)
pytest tests/unit -q

# Environment self-check (49 items, or 48 items without the fetched reference
# databases; non-zero exit if any fails)
markerfinder --check
```

### 2.4 Directory Structure

```
MarkerFinder/
├── markerfinder/              # Main package
│   ├── __init__.py            # Public API exports
│   ├── __main__.py            # CLI entry point
│   ├── _version.py            # Version and authorship (single constant)
│   ├── assertions.py          # Output-assertion registry and adjudication
│   ├── banner.py              # Startup banner display
│   ├── config.py              # Configuration dataclasses (10 config classes)
│   ├── config_loader.py       # YAML/TOML/JSON config file loading
│   ├── exceptions.py          # Custom exception hierarchy
│   ├── phases.py              # Pipeline phase registry (directory names, labels)
│   ├── pipeline.py            # Pipeline coordinator
│   ├── taxonomy.py            # Taxonomy name parsing
│   ├── validation.py          # Tree/sequence/cross-validation
│   ├── cli/                   # Argument parsing, subcommands, self-check
│   ├── models/                # Core data models
│   ├── modules/               # Five innovation path modules
│   └── utils/                 # Utility functions
├── tests/                     # Test suite (module/case counts live in the README)
│   ├── unit/                  # Unit tests (incl. the must-fail control matrix)
│   ├── integration/           # Integration tests (state recovery + product-level acceptance under mocked tools)
│   ├── data/                  # Small hand-written inputs (FASTA, Newick, taxonomy)
│   ├── fixtures/              # Reference trees and must-fail control payloads
│   └── benchmark/             # External reproduction + planted chimeras
├── validation/                # Full-feature suite: real pipeline on real genomes
│   ├── cases/                 # Acceptance cases, one file per feature area
│   ├── data/                  # Shipped fixtures + provenance (genomes are fetched)
│   ├── results/               # Archived evidence of the published run
│   └── scripts/               # Data pipeline: fetch genomes / marker sets / fixtures
├── db/                        # expected_hashes.json (the HMM profiles are fetched)
├── docs/                      # Test plan and measured test report (EN + CN)
├── scripts/                   # Utility scripts (fetch_marker_db.py, smoketest.py, …)
├── pyproject.toml             # Package configuration
├── environment.yml            # Conda environment with the external tools
├── config.example.yaml        # Example config file
├── LICENSE                    # MIT License
├── README.md                  # Landing page: links to the EN and CN documents
├── README.CN.md / README.EN.md
└── MANUAL.CN.md / MANUAL.EN.md
```

### 2.5 Step-by-Step Execution (Subcommand Mode)

MarkerFinder can be split into 4 subcommands run sequentially. State persists under the output directory (`.markerfinder/.pipeline_state.json` + `Phase5_metadata/context.json`) between steps, enabling staged debugging and reuse of intermediates.

| Subcommand | Phase | Description |
|------------|-------|-------------|
| `scan` | Phase 0+1+1.5 | quality preprocessing + marker scanning + sequence extraction (requires `-i`) |
| `filter` | Phase 2 | Phylogenetic HGT screening + marker level assignment (caches gene trees to `<output>/Phase4_trees/gene_trees/`) |
| `infer` | Phase 3 | phylogenetic inference (supermatrix + coalescent; reuses gene trees cached by `filter`) |
| `report` | Phase 4 | generates `Phase5_reports/` HTML/text reports, `Phase4_trees/` species trees, `Phase4_alignments/` partition files + `Phase5_metadata/run_config.json` |

```bash
# Run step-by-step (state auto-passes under <output>/.markerfinder/)
markerfinder scan     -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes
markerfinder filter             -o output/ -t 8
markerfinder infer              -o output/ -t 8 --gene-tree-builder fasttree
markerfinder report            -o output/

# The subcommand may also be placed after the options (still need --gtdb-markers-dir for scan)
markerfinder -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes scan
```

**Notes:**
- Without a subcommand (`markerfinder -i... -o...`), the pipeline still runs as a single aggregate `run` for backward compatibility.
- Each step validates that its prerequisite step has completed; e.g. running `filter` without first running `scan` raises a clear error.
- Gene-tree cache lives at `<output>/Phase4_trees/gene_trees/{marker_id}.nwk` and is shared between `filter` and `infer` to avoid rebuilding trees. Each entry is accompanied by `{marker_id}.nwk.input_sha256`, the SHA-256 of the (genome id, sequence) pairs the tree was built from: a cache entry is reused only when that stamp matches the sequences of the current run, so a rerun into an existing output directory (`--force`) cannot silently infer from another run's tree.

---

## 3. Configuration Reference

### 3.1 PipelineConfig

Top-level configuration object. Created from CLI arguments or a config dict.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_dir` | `str` | `""` | Input genome directory |
| `output_dir` | `str` | `"./output"` | Output directory |
| `output_prefix` | `str` | `"markerfinder"` | Output file prefix |
| `cpus` | `int` | `1` | Thread count |
| `tmp_dir` | `Optional[str]` | `None`(auto-generated per-run temp dir) | Temporary file directory. Default `None` = auto-create `markerfinder-{uuid8}` under system temp dir — **per-run isolated, auto-cleaned on exit** to avoid domtblout / concat file contention when parallel runs share `--tmp-dir`. Keep with `--keep-tmp`. Explicit path (e.g. `--tmp-dir./tmp`) is honored and not auto-cleaned. |
| `mode` | `str` | `"standard"` | Analysis mode |
| `db_dir` | `str` | `"./db"` | Database directory |
| `force` | `bool` | `False` | Force overwrite existing output files |
| `no_clobber` | `bool` | `False` | Skip existing files without overwriting |
| `sequences_path` | `Optional[str]` | `None` | Reference sequence file path |
| `mol_type` | `Optional[str]` | `None` | Molecule type (`None`=auto-detect, or `DNA`/`RNA`/`protein`) |
| `skip_length_check` | `bool` | `False` | Skip sequence length consistency check |
| `strip_annotations` | `bool` | `False` | Strip NHX annotations from tree |
| `keep_tmp` | `bool` | `False` | Keep the temporary directory after pipeline completion (for debugging). |

### 3.2 SelectionConfig

Controls Phase 1 marker gene selection.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `min_occupancy` | `float` | `0.75` | Minimum single-copy occupancy threshold |
| `max_markers` | `int` | `60` | Maximum number of markers to select |
| `min_hmm_score` | `float` | `20.0` | Minimum HMM bit score for a valid hit; Phase 0 adaptive parameters can override this (e.g., 30.0 for high-quality datasets, 15.0 for MAG-heavy datasets) |
| `strategy` | `SelectionStrategy` | `GREEDY` | Selection algorithm |
| `marker_hmm_dir` | `str` | `""` | Path to TIGRFAM/Pfam per-marker HMM directory (left empty → auto-discovered via `marker_db_source`; only used when `marker_mode="hmm"`) |
| `marker_mode` | `str` | `"gtdb_tk"` | Marker discovery mode: `gtdb_tk`=consume GTDB-TK-extracted ar53/bac120 per-marker raw sequences (requires `--gtdb-markers-dir`); `hmm`=scan input genomes with a TIGRFAM/Pfam per-marker HMM library (requires `--marker-hmm-dir` or auto-discovery). `denovo` is **not** supported |
| `marker_db_source` | `str` | `"auto"` | HMM library auto-discovery strategy when `marker_mode="hmm"` and `marker_hmm_dir` is empty: `auto`=discover under `db/gtdb_markers/{ar53,bac120}` by detected input domain; `cog`/`gtdb`=only the specified source; `none`=no auto-discovery (explicit `--marker-hmm-dir` required) |
| `gtdb_markers_dir` | `str` | `""` | Required when `marker_mode="gtdb_tk"`: GTDB-TK per-marker FASTA directory (one `.faa`/`.fasta` file per conserved protein; filename stem = marker id; GTDB-TK output is **unaligned/untrimmed** raw sequences, the software automatically runs MAFFT alignment + trimal trimming before tree building) |
| `species_tree` | `str` | `""` | Alias for the `--species-tree` CLI option (`--tree` is also accepted as a legacy alias). Path to a reference species tree (Newick) used in the HGT Phylogenetic step for RF/quartet incongruence screening. In `marker_mode="gtdb_tk"` this is typically the GTDB-TK concatenated species tree. **Recommended** to instead provide `--taxonomy-table`, letting the Phylogenetic step switch to MAD rooting + monophyly proportion screening (see §3.3 and the examples in §8). When neither `--species-tree`/`--tree` nor `--taxonomy-table` is provided, the step is skipped and every marker is recorded as `UNKNOWN` (unscreened: kept, but outside the risk grading). Provide a reference to enable screening. |
| `quality_weighted` | `bool` | `True` | Weight occupancy by genome quality |

### 3.3 HGTConfig

Controls Phase 2 HGT screening.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enable_phylogenetic` | `bool` | `True` | Enable the Phylogenetic step |
| `phylogenetic_threshold` | `float` | `0.5` | Phylogenetic-step risk threshold |
| `level_thresholds` | `Dict` | `{level1_max: 0.25, level2_max: 0.60}` | Level boundaries (normal mode). In far-distance relaxation mode `level2_max` is overridden by `level2_max_far`. |
| `adaptive_far_thresholds` | `bool` | `False` | **Far-distance adaptive relaxation** (**default off**, opt-in). When the input spans multiple high-rank taxa (≥ 2 orders, or genera ratio ≥ 0.70), phylogenetic HGT risk saturates and would drop every marker as Level 3, leaving the pipeline with no markers. When enabled, MarkerFinder relaxes `level2_max` to `level2_max_far` (default 0.95) to retain enough markers for tree building. CLI switches: `--hgt-adaptive-thresholds` (enable) / `--no-hgt-adaptive-thresholds` (explicit disable); config-file key: `hgt_adaptive_thresholds: true`. Note: `adaptive_far_thresholds` is only the HGTConfig dataclass field name and is **not** accepted as a config-file key. |
| `far_distance_genera_ratio` | `float` | `0.70` | Lower bound on genera ratio to trigger far-distance relaxation (`n_genera / n_genomes`). |
| `far_distance_min_orders` | `int` | `2` | Lower bound on number of orders to trigger far-distance relaxation. |
| `level2_max_far` | `float` | `0.95` | `level2_max` used in far-distance relaxation mode. |
| `monophyly_rank` | `str` | `"auto"` | Rank at which the gene-tree monophyly proportion is measured in the Phylogenetic step when `--taxonomy-table` is used. `auto` derives it from the gene tree's taxonomic scope (one rank below the scope; `genus` when the tips are cohesive at no rank). A named rank is measured exactly, with the scope recorded only. Choices: `auto`/`domain`/`phylum`/`class`/`order`/`family`/`genus`/`species`. |
| `monophyly_threshold` | `float` | `0.5` | Lower bound on monophyly proportion: below this value the tree is treated as phylogenetically incongruent (HGT-prone). Phylogenetic-step risk = `1 - monophyly_proportion`. |

### 3.4 PhylogeneticConfig

Controls Phase 3 phylogenetic inference.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `cpus` | `int` | `1` | Thread count for tree builders |
| `ufboot_replicates` | `int` | `1000` | UFBOOT replicates for IQ-TREE3 |
| `fast_mode` | `bool` | `False` | Use IQ-TREE3 `-fast` mode |
| `use_fasttree` | `bool` | `True` | Use FastTree2 instead of IQ-TREE3 (legacy alias, equivalent to `gene_tree_builder: fasttree`) |
| `gene_tree_builder` | `str` | `"fasttree"` | Coalescent gene-tree builder: `"fasttree"` (FastTree2 with WAG, default) or `"iqtree"` (IQ-TREE3, optional) |
| `coalescent_mode` | `str` | `"post-filter"` | `off`, `post-filter`, `always` |
| `iqtree_timeout` | `int` | `3600` | IQ-TREE3 timeout (seconds) |
| `astral_timeout` | `int` | `1800` | ASTRAL timeout (seconds) |

### 3.5 MAGConfig

Controls Phase 0 MAG preprocessing.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `checkm_results` | `str` | `None` | Pre-computed CheckM results file |
| `skip_checkm` | `bool` | `False` | Force skip CheckM and use default quality estimates |
| `min_completeness` | `float` | `50.0` | Minimum completeness for inclusion |
| `min_marker_coverage` | `float` | `0.3` | Minimum marker coverage for species (**legacy field**; use `--min-occupancy` on the CLI — the deprecated alias `--min-marker-coverage` still works and feeds marker-selection occupancy) |
| `quality_weighted` | `bool` | `True` | Weight analysis by genome quality |
| `max_markers` | `int` | `150` | Maximum markers in MAG adaptive mode |

### 3.6 OrthologConfig

Controls multi-copy gene ortholog resolution.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `bbh_evalue` | `float` | `1e-10` | BBH E-value threshold |
| `bbh_identity` | `float` | `30.0` | BBH sequence identity threshold (%) |
| `bbh_coverage` | `float` | `0.6` | BBH coverage threshold |
| `graph_clustering` | `bool` | `True` | Enable graph clustering strategy |
| `use_diamond` | `bool` | `True` | Use DIAMOND instead of BLAST+ |

### 3.7 ReportConfig

Controls report generation.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | `str` | `"./output"` | Output directory |
| `output_prefix` | `str` | `"markerfinder"` | Output file prefix |
| `report_format` | `str` | `"html"` | Report format (only `html` is supported) |

### 3.8 AlignerConfig

Controls alignment operations.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `cpus` | `int` | `1` | Thread count for alignment |

> AlignerConfig / OrthologConfig sub-configs inherit `tmp_dir` from the main `PipelineConfig.tmp_dir`, defaulting to per-run isolated temporary directory; not configured separately.

### 3.9 TaxonomyConfig

Controls taxonomy name parsing.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `taxonomy_table` | `Optional[str]` | `None` | External taxonomy table path |
| `taxonomy_format` | `str` | `"table"` | Taxonomy format: `table`, `embedded` |
| `taxonomy_source_priority` | `str` | `"table"` | Source priority: `table`, `embedded` |
| `taxonomy_delimiter_mode` | `str` | `"reverse"` | Parse mode: `reverse`, `greedy`, `segment` |
| `table_sep` | `Optional[str]` | `None` | Table separator (`None` = auto-detect) |
| `ignore_malformed` | `bool` | `False` | Skip malformed rows instead of terminating (`True`=skip, `False`=terminate) |
| `taxonomy_levels` | `Optional[str]` | `None` | Extend the built-in ranks with custom ones as `level:prefix` pairs (e.g. `kingdom:k__` for Format B, `kingdom:_k_` for Format A). The rank is parsed out of the taxonomy table into the mapping instead of being dropped; the monophyly screen's rank ladder stays `domain..species` |

---

### 3.10 Self-Check Layer: Additional Knobs (Remaining Config Classes)

| Field | Owner | Default | Notes |
|------|------|--------|------|
| `scan_stability_min` | `ReportConfig` | `0.6` | Lower bound on the Jaccard similarity of the selected marker set across threshold levels. Below it the run warns that the conclusion is recorded as `inconclusive`. Set through `--scan-stability-min`, which is really passed down to the scanner. |
| `require_evidence_coverage` | `ReportConfig` | `0.5` | When HGT evidence coverage (`n_actually_measured / n_markers`) falls below this value, `pipeline_summary.txt` prints a prominent warning at the top. Warning only: it does not change the decision rule. |
| `strict_assertions` | `ReportConfig` | `True` | The process ends with exit code 4 as soon as any `severity=fail` output assertion fires. |
| `allow_assertion_failure` | `ReportConfig` | `False` | Explicitly downgrades FAIL to a recorded warning. The bypass must be written into the report. |
| `hgt_scan` | `ReportConfig` | `False` | `--hgt-scan`: writes `{prefix}.threshold_scan.tsv`. Only stored risk values are re-evaluated; **no external tool is re-run**. |
| `composition_screen` | `ReportConfig` | `False` | composition/GC diagnostics switch, derived from `hgt_steps: "phylogenetic,composition"` (or `--hgt-steps`). There is **no CLI or config-file key with this name**. When on, it writes `Phase5_reports/{prefix}.composition.tsv`. |
| `enable_composition_screen` | `HeterogeneityConfig` | `False` | **Equivalent** to `hgt_steps: composition`: either one being true enables the diagnostics (still no CLI/YAML key, in-process only). `composition_warn_threshold` on the same object really drives the outlier test through `PipelineConfig.heterogeneity_config`; a `composition_outlier_metric` other than `rcv` now logs WARN instead of being silently ignored (see the change log). |
| `hgt_mode` | `HGTConfig` | `"risk"` | Marker decision mode: `risk` (the current weighted risk grading), `consistency` ( marker-level cross-framework consistency, with `grade == CONSISTENT` as the inclusion condition) or `hybrid` (**both gates must pass**; defines the joint semantics as the risk gate and the consistency gate being satisfied together). The default `risk` never switches on its own. Before `consistency` or `hybrid` starts, the structured gate of §5.6 must be passed: it checks observable artifacts, so a flag cannot fake it. A missing prerequisite aborts startup and lists what is missing. Once the gate passes, Phase 3 really executes the criteria and writes `excluded_profile.tsv` plus the `consistency_grade` column; the per-marker joint verdict of `hybrid` goes into the evaluation notes (`hgt_evaluation.tsv` and the decision card). |
| `consistency_stringency` | `HGTConfig` | `1` | Consistency strength level 1..5, mapping to the paper's descending thresholds. `provisional`: not yet calibrated against the external benchmark, so it must not become the default criterion beforehand. |
| `hgt_mode` (visible in the report) | `ReportConfig` | `risk` | When set to `consistency`, `pipeline_summary.txt` states that the verdict comes from the cross-framework consistency grading and that the weighted RF+quartet risk score is **used for ranking only**. If the consistency criterion cannot run in that mode (missing tree, or an illegal reference), the report says " demotion NOT IN EFFECT" and notes that the grading is still risk-based. |
| `min_informative_sites` | `HGTConfig` | `0` | Per-marker PIS (parsimony-informative sites) floor; `0` = off. Above 0, a marker whose PIS falls below the value is marked `inconclusive` even when its risk is low, so it never enters a "clean" conclusion, and the decision card records `pis_grade`/`pis_floor`/`pis_measured`. It does not change `ev.level`: dropping the marker remains a publication decision. The threshold awaits calibration. |

**Not wired up (recorded honestly; do not treat as available)**

1. **`--hgt-mode hybrid` is implemented as specified** (the joint semantics defined by: "both gates must pass"). A marker counts as released only when the risk gate (Phase 2 `excluded`, unchanged) and `grade == CONSISTENT` (Phase 3 `screen`) are both satisfied. Markers rejected by the consistency leg go to `excluded_profile.tsv`, and the per-marker verdict notes go to `hgt_evaluation.tsv` and the decision card. **Single-pass semantics**: the trees of this run are built from the risk-passing set, and the consistency leg adjudicates marker-level conclusions and artifacts; it does not pretend a second re-inference took place. The `consistency` mode **is wired into the run path** (same artifacts). When either usable tree is missing the stage reports `NOT EXECUTED` and the column stays `NA`; risk results are never reused to fake a consistency conclusion. The default `risk` mode does not enter this branch.
2. `--cog-category-map <tsv>` is implemented: two columns, `marker_id<TAB>functional_category`, skipping comments, header and single-column rows. A missing file or a file with no usable row logs WARN, the column renders as `NA`, and the report carries a note; no category value is invented.
3. The composition diagnostics artifact (`composition.tsv`) and the PIS floor depend on the `.pis` sidecar, which is only produced on the `gtdb_tk` gene-tree construction path. The `hmm` path has no sidecar, so neither artifact exists there.
4. **`markers.ids` of the must-pass baseline is still an unfilled skeleton** (fill in the 62 ribosomal-protein marker ids of your own reference set). The gate itself is implemented: with empty ids it checks **all** marker gene trees and says so in the log; with ids that do not intersect the dataset it aborts instead of reporting a pass.
5. The multi-universe factorial evaluation is out of scope for this version, per E.3.
6. **Configuration fields that no code reads** (recorded honestly; `tests/unit/test_config_surface_ratchet.py` enforces that any new field is either wired up or listed here): `PipelineConfig.verify_db` and `PipelineConfig.databases` (the database hash check only runs inside `--check`, see; the run path reads neither field); `MAGConfig.min_completeness` (the completeness floor goes through the Phase 0 adaptive occupancy path); `OrthologConfig.use_diamond` and `OrthologConfig.graph_clustering` (DIAMOND/graph-clustering ortholog resolution is not wired up, and Phase 1.5 logs an explicit skip); `HeterogeneityConfig.min_allele_freq`/`consensus_threshold`/`minor_allele_threshold`/`identity_threshold` (reserved for `MAGHeterogeneityHandler`, which is not part of the main pipeline; the main flow keeps only the best hit per marker).
   Note: `HeterogeneityConfig.enable_composition_screen` and `composition_warn_threshold` are no longer dead fields — the first is now equivalent to `hgt_steps: composition` (either one enables the diagnostics), and the second really reaches the outlier test of `composition.tsv` through `PipelineConfig.heterogeneity_config`. A `composition_outlier_metric` other than `rcv` logs WARN stating that only `rcv` is implemented, instead of being silently ignored.
7. **`--check --db-dir <directory>` now hashes the directory you name**. Previously `--check` called the hash check without arguments, so it always verified the repository's default `db/`, which means `--check --db-dir X` drew a conclusion about the **wrong object**. Behaviour is unchanged when `--db-dir` is omitted (still the repository `db/`). Note also that `--check` runs before `--config` is read, so to have it check a specific database pass `--db-dir` on the command line rather than only in a config file.

---

## 4. Pipeline Phases — Detailed Reference

### 4.1 Phase 0: Quality-Aware Preprocessing

**Purpose:** Assess genome quality and compute adaptive parameters.

**Process:**
1. Assess quality source (`Quality source`):
   - **Protein `.faa` input is required**; CheckM cannot run on protein sequences, so it is skipped by default and default quality estimates are used (`Quality source = default_no_nucleotide`).
   - User-supplied `--checkm-results` → precomputed result (`precomputed`).
   - **Explicit `--skip-checkm`** → skip CheckM and use default estimates (`skipped`).
2. Stratify genomes into quality layers (high/medium/low/contaminated)
3. Compute adaptive parameters based on quality distribution

**Quality Layer Standards (MIMAG):**

| Layer | Completeness | Contamination |
|-------|-------------|---------------|
| High | > 90% | < 5% |
| Medium | 70–90% | < 10% |
| Low | < 70% | any |
| Contaminated | any | >= 10% |

**Adaptive Parameter Rules:**

| High-quality Ratio | min_occupancy | max_markers | min_hmm_score |
|-------------------|---------------|-------------|---------------|
| > 70% | 0.75 | 60 | 30.0 |
| 30–70% | 0.55 | 80 | 20.0 |
| < 30% | 0.35 | 150 | 15.0 |

### 4.2 Phase 1: Adaptive Marker Selection

**Purpose:** Dynamically select the optimal set of marker genes for the input dataset.

**Process:**
1. **Marker loading / scanning:**
   - `gtdb_tk` mode (default): read the GTDB-TK ar53/bac120 per-marker FASTA directory given by `--gtdb-markers-dir`; each file is one marker, with the filename stem as marker id.
   - `hmm` mode: run `hmmsearch` for each genome against the TIGRFAM/Pfam per-marker HMM files in `--marker-hmm-dir` (or auto-discovered under `db/gtdb_markers/{ar53,bac120}`).
2. **State Classification:** For each (genome, marker) pair:
   - 0 valid hits (score >= threshold) → `ABSENT`
   - 1 valid hit → `SINGLE_COPY`
   - 2+ valid hits → `MULTI_COPY`
3. **Occupancy Calculation:**
   - `occupancy(marker) = n_single_copy / n_total_genomes`
   - Quality-weighted variant: `Σ(completeness × [state==SINGLE]) / Σ(completeness)`
4. **Marker Selection:** Choose optimal subset using configured strategy

> **Marker discovery mode:** The default `marker_mode="gtdb_tk"` directly consumes GTDB-TK output and performs **no pairwise comparison among the input genomes** — it neither re-searches markers with HMMER nor runs DIAMOND BBH/self-search between genomes. Markers are taken strictly from GTDB-TK's single-copy conserved proteins. GTDB-TK's ar53/bac120 output is **raw, UNALIGNED and UNTRIMMED**, so before tree building each conserved protein is aligned with **MAFFT** and gap-trimmed with **trimal** (required preprocessing for a valid gene tree, not extra marker-discovery alignment). It then builds a tree per conserved protein and judges HGT by inconsistency with the taxonomy or a user-supplied reference species tree, ultimately dropping only high-HGT markers and retaining every other conserved protein.
> `hmm` mode is for runs without a GTDB-TK pre-run: MarkerFinder uses HMMER to scan a TIGRFAM/Pfam HMM directory, identifies markers internally, and extracts the matched protein sequences. Only `gtdb_tk` and `hmm` are accepted for `--marker-mode`; `denovo` is not supported.

**Selection Strategies:**

| Strategy | Algorithm | Best For |
|----------|-----------|----------|
| `GREEDY` | Sort by single-copy occupancy descending, take top N | Standard datasets |
| `INFO_MAX` | Greedy selection with occupancy overlap penalty | Diverse datasets |
| `RATE_BALANCED` | Combined occupancy + multi-copy proportion deviation score | Deep phylogeny |
| `SPARSE_OPTIMIZED` | Include marker if any genome has a hit, then sort by occupancy | MAG datasets |

**Candidate marker pool:** In `gtdb_tk` mode the pool equals the GTDB-TK ar53/bac120 conserved-protein set; in `hmm` mode it is determined by the TIGRFAM/Pfam HMM files in `--marker-hmm-dir`.

### 4.2.5 Phase 1.5: Marker Sequence Extraction and Ortholog Resolution

After the marker set is determined, sequences are extracted from each genome's protein FASTA. The `OrthologResolver` module implements a three-level ortholog-resolution interface (BBH, graph clustering, length consensus), but the main pipeline `_resolve_orthologs` currently skips it: in both `gtdb_tk` and `hmm` modes the input markers are treated as single-copy or best-hit and no additional DIAMOND BBH resolution is performed. Multi-copy cases currently use best-bitscore deduplication. The interface is reserved for future versions.

### 4.3 Phase 2: Phylogenetic HGT Screening

**Purpose:** Identify and classify markers by HGT risk using gene-tree incongruence.

Phase 2 is controlled by `--hgt-steps`, which defaults to `"phylogenetic"`. It detects HGT signals by comparing each marker gene tree against a reference. Two alternative references are supported; if neither is available, the step is skipped and all markers are flagged as HGT-prone (Level 3).

**(a) Reference species tree (`--species-tree`, aliased as `--tree`):** When a reference species tree is provided, each marker gene tree is compared against it for incongruence:

| Metric | Weight | Method |
|--------|--------|--------|
| Normalized RF Distance | 50% | Robinson-Foulds distance / max possible RF |
| Quartet Consistency | 50% | Fraction of quartets with matching topology |

**(b) Taxonomy table (`--taxonomy-table`):** Each marker gene tree is first rooted with **MAD (Minimal Ancestor Deviation, Tria 2017)**, then the **monophyly proportion** is computed. Which rank is measured is decided by `--monophyly-rank`:

- **`auto` (the default):** the scope of a gene tree is the deepest rank at which *all* of its tips share one label, and the screen runs at the rank immediately *below* the scope — a domain-level tree is tested at phylum, a phylum-level tree at class, class→order, order→family, family→genus, genus→species. Example: if every tip of a marker's gene tree belongs to the same phylum, scope=phylum and the monophyly proportion is computed at class. When the tips are not cohesive at any rank (e.g. they span multiple domains), `auto` measures at `genus`.
- **A named rank:** the screen runs at exactly that rank; the detected scope is recorded but never applied.
- **Which taxa are tested:** a taxon enters the proportion only when it has ≥ 2 representative tips on the tree *and* the tips outside it also number ≥ 2. With a single tip on the other side, "{taxon} is monophyletic" is a claim no topology can falsify, so it is reported as not measurable instead of being counted — a denominator that ends up empty skips the screen for that marker, with the reason logged.
- **What counts as monophyletic:** the gene tree must carry the split `{taxon | rest}`, i.e. the taxon *or* its complement is a clade. The verdict is therefore root-invariant and does not depend on where the rooting rule placed the root.
- monophyly proportion = (# taxa whose split the tree carries) / (# testable taxa);
- low monophyly proportion ⇒ the gene tree conflicts with the taxonomic hierarchy ⇒ high HGT propensity. Phylogenetic-step risk = `1 - monophyly_proportion`; below `monophyly_threshold` (default 0.5) the marker is flagged as phylogenetically suspicious (high HGT risk, dropped as Level 3).
- **Rank fallback:** if the rank being measured has no testable taxon (e.g. every genus in a test set has only one species), MarkerFinder walks to other ranks until it finds a testable one, logging the rank actually used (`rank_used` in `Phase5_reports/{prefix}.hgt_evaluation.tsv`). When a *named* rank had to move, the resulting proportion is not comparable with a threshold set for the requested rank: the marker is reported as a cross-rank comparison (UNKNOWN) rather than graded. The screen is skipped only if no rank offers a testable taxon at all.

> Note: monophyly is itself root-invariant; MAD rooting is the uniform rooting method adopted by this layer so that every marker gene tree gets a consistent, unbiased root. Gene-tree tip labels use genome ids, so they must match the first column of `--taxonomy-table`.

**Level Classification:**

| Level | Risk Range | Action |
|------|-----------|--------|
| Level 1 (Clean) | < 0.25 | Priority use |
| Level 2 (Suspect) | 0.25–0.60 | Use with reduced weight |
| Level 3 (Excluded) | >= 0.60 | Remove from analysis |

**Far-distance adaptive relaxation (default off, opt-in):**

When the input spans multiple high-rank taxa (≥ 2 orders, or genera ratio ≥ 0.70), phylogenetic HGT risk saturates and would drop every marker as Level 3, leaving the pipeline with no markers at all. Once enabled, MarkerFinder relaxes `level2_max` to `level2_max_far` (default 0.95) and keeps enough markers for tree building. The CLI flag is `--hgt-adaptive-thresholds`; see §3.3 HGTConfig.

### 4.4 Phase 3: Dual-Mode Phylogenetic Inference

**Mode A: Concatenation (Supermatrix)**

1. Extract marker protein sequences from each genome
2. MAFFT `--auto` alignment per marker
3. trimAl `-automated1` trimming
4. Concatenate into supermatrix with gap (`-`) filling for missing markers
5. Write Nexus partition file
6. IQ-TREE3 with ModelFinder (`-m MFP -B 1000 -bnni`); automatically falls back to FastTree2 if IQ-TREE3 fails

**Output:** `Phase4_trees/markerfinder.species_tree_concat.newick`, `Phase4_alignments/markerfinder.partition.nex`

**Mode B: Coalescent (ASTRAL)**

1. Extract and align marker sequences (same as Mode A)
2. Build individual gene trees with **FastTree2 + WAG** by default, or IQ-TREE3 when `--gene-tree-builder iqtree` is set
3. Filter gene trees (minimum 4 tips)
4. Run ASTRAL-III for coalescent species tree

**Output:**
- `Phase4_trees/markerfinder.species_tree_concat.newick` (Mode A) — always produced by default.
- Coalescent tree file, named by actual method:
  - `Phase4_trees/markerfinder.species_tree_astral.newick` — ASTRAL-III result (only when `--coalescent-mode` is not `off` and ASTRAL-III exits cleanly);
  - ASTRAL-III failure: **no empty tree file is written**; the failure is logged as ERROR and `Phase4_trees/markerfinder.gene_trees.newick` still contains the per-marker gene trees.
- `Phase4_trees/markerfinder.gene_trees.newick` — per-marker gene tree collection.

The actual source is recorded in `Phase5_reports/markerfinder.pipeline_summary.txt` under the `Species tree source` field.

**Conflict Detection**

When both modes are run, MarkerFinder computes the normalized Robinson-Foulds distance and quartet agreement between the concatenation tree and the coalescent tree. The current implementation reports the degree of topological conflict and a recommendation; it does **not** classify individual conflicts as HGT, ILS, or method bias.

**Tree Recommendation Logic:**

| Normalized RF | Recommendation |
|---------------|---------------|
| < 0.1 | High confidence — both methods agree |
| 0.1–0.3 | Medium — prefer the coalescent tree (more robust to ILS) |
| >= 0.3 | Low — strong conflict, examine gene trees |

### 4.5 Phase 4: Report Generation

Generates two types of output:
1. **Static HTML report** — self-contained HTML summary tables, saved to `Phase5_reports/markerfinder.report.html` (Plotly interactive visualizations and PDF are not supported)
2. **Plain text outputs** — TSV, Newick, Nexus files + `Phase5_metadata/run_config.json` parameter snapshot (see Section 7)

---

## 5. External Tool Integration

### 5.1 HMMER (hmmsearch)

- **Purpose:** Search protein sequences against TIGRFAM/Pfam per-marker HMM profiles (`hmm` mode)
- **Command:** `hmmsearch --noali --cpu N --domtblout <out> <hmm_db> <fasta>`
- **Output format:** DOMTABLOUT (domain table)
- **Key fields parsed:** target name, query name, bit score, E-value, alignment coordinates (hmm_from, hmm_to, ali_from, ali_to)

### 5.2 MAFFT

- **Purpose:** Multiple sequence alignment
- **Command:** `mafft --auto --quiet --thread N <input>`
- **Mode:** Auto-selects algorithm based on input size

### 5.3 trimAl

- **Purpose:** Automated alignment trimming
- **Command:** `trimal -automated1 -in <input> -out <output>`
- **Method:** Automated1 heuristic (gap threshold based on alignment statistics)

### 5.4 IQ-TREE3

- **Purpose:** Maximum likelihood tree inference
- **Command:** `iqtree3 -s <aln> -m MFP -B 1000 -bnni -nt N -pre <prefix>`
- **Features:** ModelFinder Plus for model selection, UFBOOT for branch support

### 5.5 ASTRAL-III

- **Purpose:** Coalescent species tree from gene trees
- **Command:** `astral -i <gene_trees> -o <output> -t 2`
- **Output:** Species tree with local posterior probabilities (LPP)

### 5.6 FastTree2

- **Purpose:** Fast approximate ML tree (for gene trees)
- **Command:** `FastTree -quiet -out <output> <alignment>`

### 5.7 DIAMOND

- **Purpose:** Fast protein homology search. Used by the `OrthologResolver` BBH ortholog-resolution interface (currently skipped by the main pipeline; interface reserved for future use)
- **Command:** `diamond blastp --query <faa> --db <db> --evalue <e> --max-target-seqs N --outfmt 6...`

### 5.8 CheckM

- **Purpose:** Genome quality assessment (completeness, contamination)
- **Command:** `checkm lineage_wf -t N -x fna --tablename checkm.tsv <input_dir> <output_dir>`
- **Output:** TSV with Bin Id, Completeness, Contamination columns
- **Input note:** MarkerFinder requires protein input (`.faa`), so CheckM is not run by default because completeness prediction is unreliable without nucleotide sequences. Default quality estimates are used instead (`Quality source = default_no_nucleotide`).
  - Supply a pre-computed quality file via `--checkm-results` to override; `Quality source = precomputed`.
  - Use `--skip-checkm` to explicitly skip CheckM and use default estimates; `Quality source = skipped`.

---

## 6. Database Setup

### 6.1 `hmm` Mode: TIGRFAM/Pfam Per-Marker HMM Directory

In `hmm` mode, MarkerFinder needs a directory of **one HMM file per marker** (e.g. `GTDB-Tk-214-Markers/pfam/individual_hmms/` and `tigrfam/individual_hmms/`). The marker id is the filename stem; suffixes `.HMM` and `.hmm` are accepted. Merged libraries (e.g. `Pfam-A.hmm`) are excluded.

```bash
# Option 1: explicitly provide the per-marker HMM directory
markerfinder -i genomes/ -o output/ --marker-mode hmm \
  --marker-hmm-dir GTDB-Tk-214-Markers/pfam/individual_hmms

# Option 2: assemble db/gtdb_markers/{ar53,bac120} for auto-discovery. The
# profiles are third-party models and are NOT committed to this repository, so
# build the directory from a GTDB-Tk install (or from files you already hold):
python scripts/fetch_marker_db.py --auto
python scripts/fetch_marker_db.py --set bac120 --from-hmm .../gtdbtk_bac120.a.hmm
python scripts/fetch_marker_db.py --set ar53  --from-dir .../ar53_marker_genes
# --record additionally pins the resulting directory hash in
# db/expected_hashes.json, which is what `markerfinder --check` compares against.

# Verify
hmmstat db/gtdb_markers/bac120/*.HMM | head
markerfinder --check                       # reports computed vs expected hash
```

### 6.2 `gtdb_tk` Mode: GTDB-TK Per-Marker FASTA (Required)

This mode does not use `db/` at all: it consumes the per-marker FASTA files that
`gtdbtk align` wrote for **your** genomes, passed with `--gtdb-markers-dir`.

```bash
# Run GTDB-Tk once over the input genomes (it does the alignment and extraction)
gtdbtk align --genome_dir genomes/ --out_dir gtdbtk_out/ --cpus 8

# MarkerFinder reads gtdbtk_out/align/marker_genes/{bac120,ar53}/<marker>.faa
markerfinder -i genomes/ -o output/ --marker-mode gtdb_tk \
  --gtdb-markers-dir gtdbtk_out/align/marker_genes
```

Each marker is one FASTA whose filename stem is the marker id; the GTDB-Tk
version you align with is recorded per run (see 6.4).

### 6.3 gtdb_tk Mode Example (Phylogenetic step via MAD rooting + monophyly proportion)

```bash
# gtdb_tk mode + taxonomy table (Phylogenetic step uses MAD rooting + monophyly proportion, no species tree needed)
markerfinder -i genomes/ -o output/ \
  --marker-mode gtdb_tk \
  --gtdb-markers-dir gtdb_out/ \
  --taxonomy-table gtdb.summary.tsv \
  --monophyly-rank auto \
  --monophyly-threshold 0.5
```

### 6.4 Version Locking

After each run, MarkerFinder writes `Phase5_metadata/run_config.json` to the output directory, recording the full parameter snapshot and database version hashes for reproducibility:

```json
{
  "markerfinder_version": "0.1.0",
  "timestamp": "2025-07-02T10:30:00",
  "run_duration_seconds": 1234.5,
  "parameters": { "mode": "standard", "min_occupancy": 0.75, ... },
  "database_versions": {
    "gtdb_markers": {"name": "gtdb_markers", "path": "db/gtdb_markers/bac120", "hash": "abc123...", "exists": true}
  }
}
```

To reproduce a previous run:

```bash
markerfinder -i genomes/ -o output/ --config previous_run/Phase5_metadata/run_config.json
```

---

## 7. Output Files — Complete Reference

The output directory is organized into phase-prefixed subdirectories to reflect the pipeline steps:

- `Phase5_reports/` — HTML/text reports and TSV summaries
- `Phase4_trees/` — species trees and gene tree collections
- `Phase4_trees/gene_trees/` — single-gene tree cache, one `{marker_id}.nwk` per marker plus its `{marker_id}.nwk.input_sha256` provenance stamp
- `Phase4_alignments/` — Nexus partition files
- `Phase5_metadata/` — `run_config.json`, `context.json`
- `.markerfinder/` — internal state `.pipeline_state.json` for step-by-step execution

### 7.1 HTML Report (`Phase5_reports/markerfinder.report.html`)

Self-contained static HTML file with:

- Execution summary (genomes, markers, runtime)
- HGT level distribution (L1/L2/L3 counts)
- File path notes for all output files

> Only a static HTML report is currently supported; interactive Plotly charts and PDF are not available.

### 7.2 Marker Summary (`Phase5_reports/markerfinder.marker_summary.tsv`)

```tsv
marker_id	occupancy_score	marker_quality_score	marker_quality_level
COG0049	0.9500	0.8234	level_1
COG0085	0.8800	0.7856	level_1
```

### 7.3 HGT Evaluation (`Phase5_reports/markerfinder.hgt_evaluation.tsv`)

```tsv
marker_id	overall_risk	hgt_risk_level	hgt_evidence_confidence
COG0049	0.1200	level_1	high
```

### 7.4 Species Tree (`Phase4_trees/markerfinder.species_tree_concat.newick`)

Concatenation species tree in standard Newick format with branch support values.

### 7.5 Partition File (`Phase4_alignments/markerfinder.partition.nex`)

```nexus
#nexus
begin sets;
  charset COG0049 = 1-300;
  charset COG0085 = 301-600;
end;
```

### 7.6 Pipeline Summary (`Phase5_reports/markerfinder.pipeline_summary.txt`)

Human-readable text file with key statistics. Sample content:

```
============================================================
MarkerFinder Pipeline Summary
============================================================

Genomes: 6
Markers: 60
HGT L1: 0 L2: 30 L3: 30
Species tree source: concat
Quality source: default_no_nucleotide
Runtime: 57.8s
```

- `Species tree source`: `concat` (Mode A supermatrix) / `astral` (Mode B ASTRAL-III succeeded) / `none` (ASTRAL-III failed or coalescent not run). The current implementation does **not** generate `consensus` or `first_gene_tree` files on ASTRAL-III failure. See §4.3.
- `Quality source`: `precomputed` (user-supplied via `--checkm-results`) / `default_no_nucleotide` (protein `.faa` input skips CheckM proactively) / `skipped` (`--skip-checkm` explicitly skipped).

### 7.7 Coalescent Species Tree (`Phase4_trees/markerfinder.species_tree_{source}.newick`)

- **Concatenation species tree** `Phase4_trees/markerfinder.species_tree_concat.newick`: Mode A (Supermatrix) output, inferred with IQ-TREE3 by default and automatically falls back to FastTree2 if IQ-TREE3 fails; **always produced** by default.
- **Coalescent species tree** `Phase4_trees/markerfinder.species_tree_astral.newick`: output only when `--coalescent-mode` is not `off` and ASTRAL-III exits cleanly; on failure no empty tree file is written.

### 7.8 Gene Tree Collection (`Phase4_trees/markerfinder.gene_trees.newick`)

Per-marker gene trees in standard multi-newick format (one tree per line), **no** `>gene_id` separators. Generated when coalescent inference is enabled. Single-gene tree cache lives at `Phase4_trees/gene_trees/{marker_id}.nwk`.

### 7.9 Run Configuration Snapshot (`Phase5_metadata/run_config.json`)

JSON-format full parameter snapshot, including all configuration parameters, database version SHA256 hashes (primarily `gtdb_markers` or `marker_hmm_dir`), software version, and timestamp, for reproducibility (see Section 6.4).

Path fields are recorded as the **absolute locations the run actually used**. A
value relative to the output directory cannot be interpreted later without knowing
which directory it was relative to, and reading it against the current working
directory sends a replay somewhere nobody asked for. Feeding the file back through
`--config <output>/Phase5_metadata/run_config.json` therefore restores the marker
directory, taxonomy table and risk bands, restores the *requested* marker budget
and occupancy floor rather than the values Phase 0 adapted them into, and drops the
previous run's temporary directory (scratch belongs to the run doing the replay).
Fields the command line cannot express (`database_versions`, `databases`, …) are
listed in a warning rather than silently ignored.

---

## 8. Advanced Usage

### 8.1 Using Pre-computed CheckM Results / Skipping CheckM

```bash
# Pre-computed CheckM results
markerfinder -i genomes/ -o output/ -t 8 \
    --checkm-results /path/to/quality_report.tsv

# Explicitly skip CheckM (e.g. quick tests on protein .faa input)
markerfinder -i genomes/ -o output/ -t 8 \
    --skip-checkm
```

The CheckM results file should be a TSV with columns: `Bin Id`, `Completeness`, `Contamination`.
When `--skip-checkm` is used, `pipeline_summary.txt` reports `Quality source: skipped`.

**Not-measured declaration.** Every metric this run could not produce because an external tool or dependency is unavailable is listed at the top of `pipeline_summary.txt` — metric name, cause, affected markers — and the same block is written to the log at WARNING level. Those metrics render as `NA` in the tables and never as a neutral 0 or 0.5. `Evidence coverage` only answers *how many* markers were measured; this declaration answers *which* metric is missing and *why*, and the two are not interchangeable.

### 8.2 Custom HGT Thresholds

```bash
markerfinder -i genomes/ -o output/ -t 8 \
    --hgt-threshold 0.30

# Disable far-distance adaptive relaxation (force normal level boundaries 0.25 / 0.60)
markerfinder -i genomes/ -o output/ -t 8 \
    --no-hgt-adaptive-thresholds

# Opt IN to far-distance adaptive relaxation (relaxes level2_max to 0.95; off by default)
markerfinder -i genomes/ -o output/ -t 8 \
    --hgt-adaptive-thresholds
```

> Far-distance relaxation is **disabled by default**: pass `--hgt-adaptive-thresholds` (or set `hgt_adaptive_thresholds: true` in a config file) to let MarkerFinder relax `level2_max` to 0.95 on **far-distance datasets** (≥ 2 orders or genera ratio ≥ 0.70), so that phylogenetic HGT risk does not saturate and drop every marker as Level 3 (see §3.3 HGTConfig).

### 8.3 Large Dataset Optimization

```bash
markerfinder -i genomes/ -o output/ -t 32 \
    --mode standard \
    --gene-tree-builder fasttree \
    --coalescent-mode post-filter
```

Key optimizations:

- `--gene-tree-builder fasttree`: Use FastTree2 + WAG for gene trees (10× faster; default)
- `--coalescent-mode post-filter`: Run coalescent only after HGT filtering

### 8.4 Skipping Coalescent Inference

```bash
markerfinder -i genomes/ -o output/ -t 8 --coalescent-mode off
```

### 8.5 Using Config Files

Load YAML, TOML, or JSON config files via `--config` to reduce repetitive CLI arguments:

```bash
markerfinder -i genomes/ -o output/ --config config.yaml
```

Example config file (`config.yaml`):

```yaml
threads: 8
mode: standard
hgt_steps: "phylogenetic"
hgt_threshold: 0.25
coalescent_mode: post-filter
ufboot: 1000
```

**Priority rules:** CLI args > Config file > Built-in defaults. Explicitly specified CLI arguments always override config file values.

### 8.6 Output File Management

```bash
# Force overwrite existing output
markerfinder -i genomes/ -o output/ --force

# Skip existing files (no overwrite)
markerfinder -i genomes/ -o output/ --no-clobber
```

Default behavior: If the output directory already contains files, MarkerFinder will exit with an error, prompting you to use `--force` or `--no-clobber`.

### 8.7 Temporary Directory Management

MarkerFinder automatically generates a per-run temporary root directory at `markerfinder-{uuid8}` under the system temp dir, auto-cleaned on exit.

```bash
# Keep temporary directory for debugging (inspect .domtblout / .faa / .aln files)
markerfinder -i genomes/ -o output/ --keep-tmp

# Explicit temporary directory (not auto-cleaned)
markerfinder -i genomes/ -o output/ --tmp-dir ./my_tmp
```

### 8.8 Backward Compatibility with Original MarkerFinder

For users migrating from the original MarkerFinder script:

- Use `--mode conservative` for behavior closest to the original fixed marker set
- The original CSCG sets 16, 27, 37 and 38 are included in the candidate COG pool

---

## 9. Troubleshooting

### 9.1 Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `NoMarkerAvailableError` | No markers meet occupancy threshold | Use `--mode mag_adaptive` or `--mode expanded` |
| `hmmsearch not found` | HMMER not installed | `conda install -c bioconda hmmer` |
| `HMM database not found` | Database path not configured | Set `--marker-hmm-dir` or place TIGRFAM/Pfam HMM files under `db/gtdb_markers/{ar53,bac120}`; or use `gtdb_tk` mode (`--gtdb-markers-dir`) |
| `gtdb-markers-dir not provided` | `gtdb_tk` mode missing required argument | Provide the GTDB-TK per-marker FASTA directory |
| `IQ-TREE3 failed` | Insufficient memory or alignment issues | Use `--gene-tree-builder fasttree` or reduce `--ufboot` |
| `ASTRAL failed` | Gene trees have < 4 taxa | Ensure at least 4 genomes in input |
| `No genomes found` | Wrong directory or file extension | Ensure the input directory contains `.faa` protein files (nucleotide-only input is unsupported) |

### 9.2 Performance Issues

| Symptom | Diagnosis | Solution |
|---------|-----------|----------|
| Slow Phase 1 | hmmsearch bottleneck | Increase `--threads`, use SSD for tmp |
| Slow Phase 2 Phylogenetic step | IQ-TREE3 per marker | Use `--gene-tree-builder fasttree` to build gene trees with FastTree2 + WAG, or reduce the marker set |
| Slow Phase 3 | Large supermatrix | Use `--gene-tree-builder fasttree` |
| High memory | Large alignment | Reduce `--ufboot`, use `--gene-tree-builder fasttree` |

### 9.3 Degraded Mode Behavior

When databases are unavailable, MarkerFinder degrades gracefully:

| Missing Resource | Phase 1 | Phase 2 | Phase 3 |
|-----------------|---------|---------|---------|
| HMM database | All ABSENT → error (or use `gtdb_tk` mode / explicit `--marker-hmm-dir`) | — | — |
| CheckM | Default estimates | — | — |

### 9.4 Exit Codes

MarkerFinder uses standardized exit codes for scripting and workflow integration:

| Exit Code | Meaning | Typical Scenario |
|-----------|---------|------------------|
| `0` | Success | Pipeline completed normally |
| `1` | Runtime error | External tool failure, out of memory, data format error |
| `2` | Argument error | Invalid CLI arguments, missing required arguments |
| `3` | Data error | Input file not found, unrecognized format, tree validation failure, taxonomy table contains unparseable rows |
| `4` | **Output assertion failed / must-pass not met** | Any `severity=fail` output assertion fires, or a relation that must hold in the taxonomy must-pass set does not hold. `--allow-assertion-failure` downgrades this to a warning, and the bypass is recorded in the report |
| `5` | **Inconclusive decision** | The pipeline finished but the recommendation layer cannot deliver a resolved tree: when the conflict metrics are unavailable, `recommend_tree` returns `recommended_tree=None` / `confidence="inconclusive"` instead of handing back a tree as usual |
| `130` | User interrupt | User pressed Ctrl+C (SIGINT) |

---

## 10. Algorithm Details

### 10.1 Robinson-Foulds Distance

```
RF(T1, T2) = number of different bipartitions
normalized_RF = RF / (2 × (n - 3))
```

- 0.0: Identical topologies
- 0.3–0.5: Moderate incongruence
- > 0.5: Strong incongruence

### 10.2 Marker Quality Score

```
score = 0.25 × hmm_norm + 0.25 × occupancy + 0.10 × length_norm 
      + 0.20 × informativeness + 0.20 × (1 - hgt_risk)
```

| Score Range | Level |
|-------------|------|
| >= 0.8 | Level 1 (high quality) |
| 0.5–0.8 | Level 2 (medium quality) |
| < 0.5 | Level 3 (low quality) |

---

## Self-Check Artifacts Quick Reference

| Artifact | Contents | Work package |
|------|------|--------|
| `Phase5_reports/{prefix}.mustpass.tsv` | The `summary`/`violation`/`not_checked` row kinds of the taxonomy must-pass gate (written when `--taxonomy-mustpass` is on) | |
| `Phase5_reports/{prefix}.assertions.tsv` | `assertion_id/name/severity/result/detail/provisional` for every assertion | |
| `Phase5_reports/{prefix}.threshold_scan.tsv` | The `(L1,L2,L3,UNKNOWN)` tuple per level + `n_at_boundary` + the cross-level Jaccard matrix | |
| `Phase5_evidence/decision_<n>.json` | Per-marker decision card (machine-readable, `schema_version: 2`) | |
| `Phase5_reports/{prefix}.excluded_profile.tsv` | `marker_id`, `grade`, `stringency`, `reasons`, `pis`, `functional_category`, `category_note` of the excluded markers, plus each leg's `concat_strength`/`concat_quartets` and `coalescent_strength`/`coalescent_quartets` (strength = the consistency strength, quartets = how many four-tip combinations that consistency was computed from; anything unmeasurable is `NA` and never replaced by 0) | |
| Columns appended to the right of `hgt_evaluation.tsv` | `detector/rank_used/n_total/n_mono/rf/rf_state/quartet/quartet_state/…/assertion_ids_fired` (existing column order unchanged, `None` rendered as `NA`) | |
| Columns appended to the right of `marker_summary.tsv` | `pis` / `effective_columns` / `consistency_grade` (column order unchanged, appended right; `consistency_grade` is `NA` when the gate never ran). The two extra columns listed in, `hgt_risk_score` and `info_rank`, are **not implemented** (the risk score lives in `hgt_evaluation.tsv` and is not duplicated) | |
| Top of `pipeline_summary.txt` | Assertion PASS/WARN/FAIL counts, evidence coverage, **the statement of unmeasured metrics**, and the far-mode effectivity statement | |

New exit codes: **4** = `EXIT_ASSERTION_FAILED` (an output value is implausible, or the taxonomy must-pass set did not hold); **5** = `EXIT_INCONCLUSIVE` (the pipeline finished but the recommendation layer could not deliver a resolved tree). `0/1/2/3/130` keep their meaning.

---

## Glossary

This table shares one source with appendix A of.

| Term | Usage in this project | Caveat |
|---|---|---|
| **marker / gene** | marker = the orthologous marker that gets selected; gene = a copy of that marker in a given genome | New fields always use `marker_id` |
| **risk** | `overall_risk ∈ [0,1]`, a weighted combination of the RF and quartet signals | Not a probability; never read as an "HGT probability" |
| **incongruence** | The disagreement itself between a gene tree and the reference (species tree / taxonomy framework), measured by `normalized_rf` and `quartet_agreement`; it is an **observation** | Kept strictly apart from `risk`: risk is the weighted score used for **ranking**, neither a probability nor an HGT verdict. Incongruence may also come from ILS, assembly fragmentation or long-branch attraction, so "incongruent" does not mean "HGT happened" |
| **consistent / concordant** | `consistent` refers specifically to the marker-level two-framework same-side gate; quartet agreement is a dataset-level proportion of concordant topologies | Not interchangeable: they measure different levels |
| **NOT_MEASURABLE** | A missing dependency, tool or data makes measurement impossible in principle | Distinct from NOT_APPLICABLE; the two lead to different user actions |
| **NOT_APPLICABLE** | The component should not exist on this decision path at all (e.g. the monophyly path does not measure RF) | Not an error, but it must not be filled with 0 either |
| **REJECTED** | Measured, but the reference object failed the legality check | The defence against error A — the merged-reference defect of the retracted paper Steenwyk & King 2025, *Science* 390:751-756, doi:10.1126/science.adw9456; see the "Conceptual Origin" section of the README |
| **UNKNOWN** | `MarkerLevel.UNKNOWN`: the HGT layer did not screen the marker, which is kept but excluded from the risk narrative | A different thing from `inconclusive` (recommendation / consistency layer) |
| **inconclusive** | Evidence exists but is too weak to draw a conclusion (recommendation / consistency layer enum) | `MarkerLevel` does not gain this value |
| **Evidence coverage** | `n_actually_measured / n_markers` | Below `--require-evidence-coverage` the report warns at the top |
| **must-fail control** | A negative fixture that proves a check really can turn red | Without it, a check is equivalent to not existing: every negative scan in this suite carries a planted positive control |
| **Reference tree** | The tree used to score a gene tree (species tree / concatenated tree / ASTRAL tree / the tree itself after MAD rooting) | The reference object must itself pass the legality check |
| **far mode** | `adaptive_far_thresholds`: relaxes level2_max from 0.60 to 0.95 on far-distance datasets | When it takes effect, both the report header and the affected marker cards say so |

---

*End of Manual*
