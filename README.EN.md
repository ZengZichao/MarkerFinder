# MarkerFinder

**Adaptive HGT-Aware Phylogenomic Pipeline for Prokaryotic Marker Gene Selection**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/ZengZichao/MarkerFinder/releases)
[![Python 3.10–3.12](https://img.shields.io/badge/python-3.10--3.12-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-988%20collected-brightgreen.svg)](#testing)

[中文文档](README.CN.md) | [详细手册（EN）](MANUAL.EN.md) | [详细手册（CN）](MANUAL.CN.md)

---

## Overview

MarkerFinder is a phylogenomic pipeline for prokaryotic marker gene selection and species tree inference. It targets three common challenges in existing methods:

1. **Rigid marker sets** — MarkerFinder dynamically selects optimal markers adapted to your dataset
2. **HGT contamination** — Phylogenetic-only horizontal gene transfer screening removes phylogenetically misleading markers
3. **MAG fragility** — Quality-aware adaptive parameters handle incomplete metagenome-assembled genomes

**Biological domain:** Prokaryotic systematics, phylogenomics, metagenome-assembled genome (MAG) classification, and microbial taxonomy.

---

## Pipeline Architecture

```
Genome Input (.faa)
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 0: Quality-Aware Preprocessing               │
│  CheckM assessment → Quality stratification         │
│  → Adaptive parameter computation                    │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 1: GTDB-TK Marker Loading & Selection          │
│  Consume GTDB-TK ar53/bac120 per-marker seqs → occupancy │
│  → Occupancy-based marker subset selection            │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 2: HGT Screening                              │
│  Phylogenetic step: RF distance, quartet consistency, │
│  MAD rooting + monophyly proportion                 │
│  → Level 1 (clean) / Level 2 (suspect) / Level 3 (excl)│
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 3: Dual-Mode Phylogenetic Inference           │
│  Mode A: MAFFT → trimAl → IQ-TREE3 (Supermatrix, default; FastTree2 fallback) │
│  Mode B: MAFFT → trimAl → Gene Trees → ASTRAL-III   │
│  → Conflict detection (RF distance, quartet)         │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 4: Report Generation                          │
│  Static HTML report + plain text outputs             │
│  (TSV, Newick, Nexus) + Phase5_metadata/run_config.json     │
└─────────────────────────────────────────────────────┘
```

---

## Requirements

### Python Dependencies

| Package | Minimum Version | Purpose |
|------------|-----------------|---------------------------------------------------------------------|
| Biopython | >= 1.81 | Sequence parsing |
| PyYAML | >= 6.0 | Config file parsing |
| ete3 | >= 3.1 | Phylogenetic tree ops (Python 3.10–3.12; 3.13 and later need a fallback) |
| tomli | >= 2.0 | TOML config parsing (Python < 3.11) |

> Note: `pandas`, `numpy`, `scipy` and `rich` were declared as dependencies before
> this rework, but the package never imported them. They have been removed from the runtime
> dependencies. `--check` still reports them as information only, without scoring them.

### Environment and Verification (Supported Interpreters)

`pyproject.toml` declares `requires-python = ">=3.10,<3.13"`. **ete3 cannot be imported on
Python 3.13 and later**: its `ete3/webplugin/webapp.py` still does `import cgi`, and the
standard-library `cgi` module was removed in 3.13. The ete3 path for RF/quartet is therefore
only available on 3.10–3.12. On 3.13 and later, `--check` reports "interpreter out of range" and
exits non-zero, while the pure standard-library split-set fallback (`method="splits-python"`)
takes over.

The recommended route is the `environment.yml` at the repository root, which also installs
mafft, trimal, FastTree, IQ-TREE and ASTRAL:

```bash
conda env create -f environment.yml
conda run -n markerfinder markerfinder --check
conda run -n markerfinder pytest -q
```

**When installing with pip instead**: the `ete3` wheel on PyPI does not declare `six`, so
`pip install ete3` alone raises `ModuleNotFoundError: No module named 'six'` on `import ete3`.
Install `pip install ete3 six`. The conda build of `ete3` does not have this problem.

**Measured results inside the supported range** (the same code, one full run per interpreter with
ete3 installed):

| Interpreter | Result |
|---|---|
| CPython 3.10.20 / 3.11.15 / 3.12.13 with ete3 | **0 skipped / 0 failed** overall (the case total lives in the badge above, refreshed with every commit; it is not repeated here); `--check` reports 49 items (48 items on a checkout without the fetched reference databases, where the hash comparison is INFO) with 0 FAIL and exit code 0 |
| CPython 3.12.13 | ete3 prints many `SyntaxWarning: invalid escape sequence` lines on this version. The warnings do not come from the code under test, but they distract the user |
| CPython 3.14.6 (**outside the declared range**) | ete3 cannot be imported: the ete3-dependent modules skip explicitly with `NOT EXECUTED` in the reason, everything else passes and the exit code is 0; `--check` reports "interpreter out of range" and exits non-zero |

The ete3 differential tests (comparing the pure-Python split-set RF/quartet results with ete3 tree
by tree) **only really execute on 3.10–3.12 with ete3**. Without ete3 they skip explicitly with
`NOT EXECUTED` in the reason, which never counts as a pass.

The upper bound of the declared range is itself under test (`test_supported_range_is_earned.py`).
If ete3 ever became importable on 3.13 and later, that test fails and asks the maintainer to
**decide the `<3.13` bound again** rather than keep it from memory.

On Python 3.14 (outside the declared range), running `pytest tests` directly skips the two
ete3-dependent modules during collection with an explicit `NOT EXECUTED` marker rather than
crashing and reporting **no result at all** (exit code 2); 3.11 and 3.12 are brought into the
measured range as well.

### External Tools

| Tool | Minimum Version | Purpose |
|------------|-----------------|----------------------------|
| MAFFT | >= 7.520 | Multiple sequence alignment (per-marker trimming) |
| trimAl | >= 1.4 | Alignment trimming |
| IQ-TREE3 | >= 3.0.0 | Maximum likelihood tree |
| ASTRAL-III | >= 5.7.0 | Coalescent species tree |
| FastTree2 | >= 2.1.11 | Fast gene-tree building |
| DIAMOND | >= 2.1.0 | BBH ortholog resolution (skipped in current main flow; interface reserved) |
| CheckM | >= 1.2 | Genome quality assessment (optional) |

### Databases

| Database | Description | Path | Required? |
|-------------|------------------------------------|-------------------|----------|
| GTDB-TK per-marker FASTA | bac120 (bacteria) + ar53 (archaea) conserved proteins, one file per marker (the file name stem is the marker id) | `--gtdb-markers-dir` | Required for `gtdb_tk` mode — MarkerFinder directly consumes GTDB-TK output |
| TIGRFAM/Pfam per-marker HMMs | Individual `.HMM`/`.hmm` profiles used in `hmm` mode; auto-discovered under `--db-dir`/gtdb_markers/{ar53,bac120} when `--marker-hmm-dir` is omitted. **Not shipped in this repository** — they are third-party models, so assemble them locally with `python scripts/fetch_marker_db.py` (see [`db/README.md`](db/README.md)) | `--marker-hmm-dir` | Required for `hmm` mode (or auto-discovered) |

> The default marker discovery mode is **`gtdb_tk`**; `hmm` mode is provided for runs without a GTDB-TK pre-run. `--marker-mode` only accepts `gtdb_tk` and `hmm`; `denovo` is **not** supported.
> HGT screening uses only the **Phylogenetic** step; **no large external protein database** is required (NCBI nr, 40–80 GB, for example): the pipeline runs fully self-contained.
> Database versions are locked via SHA256 hashes. The run primarily records the `gtdb_markers` directory (in `gtdb_tk` mode) or `marker_hmm_dir` (in `hmm` mode) in `Phase5_metadata/run_config.json` for reproducibility.

---

## Installation

### Option 1: Conda (Recommended)

```bash
# Create environment
conda create -n markerfinder python=3.10 -y
conda activate markerfinder
# Or create the environment directly: conda env create -f environment.yml

# Install external tools
conda install -c bioconda -c conda-forge \
    mafft trimal iqtree astral fasttree \
    diamond checkm-genome ete3

# 1. Use GTDB-TK to extract conserved-protein per-marker FASTA and species tree
#    (run once; MarkerFinder does not invoke GTDB-TK internally)
#    https://ecogenomics.github.io/GTDBTk/

# Install MarkerFinder
pip install markerfinder

# Or install from source
git clone https://github.com/ZengZichao/MarkerFinder.git
cd markerfinder
pip install -e ".[dev]"
```

### Option 2: From Source

```bash
git clone https://github.com/ZengZichao/MarkerFinder.git
cd markerfinder
pip install -e ".[dev]"
```

### Verify Installation

```bash
# Check CLI
markerfinder --version
# Expected: markerfinder 0.1.0
# The version is a constant in the package, so an unpacked source tree, an
# editable install and a built wheel all report the same string.

# Check external tools
mafft --version                # MAFFT
iqtree3 --version              # IQ-TREE3
diamond version                # DIAMOND

# Run unit tests
pytest tests/unit -q
```

---

## Quick Start

### Minimum Runnable Example

```bash
# 1. Prepare input directory with protein FASTA files (.faa)
mkdir -p example/genomes
# Each file should contain protein sequences for one genome
# The example below uses gtdb_tk mode; --gtdb-markers-dir must be prepared with GTDB-Tk

# 2. Run MarkerFinder with standard mode
markerfinder \
    -i example/genomes \
    -o example/output \
    -t 8 \
    --mode standard \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes

# 3. Check results
ls example/output/
# Phase5_reports/markerfinder.report.html       — Static HTML report
# Phase5_reports/markerfinder.marker_summary.tsv — Marker gene summary table
# Phase4_trees/markerfinder.species_tree_concat.newick — Species tree (Newick)
```

### MAG Analysis Example

```bash
markerfinder \
    -i mags/ \
    -o output/ \
    -t 8 \
    --mode mag_adaptive \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --checkm-results checkm_quality_report.tsv
```

### Skipping CheckM (quick tests or when CheckM results are unavailable)

```bash
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --mode standard \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --skip-checkm
```

> `--skip-checkm` bypasses CheckM quality assessment and uses default quality estimates so the pipeline can continue. `pipeline_summary.txt` will report `Quality source: skipped`.

### GTDB-TK Conserved-Protein Mode (gtdb_tk)

No reliance on MarkerFinder's bundled HMM library, and **no genome-internal pairwise
comparison** (no DIAMOND blastp self-search): it consumes the ar53 (archaea) / bac120
(bacteria) per-marker raw sequences that GTDB-TK extracts for your genome set as the
candidate marker set. Note that GTDB-TK's ar53/bac120 output is **unaligned/untrimmed**
raw sequence, so the software first runs MAFFT alignment + trimal trimming per marker
before building the gene tree. The HGT Phylogenetic step can be judged in two ways:

Because MarkerFinder requires protein sequences (`.faa`), HGT screening is performed with the Phylogenetic step only.

- **Recommended:** supply `--taxonomy-table` — each marker gene tree is rooted with
  MAD (Minimal Ancestor Deviation), then the monophyly proportion is computed at the
  rank **one level below the gene tree's automatically inferred taxonomic scope**
  (domain→phylum, phylum→class, class→order, order→family, family→genus, genus→species);
  risk = `1 − monophyly_proportion`, no species tree needed. Naming a rank with
  `--monophyly-rank` measures at exactly that rank instead (default `auto`). A taxon is
  counted only when both sides of its split carry ≥ 2 representatives, and a taxon counts
  as monophyletic when the tree carries the split either way round — the verdict does not
  depend on where MAD put the root.

- Or supply `--species-tree` — compare each gene tree against the GTDB-TK concatenated
  species tree (RF/quartet incongruence ⇒ HGT).

```bash
# 1) Extract conserved-protein sequences + concatenated species tree with GTDB-TK
gtdb-tk align    --genome_dir genomes/ --out_dir gtdbtk_out/ --cpus 8 --extension gz
gtdb-tk classify --genome_dir genomes/ --out_dir gtdbtk_out/ --cpus 8

# 2a) Recommended: supply a taxonomy table; Phylogenetic step uses MAD rooting + monophyly (no species tree)
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --mode standard \
        --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --taxonomy-table gtdbtk_out/classify/gtdbtk.bac120.summary.tsv \
    --monophyly-rank auto \
    --monophyly-threshold 0.5 \
    --hgt-steps phylogenetic

# 2b) Or: supply the GTDB-TK concatenated species tree as the Phylogenetic step reference
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --mode standard \
        --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --species-tree gtdbtk_out/classify/gtdbtk.bac120.classify.tree \
    --hgt-steps phylogenetic
```

> In this mode the marker set is exactly the GTDB-TK conserved-protein set. Because
> GTDB-TK output is unaligned, MAFFT alignment + trimal trimming are run before tree
> building. Only markers judged high-HGT are dropped; everything else is retained. For a
> single domain, supply just the corresponding bac120 (bacteria) or ar53 (archaea)
> directory and taxonomy table / species tree.

### Step-by-Step

MarkerFinder can also run as 4 subcommands (`scan`/`filter`/`infer`/`report`) sequentially, with state auto-passed via the output directory:

```bash
markerfinder scan     -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes
markerfinder filter             -o output/ -t 8
markerfinder infer              -o output/ -t 8 --gene-tree-builder fasttree
markerfinder report             -o output/
```

See [MANUAL.EN.md](MANUAL.EN.md) for details.

---

## Input Format

MarkerFinder accepts **protein FASTA files** (`.faa`, `.fasta`, `.fa`) in the input directory. You can explicitly skip CheckM at any time with `--skip-checkm`.

```
input_dir/
├── genome_1.faa     # Predicted protein sequences for genome 1
├── genome_2.faa     # Predicted protein sequences for genome 2
├── mag_bin1.faa     # MAG protein sequences (auto-detected as MAG)
└── ...
```

**File naming convention:**
- Files containing `mag` or `bin` in the name are automatically classified as MAGs
- Files containing `sag` are classified as single-amplified genomes (SAGs)
- All other files are classified as isolate genomes

**Sequence format:** Standard FASTA format. Each file should contain one protein sequence per gene, with headers in the format `>sequence_id [optional description]`.

---

## Key Parameters

### Basic Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-i, --input` | *required* | Input directory containing `.faa` files |
| `-o, --output` | *required* | Output directory |
| `-t, --threads` | `min(8, CPU count)` | Number of CPU threads |
| `--mode` | `standard` | Analysis mode: `conservative`, `standard`, `expanded`, `mag_adaptive` |
| `--config` | — | YAML/TOML/JSON config file path (CLI args override config file) |

### Tree & Taxonomy Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tree` | — | Alias for `--species-tree`; reference species tree (Newick/Nexus) for the HGT Phylogenetic-step RF/quartet comparison |
| `--sequences` | — | Reference sequence file (FASTA) |
| `--taxonomy-table` | — | External taxonomy table (TSV/CSV) |
| `--taxonomy-format` | `table` | Taxonomy format: `table`, `embedded` |
| `--taxonomy-source-priority` | `table` | Taxonomy source priority: `table`, `embedded` |
| `--taxonomy-delimiter-mode` | `reverse` | Format A parsing strategy: `reverse`, `greedy`, `segment` |
| `--table-sep` | — | Force table separator (auto-detected if not set) |
| `--multi-tree-mode` | `ask` | Multi-tree handling: `ask`, `split`, `first`, `last`, `random` |
| `--strip-annotations` | `false` | Strip NHX annotations from tree |
| `--mol-type` | `auto` | Molecule type: auto-detect or `DNA`, `RNA`, `protein` |
| `--skip-length-check` | `false` | Skip sequence length consistency check |
| `--no-cross-check` | `false` | Skip tree-sequence cross-validation |
| `--ignore-malformed` | `false` | Skip malformed rows/files instead of terminating (default: terminate) |
| `--taxonomy-levels` | — | Extend the built-in ranks with custom ones, as `level:prefix` pairs (e.g. `kingdom:k__` for Format B, `kingdom:_k_` for Format A). The rank is then parsed out of `--taxonomy-table` and reaches the taxonomy mapping downstream, instead of being dropped with an "Unknown level prefix" warning. The monophyly screen's rank ladder stays `domain..species`, so a custom rank is carried, not screened |

### HGT Filtering Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hgt-steps` | `phylogenetic` | Comma-separated list of active HGT screening steps. `phylogenetic` = RF/quartet against a reference species tree, or MAD rooting + monophyly proportion with a taxonomy table. `composition` = RCV/GC composition-bias diagnostics, written as PARALLEL evidence columns in `Phase5_reports/{prefix}.composition.tsv` and **never** merged into the risk score. Default: `phylogenetic` only |
| `--marker-mode` | `gtdb_tk` | Marker discovery mode: `gtdb_tk` (default, consumes GTDB-TK ar53/bac120 per-marker FASTA) or `hmm` (TIGRFAM/Pfam HMM scanning). Only `gtdb_tk` and `hmm` are accepted; `denovo` is not supported |
| `--gtdb-markers-dir` | — | **Required** for `gtdb_tk` mode: GTDB-TK per-marker FASTA directory (one file per conserved protein; filename stem = marker id; GTDB-TK output is unaligned/untrimmed raw sequence, the software runs MAFFT alignment + trimal trimming automatically before tree building) |
| `--marker-hmm-dir` | — | **Required** for `hmm` mode when auto-discovery is disabled: directory of TIGRFAM/Pfam per-marker HMM profiles (one `.HMM`/`.hmm` file per marker, filename stem = marker id). Auto-discovered under `--db-dir` when omitted |
| `--species-tree` | (empty) | GTDB-TK concatenated species tree (Newick) used as the HGT Phylogenetic-step RF/quartet reference. **Recommended** to use `--taxonomy-table` instead so the Phylogenetic step uses MAD rooting + monophyly proportion |
| `--taxonomy-table` | (empty) | Recommended external taxonomy table (TSV/CSV, first column = genome id), used with MAD rooting + monophyly proportion for the Phylogenetic step |
| `--monophyly-rank` | `auto` | Rank at which the Phylogenetic-step monophyly proportion is measured when `--taxonomy-table` is used. `auto` (default) derives it from the gene tree's taxonomic scope (one rank below the scope; `genus` when the tips are not cohesive at any rank). Naming a rank — `domain`/`phylum`/`class`/`order`/`family`/`genus`/`species` — measures at exactly that rank; if no taxon there has an informative split (≥ 2 representatives on both sides) the walk moves upward and the result is reported as a cross-rank comparison (graded UNKNOWN), never as a proportion comparable with `--monophyly-threshold` |
| `--monophyly-threshold` | `0.5` | Lower bound on monophyly proportion; below this value the tree is phylogenetically incongruent (HGT-prone). Phylogenetic-step risk = `1 − monophyly_proportion` |
| `--hgt-threshold` | `0.25` | HGT risk threshold for level classification |
| `--hgt-adaptive-thresholds` | `false` | Enable far-distance adaptive relaxation (**default off**; when on, `level2_max` is relaxed to 0.95 on cross-family / cross-order datasets where phylogenetic HGT risk saturates) |
| `--no-hgt-adaptive-thresholds` | `false` | Explicitly disable far-distance adaptive relaxation (same as the default; placeholder flag) |

### Phylogenetic Inference Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ufboot` | `1000` | UFBOOT replicates for IQ-TREE3 |
| `--fast-tree` | `false` | Use FastTree2 instead of IQ-TREE3 for coalescent gene trees (legacy alias, equivalent to `--gene-tree-builder fasttree`, now also the default) |
| `--gene-tree-builder` | `fasttree` | Coalescent gene-tree builder: `fasttree` (FastTree2 + WAG, fast, default) or `iqtree` (IQ-TREE3, thorough, via `--gene-tree-builder iqtree`) |
| `--coalescent-mode` | `post-filter` | Coalescent inference: `off`, `post-filter`, `always` |

### MAG & Output Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--checkm-results` | — | Pre-computed CheckM quality report |
| `--skip-checkm` | `false` | Force skip CheckM and use default quality estimates |
| `--min-occupancy` | (empty) | Explicit minimum marker-occupancy floor (0-1) for marker selection; overrides the Phase-0 adaptive threshold when set |
| `--min-marker-coverage` | (empty) | **Deprecated** alias of `--min-occupancy` (kept for compatibility) |
| `--db-dir` | `./db` | Database directory |
| `--report-format` | `html` | Report format (only `html` is supported) |
| `--force` | `false` | Force overwrite existing output files |
| `--no-clobber` | `false` | Skip existing files without overwriting |
| `--save-intermediates` | `false` | Copy useful intermediate files to `<output>/Phase4_intermediate/` |
| `--redo` | `false` | In stepwise mode, force re-run an already-completed step |
| `--resume` | `false` | In stepwise mode, auto-run missing prerequisite steps and continue |
| `--tmp-dir` | per-run `markerfinder-{uuid8}` under system temp | Temporary file directory (auto-cleaned by default; an explicit path is not auto-cleaned) |
| `--keep-tmp` | `false` | Keep temporary directory after pipeline completion (for debugging) |
| `--log-file` | — | Log output file path |
| `-v, --verbose` | `0` | Verbose output (stackable: `-v` for INFO, `-vv` for DEBUG) |

### Analysis Modes

| Mode | Strategy | Min Occupancy | Max Markers | Use Case |
|------|----------|---------------|-------------|----------|
| `conservative` | InfoMax | 0.90 | 30 | High-quality isolates, strict filtering |
| `standard` | Greedy | 0.75 | 60 | Mixed isolate datasets |
| `expanded` | RateBalanced | 0.50 | 120 | Diverse datasets, tolerant of missing data |
| `mag_adaptive` | SparseOptimized | 0.40 | 150 | MAG-heavy datasets |

---

## Output Files

### Static HTML Report

- `Phase5_reports/markerfinder.report.html` — Self-contained static HTML report with:
  - Execution summary dashboard
  - Marker summary and HGT risk tables
  - File path notes for all output files
  - Species tree source note

> Only a static HTML report is currently supported; Plotly interactive visualizations and PDF are not available.

### Plain Text Outputs

Output files are organized into phase-prefixed subdirectories to reflect the pipeline steps:

| File | Format | Description |
|------|--------|-------------|
| `Phase5_reports/markerfinder.marker_summary.tsv` | TSV | Selected markers with occupancy, quality, level |
| `Phase5_reports/markerfinder.hgt_evaluation.tsv` | TSV | HGT risk scores and level per marker |
| `Phase5_reports/markerfinder.pipeline_summary.txt` | TXT | Human-readable run summary |
| `Phase5_reports/markerfinder.report.html` | HTML | Static HTML report |
| `Phase4_trees/markerfinder.species_tree_concat.newick` | Newick | Concatenation species tree (IQ-TREE3 by default; FastTree2 fallback) |
| `Phase4_trees/markerfinder.species_tree_astral.newick` | Newick | Coalescent species tree (output only when ASTRAL-III succeeds; on failure no empty tree file is written, only an ERROR log entry) |
| `Phase4_trees/markerfinder.gene_trees.newick` | Newick | Per-marker gene tree collection in standard multi-newick format (one tree per line), generated when coalescent is enabled |
| `Phase4_trees/gene_trees/{marker_id}.nwk` | Newick | Single-gene tree cache (reused between `filter` and `infer`) |
| `Phase4_alignments/markerfinder.partition.nex` | Nexus | Partition file for model selection |
| `Phase5_metadata/run_config.json` | JSON | Full parameter snapshot and database version hashes (`gtdb_markers` or `marker_hmm_dir`) for reproducibility |
| `.markerfinder/.pipeline_state.json` | JSON | Explicit step index for step-by-step execution |
| `Phase4_intermediate/markers/{marker_id}.{faa,aln,aln.trim}` | FASTA | Marker sequences and alignments saved when `--save-intermediates` is enabled |
| `Phase4_intermediate/supermatrix/markerfinder.concat.fasta` | FASTA | Concatenated alignment saved when `--save-intermediates` is enabled |
| `Phase4_intermediate/quality/checkm.tsv` | TSV | CheckM quality results saved when `--save-intermediates` is enabled |

> **Phase prefix legend**: MarkerFinder creates the `Phase4_*` (phylogenetic inference and intermediates) and `Phase5_*` (reports and metadata) subdirectories in the output directory. Phase 0 to Phase 3 are internal logical phases and are **not** emitted as directories.

---

## Resource Requirements

| Dataset Size | CPUs | Memory | Estimated Time |
|-------------|------|--------|----------------|
| 10 genomes, 60 markers | 4 | 4 GB | 15–30 min |
| 50 genomes, 80 markers | 8 | 8 GB | 1–2 hours |
| 200 genomes, 100 markers | 16 | 16 GB | 3–6 hours |
| 500+ genomes | 32+ | 32 GB+ | 12+ hours |

> The Phylogenetic step is the most computationally intensive part of HGT screening. Every candidate marker goes through it.

---

## Testing

```bash
# Run all unit tests
pytest tests/unit -q

# Run with coverage
pytest tests/ --cov=markerfinder --cov-report=html

# Run specific test module
pytest tests/unit/test_marker_selection.py -v

# Environment self-check (49 items, or 48 items without the fetched reference
# databases; exit 0 only when none fails)
markerfinder --check

# Full-feature validation layer: the real pipeline on real genomes
python validation/run_validation.py --all -n 8
```

**Test coverage (measured, not estimated):** the suite has 88 modules and collects
**988** cases — 920 in `tests/unit/`, 17 in `tests/integration/` and 51 in
`tests/benchmark/` (the three directories sum to the badge; both sides are re-checked
by `tests/unit/test_docs_numbers_are_current.py`).
On the supported interpreters (CPython 3.10.20, 3.11.15 and 3.12.13, all with ete3)
the run finishes with **0 failed and 0 skipped**. The badge carries the collected
total (`python3 -m pytest tests --collect-only -q`) and is refreshed with every
commit. The repository ships no CI pipeline definition (there is no
`.github/`), so the badge states a case count, not an external pipeline status.

Coverage spans core configuration, models, marker-selection strategies, the HGT decision engine,
GTDB-TK marker loading, MAD rooting, taxonomy parsing, tree and sequence validation, and report
generation. Three specialised groups sit on top of that:

- behaviour-regression tests for the consolidated modules (`tests/unit/test_module_regressions.py`)
- the must-fail control matrix
- product-level acceptance (`tests/integration/test_products_acceptance.py`), which runs the
  pipeline under mocked external tools and reads the artifacts back

Both `tests/integration/` and `tests/benchmark/` hold real cases: the benchmark directory ships
the fetch, planted-chimera, metric, provenance and reverse-ablation scripts together with their
contract tests. What is still missing is end-to-end numeric correctness against a real
phylogenetics toolchain on real genomes; see the limitations note in the manual.

> **Python version note:** the ete3 dependency is currently incompatible with Python ≥ 3.13 (the standard-library `cgi` module was removed). Please use Python 3.10–3.12.

---

## Reproducibility

MarkerFinder records a complete parameter snapshot in `Phase5_metadata/run_config.json` after each run, including:

- All configuration parameters
- Database version hashes (SHA256)
- Software version
- Timestamp

To reproduce a previous run:

```bash
markerfinder \
    -i genomes/ \
    -o output/ \
    --config previous_run/Phase5_metadata/run_config.json
```

---

## Authors

**Zengzichao (曾子超)** — [ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X),
[zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn)

---

## Citation

If you use MarkerFinder in your research, please cite:

```bibtex
@software{markerfinder2026,
  author       = {Zengzichao},
  title        = {MarkerFinder: Adaptive HGT-Aware Phylogenomic Pipeline},
  year         = {2026},
  version      = {0.1.0},
  url          = {https://github.com/ZengZichao/MarkerFinder},
  orcid        = {0000-0001-6553-970X},
  license      = {MIT}
}
```

---

## Contributing

Contributions are welcome. Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

Please report bugs and feature requests through [GitHub Issues](https://github.com/ZengZichao/MarkerFinder/issues).

---

## Contact

- **Issues & Discussions:** [GitHub Issues](https://github.com/ZengZichao/MarkerFinder/issues)
- **Email:** [zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn)
- **ORCID:** [0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

MarkerFinder builds upon and gratefully acknowledges the following tools and databases:

- [MAFFT](https://mafft.cbrc.jp/) — Multiple sequence alignment tool
- [IQ-TREE3](http://www.iqtree.org/) — Efficient phylogenomic software
- [ASTRAL-III](https://github.com/smirarab/ASTRAL) — Optimal species tree estimation
- [TIGRFAM](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) / [Pfam](https://pfam.xfam.org/) — Protein family and HMM marker libraries
- [GTDB](https://gtdb.ecogenomic.org/) — Genome Taxonomy Database
- [CheckM](https://github.com/Ecogenomics/CheckM) — Genome quality assessment
