# MarkerFinder

**Adaptive HGT-Aware Phylogenomic Pipeline for Prokaryotic Marker Gene Selection**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/ZengZichao/MarkerFinder/releases)
[![Python 3.10～3.12](https://img.shields.io/badge/python-3.10--3.12-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-988%20collected-brightgreen.svg)](#testing)
<!-- The badge count is the live collection total of the fast suite under tests/;
     it is re-measured on every run by tests/unit/test_docs_numbers_are_current.py,
     so a stale number here is a failing test, not a forgotten edit. There is no
     external CI pipeline: the badge reports collected cases, nothing more. -->

[English README](README.EN.md) | [中文说明](README.CN.md) | [Manual (EN)](MANUAL.EN.md) | [Manual (CN)](MANUAL.CN.md)

> Documentation in this repository ships in both languages, with the **English
> version as the primary text**. `tests/unit/test_docs_are_bilingual.py` fails if
> a document exists in only one language, if an "English" file is written in
> Chinese, or if the two versions diverge in shape.

---

## Overview

MarkerFinder is a phylogenomic pipeline for prokaryotic marker gene selection
and species tree inference. It targets three common challenges in existing
methods:

1. **Rigid marker sets** — MarkerFinder dynamically selects optimal markers
   adapted to your dataset
2. **HGT contamination** — phylogenetic-only horizontal gene transfer screening
   removes phylogenetically misleading markers
3. **MAG fragility** — quality-aware adaptive parameters handle incomplete
   metagenome-assembled genomes

## Conceptual Origin

MarkerFinder's integrative (concatenation + coalescence) marker-consistency
idea is methodologically derived from Steenwyk & King (2025), *Science* 390,
751–756, doi:10.1126/science.adw9456 — a paper retracted at the authors' own
request on 2026-02-05 (retraction notice: *Science* 391, 564,
doi:10.1126/science.aef5589). The retraction concerns the paper's biological
conclusion, not its methodological framework; MarkerFinder borrows only the
latter and implements, evaluates and defends it independently. See
[README.EN.md § Conceptual Origin](README.EN.md#conceptual-origin-with-retraction-note)
([中文](README.CN.md#概念来源含撤稿说明)) for the full citation and the
element-by-element mapping into this codebase.

## Quick Start

```bash
# Install (into an environment that also has the external tools)
conda env create -f environment.yml
pip install -e ".[dev]"

# Run
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --mode standard \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes
```

## Documentation

| Document | English | Chinese |
|---|---|---|
| Full README (features, parameters, outputs) | [README.EN.md](README.EN.md) | [README.CN.md](README.CN.md) |
| Manual (formats, schemas, troubleshooting) | [MANUAL.EN.md](MANUAL.EN.md) | [MANUAL.CN.md](MANUAL.CN.md) |
| Reference databases (what is fetched, not shipped) | [db/README.md](db/README.md) | [db/README.md](db/README.md) |
| Test plan & validation methodology | [docs/TESTING.EN.md](docs/TESTING.EN.md) | [docs/TESTING.CN.md](docs/TESTING.CN.md) |
| Test report (measured results) | [docs/TEST-REPORT.EN.md](docs/TEST-REPORT.EN.md) | [docs/TEST-REPORT.CN.md](docs/TEST-REPORT.CN.md) |
| Validation suite (data, cases, archived results) | [validation/README.EN.md](validation/README.EN.md) | [validation/README.CN.md](validation/README.CN.md) |
| External benchmark & positive controls | [tests/benchmark/README.EN.md](tests/benchmark/README.EN.md) | [tests/benchmark/README.CN.md](tests/benchmark/README.CN.md) |

## Testing

The repository carries three test layers. The first two are fast and
self-contained; the third runs the real pipeline on real data and is the
acceptance evidence published with the release.

| Layer | What it does | How to run |
|---|---|---|
| `tests/unit/`, `tests/integration/` | logic, parsing, schemas, wiring; external tools mocked | `pytest tests` |
| `tests/benchmark/` | contract tests for the external-benchmark driver (no GB download) | `pytest tests/benchmark` |
| `validation/` | **full-feature, full-workflow** runs on real genomes, with measured detection metrics | `python validation/run_validation.py --all -n 8` |

`validation/capabilities.py` + `validation/cases/test_v90_capability_coverage.py`
make coverage checkable: every option, subcommand, exit code and output product
the software advertises is enumerated from the code, and any entry without a
case fails the suite. The exported matrix lives in
`validation/results/capability_matrix.tsv`.

## Author

**Zengzichao (曾子超)** — [zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn) ·
[ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

## License

[MIT](LICENSE)
