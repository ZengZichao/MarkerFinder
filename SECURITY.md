# Security Policy

## Supported versions

MarkerFinder is a young project; only the latest release receives fixes.

| Version | Supported |
|---------|-----------|
| 0.1.x   | yes       |
| < 0.1   | no        |

## Reporting a vulnerability

Please do **not** open a public issue for a security report.

Use GitHub's private vulnerability reporting on this repository
(**Security → Report a vulnerability**), which is enabled. You will be able to
discuss the details privately, and a fix will be credited in the changelog
once it ships. If private reporting is unavailable for some reason, contact
the maintainer directly (see `CITATION.cff` for the author of record) and
include the word "security" in the subject.

## Scope

MarkerFinder is an offline research pipeline. It runs no network services and
holds no secrets, so the realistic exposure is input handling:

- Crafted input files — Newick/NHX trees, FASTA proteomes, taxonomy tables,
  configuration files — that could trigger a crash, an unbounded resource
  use, or a write outside the output directory.
- Shell-out behavior: the pipeline invokes third-party binaries
  (`hmmsearch`, `mafft`, `trimal`, `FastTree`, `iqtree3`, `astral`,
  `datasets`, `checkm`) with paths derived from configuration and file
  names. Reports about argument injection or unexpected execution via those
  paths are in scope.

## What is not in scope

- Memory-safety issues inside the third-party binaries themselves; report
  those upstream (HMMER, MAFFT, EBI, IQ-TREE, ASTRAL, NCBI Datasets).
- Results that are scientifically wrong but produced from valid inputs
  without a safety impact — those are ordinary bugs, please open a normal
  issue.

## What to include

The version (`markerfinder --version`), the exact command line, the input
files (or a minimal reproduction), and the log tail. A `markerfinder --check`
output helps rule out environment problems.
