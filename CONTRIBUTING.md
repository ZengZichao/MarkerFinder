# Contributing to MarkerFinder

Thank you for considering a contribution. This document describes how the
repository is wired, so that a change lands green on the first try.

## The two test layers

The repository runs two deliberately separate suites; knowing which one your
change belongs in is most of the work.

| | `tests/` (fast layer) | `validation/` (acceptance layer) |
|---|---|---|
| Runs on | every push and pull request (`ci.yml`) | nightly schedule + manual dispatch (`validation-nightly.yml`) |
| Needs | pure Python only | real genomes (fetched from NCBI) and the external tools on `PATH` |
| Gate | the `ci-required` branch-protection check | never a pull-request gate |

- `tests/unit`, `tests/integration`, `tests/benchmark` — run in seconds:

  ```bash
  python -m pytest tests/unit -q
  python -m pytest tests/integration -q
  python -m pytest tests/benchmark -q
  ```

- `validation/` — 153 cases running the real pipeline end to end:

  ```bash
  python validation/run_validation.py --prepare   # first run: fetches 12 RefSeq proteomes and builds fixtures
  python validation/run_validation.py -n 4        # everything except the slow planted-chimera case
  ```

## Local setup

```bash
conda env create -f environment.yml
conda activate markerfinder
pip install -e ".[dev]"
python -m markerfinder --check
```

The acceptance layer additionally needs `hmmsearch`, `mafft`, `trimal`,
`FastTree`, `astral` and (for the v09 support-scale case) `iqtree3` on
`PATH`; `environment.yml` installs exactly that set.

### The third-party HMM library

The TIGRFAM/Pfam profiles behind `--marker-mode hmm` are **not redistributed**
(see `db/README.md`). On a fresh checkout the hmm-mode validation cases
report themselves as **skipped** rather than failed. To run them for real:

```bash
python3 scripts/fetch_marker_db.py --auto      # assembles db/gtdb_markers from a GTDB-Tk install
python validation/scripts/02_build_marker_sets.py   # builds validation/data/hmms/ and the marker fixtures
```

## Conventions a pull request is held to

1. **Bilingual documentation.** Every document ships as an `*.EN.md` /
   `*.CN.md` pair, English first. `tests/unit/test_docs_are_bilingual.py`
   enforces pairing, parallel heading shape, and that the README index links
   both languages in EN-first order. A docs change is not done until both
   halves are updated.
2. **The changelog.** User-visible changes get an entry under `[Unreleased]`
   in `CHANGELOG.md`.
3. **Honest numbers.** Claims in the README and the validation report are
   checked against the shipped evidence
   (`tests/unit/test_docs_numbers_are_current.py`,
   `tests/unit/test_validation_report_is_current.py`). If your change alters
   a measured quantity, re-run the relevant validation cases and refresh
   `validation/results/` — the committed evidence is the source those tests
   read.
4. **Skips must say why.** A case that cannot be measured in a given
   environment (missing optional binary, unfetched third-party library)
   calls `pytest.skip` with the reason and the command that would make it
   measurable. A bare failure is reserved for a real regression.
5. **Declare what you cover.** Acceptance cases declare the capability they
   exercise with `@pytest.mark.capability(...)`; the coverage matrix is built
   from those declarations, so an undeclared case is invisible to it.

## Commit and pull-request style

- Titles follow the existing pattern: `type(scope): summary`
  (`fix(validation): …`, `feat(cli): …`, `docs: …`).
- One logical change per pull request; the fast CI runs the whole matrix on
  CPython 3.10–3.14, so unrelated churn only widens the review surface.

## Reporting bugs and security issues

Use the issue templates. Security-relevant reports go through the private
channel described in `SECURITY.md` — never a public issue.
