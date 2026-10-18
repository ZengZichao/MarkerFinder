# AI assistance statement

This file records how a large-language-model (LLM) assistant was used in this
project, so that the statement required by the journal can be reproduced
accurately at submission and checked afterwards.

## What was done with an AI assistant

A commercial large-language-model coding assistant was used to help **revise
parts of the source code and parts of the manuscript text** of MarkerFinder
(v0.1.0). It was not used to produce the research data, to run the analyses, or
to draw the scientific conclusions.

**TODO (author): confirm the version range of the assisted work.** The statement above
covers v0.1.0. Version 0.1.1 adds the `--marker-preset` strategy fix, its regression
test, and a documentation sync; whether any of that was drafted with assistance is a
fact about the project history that only the author can state. Extend the range (and
the checklist) if it was, and leave it at v0.1.0 if it was not.

**TODO (author): name the model and its version here.** The journal's rule is that the
application name and version are reported **in the Methods section** of the paper, not
merely told to the editors; a submission that only says "reported to the editors"
leaves the Methods non-compliant, including in the accepted version. The two values
below and the Methods statement must therefore carry the same text.

* Application name: `________________` — TODO, author only
* Version / release identifier: `________________` — TODO, author only
* Period of use: `________________` — TODO, author only

This file deliberately does not guess those values. Fill them in here **and** in the
manuscript's Methods and the Code Checklist, keeping all three identical; do not
submit while any of them is still blank.

## Record: the 0.1.1 repair round (recorded by the assistant that did it)

Unlike the v0.1.0 work above, the involvement in the v0.1.1 commits is something the
assistant can state from its own evidence rather than reconstruct from memory:

* **Date:** 2026-10-08 (commits `27eefbe`, `8b653a4`, `90a1b00`, `c298cc6`, `38db0f4`
  and `a5fa98f`).
* **Application:** Qoder CN agent CLI, application version **0.4.3** (read from the
  installed bundle's `CFBundleShortVersionString`). The bundle does not surface the
  underlying model identifier, so if the journal wants that too it must come from the
  account's own usage record.
* **Scope:** the `--marker-preset` strategy-dispatch fix and its regression test, the
  config-key warning path, `config.example.yaml`, the database-hash notes and self-test
  wording, `scripts/make_review_package.py` and its new tests,
  `scripts/add_llm_markers.py` and this file, `worked_example/README.md`, and the
  README/MANUAL/TESTING/TEST-REPORT synchronisation. Two later items in the same round
  were fixes of mistakes made during it: `worked_example/provenance_manifest.json` (the
  Supporting Information cited a provenance record that had never existed, so it was
  built from the six downloaded `.faa` files, with counts cross-checked against the
  table before writing) and the runtime assembly of the synthetic e-mail / ORCID /
  absolute-path fixtures in `tests/unit/test_make_review_package.py` (the packager had
  been scrubbing its own test data, which made five cases fail inside the shipped
  package while passing in the repository). Manuscript, SI and figure text in
  the accompanying paper were revised in the same session.
* **Verification:** every code change above was kept only after `python -m pytest tests`
  passed (currently **1039 passed, 1 skipped** — the skip is `NOT EXECUTED: no mypy`),
  and the packager's output was re-scanned for identifiers independently of its own audit.
  The package was additionally **executed**, not just compared: unpacked, `--version`
  reports `markerfinder 0.1.1`, `--check` gives 49 PASS / 0 FAIL, `pytest tests/unit`
  gives 971 passed / 1 skipped, and the documented worked-example command runs to
  completion with the same marker counts, risk statistics and alignments as the
  repository tree (the concatenated tree differs only in the seventh decimal of branch
  lengths, which `validation/cases/test_v15_reproducibility.py` already documents as
  expected).

This is a record of what the tool did, **not** an author attestation: whether any part of
the v0.1.0 work was assisted, and the model name/version required by the journal, remain
the author's to state. The marker-line convention below has still not been applied to any
file, because the file list in `ai-annotation-checklist.md` is empty by design.

## What the authors are responsible for

Every AI-assisted change was **executed and reviewed against the test suite**
before it was kept: the fast suite under `tests/` (`python -m pytest tests`), and,
for changes that touch a user-visible capability, the end-to-end validation suite
under `validation/`. The corresponding author reviewed the changes and takes
responsibility for the correctness of all code and text in this repository,
including the parts drafted with assistance. No result reported in the manuscript
rests on unreviewed assistant output.

## Annotation convention

Affected files carry **one marker line** at the top of the module (immediately
after the module docstring), or immediately above the affected function for a
change that is narrower than a whole module:

```python
# LLM-assisted revision; reviewed against tests (see AI_ASSISTANCE.md)
```

The convention is deliberately a single line rather than a per-hunk log:

* it is greppable (`grep -rn "LLM-assisted revision" markerfinder`, or
  `python scripts/add_llm_markers.py inventory` for the whole tree), so the set of
  assisted files is auditable without relying on anyone's memory;
* it never claims more than "revision assisted" -- authorship, review and
  responsibility stay with the human authors, and the tests are the evidence;
* it does not label a module as "AI-generated". Which files the assistant touched
  is a factual question about the project history, and **only the author can
  answer it**. The list is collected in
  [ai-annotation-checklist.md](ai-annotation-checklist.md) and applied by
  [scripts/add_llm_markers.py](scripts/add_llm_markers.py), which marks *only* the
  files the author lists and refuses to guess.

## Current state of the annotations

As of this writing the list in
[ai-annotation-checklist.md](ai-annotation-checklist.md) is **empty**, and so is the
code: `python scripts/add_llm_markers.py inventory` reports `0 marked file(s)`. The
mechanism exists and is tested (`tests/unit/test_add_llm_markers.py`); the attestation
does not yet. That is an author action, not a tooling gap, and it must be resolved
before the manuscript claims that assisted code is annotated. After filling the list
in, run:

```bash
python scripts/add_llm_markers.py apply --dry-run            # what would be marked
python scripts/add_llm_markers.py apply                      # write the markers
python scripts/add_llm_markers.py check --require-nonempty   # CI gate
```

## Where this statement is echoed

| Place | What it says |
|---|---|
| Manuscript, Methods | AI use, **application name and version, usage period**, review method |
| Code Checklist (editors only) | Which parts were AI-assisted, and how they were verified |
| This file | The convention, the review route, the responsibility statement, and the TODO for the model name |
| [ai-annotation-checklist.md](ai-annotation-checklist.md) | The author's list of assisted targets, and the tool that applies it |

Keep these consistent: the Methods must carry the model name, version and usage period
the journal requires (the TODO above), this file and the Methods must state the same
values, and the checklist must match the marker lines actually present in the code.
A "reported to the editors" formulation is **not** sufficient: the requirement is that
the Methods reports the application name and version.
