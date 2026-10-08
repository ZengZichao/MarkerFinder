# LLM-assistance annotation checklist

This is the author-completed record behind the annotation convention in
[AI_ASSISTANCE.md](AI_ASSISTANCE.md). It exists because that file says, correctly:

> Which files the assistant touched is a factual question about the project history,
> and **only the author can answer it**.

So this document is a template with an **empty** answer, and the tooling
([`scripts/add_llm_markers.py`](scripts/add_llm_markers.py)) is deliberately built so
that an empty answer is safe: it marks nothing, guesses nothing, and says so.

## Status

| Item | State |
|---|---|
| Annotation convention defined | Done — `AI_ASSISTANCE.md` |
| Mechanism to apply it | Done — `scripts/add_llm_markers.py` (`apply` / `check` / `inventory`) |
| Mechanism tested | Done — `tests/unit/test_add_llm_markers.py` |
| Marker lines actually present in the code | **None** (`python scripts/add_llm_markers.py inventory` reports `0 marked file(s)`) |
| List of assisted files filled in | **NOT DONE — author only** |

The convention is therefore *implemented but not yet populated*. No file in this
repository currently carries the marker line, which means either the assistant
touched nothing that survives in the tree, or the list below has not been filled in.
Only the author can tell which, and the manuscript cannot state an annotated-code
claim until they do.

## How to fill this in

1. Reconstruct the set of files the assistant actually helped revise. Use whatever
   records exist (assistant session logs, your own notes, review history). Do **not**
   use `git log` authorship or file heuristics as a substitute — the point of the
   convention is a human attestation, not an inference.
2. Add one line per target in the fenced block below.
   * `path/to/module.py` — the whole module was assisted. The marker goes immediately
     after the module docstring.
   * `path/to/module.py::some_function` — the change was narrower than a module. The
     marker goes immediately above that function or class, ahead of its decorators.
3. Apply and verify:

```bash
python scripts/add_llm_markers.py apply --dry-run   # show what would be marked
python scripts/add_llm_markers.py apply             # write the marker lines
python scripts/add_llm_markers.py check             # CI gate: markers == this list
python scripts/add_llm_markers.py inventory         # what carries the marker now
```

4. Re-run the test suite (`python -m pytest tests`), because the convention promises
   that every assisted change was reviewed against it.
5. Keep the three disclosure places in step: this list, `AI_ASSISTANCE.md`, and the
   manuscript's Methods statement (see the TODO there for model name/version/period).

## Author-completed list of AI-assisted targets

The block below is the file `add_llm_markers.py` reads. Everything in it that is a
target must be a line of the form `path.py` or `path.py::symbol`; headings, table
rows, bullets, quotes and `#`-commented examples are ignored, so this document stays
readable. It ships **empty on purpose**: an empty list declares nothing, and the tool
refuses to fill it in by guessing.

```text
# ONE TARGET PER LINE. AUTHORED BY THE CORRESPONDING AUTHOR ONLY.
#
# Fill in the files the assistant actually helped revise. Until then this block is
# intentionally empty and the repository makes no annotated-code claim.
#
# Examples (delete the leading '# ' to activate one):
# markerfinder/modules/hgt_filter.py
# markerfinder/pipeline.py::run_pipeline
```

## CI

`python scripts/add_llm_markers.py check` exits

* `0` when every declared target carries the marker and nothing else does (including
  the vacuous case of an empty list, which prints an explicit reminder that it
  declares nothing);
* `1` when a declared marker is missing, or a marker exists on a file this list does
  not declare;
* `2` when the list itself is unusable — absent, a path outside the repository, a
  missing target file, or a named symbol that does not exist.

Use `check --require-nonempty` in the job that must not go green before the author has
filled the list in.

## Reviewer note

The absence of marker lines is disclosed rather than papered over. A reviewer can
reproduce all four commands above in seconds; `inventory` is the one that shows the
current state, and it is the reason this file exists in the package.
