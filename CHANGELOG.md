# Changelog

All notable changes to MarkerFinder are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **The bilingual guard now checks content, not only structure.**
  `tests/unit/test_docs_are_bilingual.py` compares, for every EN/CN pair, the
  measured figures (decimals, percentages, data sizes), the command-line flags
  inside fenced blocks, and the ASCII-only inline code spans -- the parts of a
  document that have no translation. Run against the tree as it stood before
  the translation round, 10 of the 18 pair-checks it performs would have
  failed. Each check carries a planted-positive control of its own: a regex
  that matches nothing would otherwise report clean forever.
- `tests/unit/test_config_surface_ratchet.py` requires a reserved knob to be
  named in **both** halves of the MANUAL. It used to read `MANUAL.CN.md` only,
  so the not-wired list could drift in the authoritative English page with the
  suite green.

### Fixed

- Four references that pointed at nothing, or at the wrong thing: the module
  docstring of `tests/unit/test_module_regressions.py` cited
  `MarkerFinder-代码审查.md`, a file that was never in this repository;
  `tests/benchmark/expected/accessions.yaml` sent the reader to a `README.md`
  that exists here only as `.EN.md` / `.CN.md`; and
  `tests/unit/test_hybrid_conjunction.py` asserted the Chinese wording of the
  hybrid spec inside `pipeline.py`, so translating one comment would have
  failed a test that is about the decision being on the page, not about the
  script it is written in.
- The last two defects the audit found: the Chinese halves of `docs/TESTING`
  and `docs/TEST-REPORT` quoted the literal `import cgi` where the English
  names the `cgi` module. Both statements are true; the halves now quote the
  same identifier, which is what the new check compares.

### Changed

- The README badge and suite-count prose are refreshed to the live collection
  (1083 cases: 1015 unit + 17 integration + 51 benchmark), as the single-source
  rule requires, and `docs/TESTING` states the same total. The 0.1.1
  verification record in `AI_ASSISTANCE.md` is now dated to its own commit
  rather than claiming to match whatever the README says today.

## [0.1.1] - 2026-10-18

This release is user-visible in two ways: it widens the supported interpreter
range, and it adds a required runtime dependency. Both were previously wrong in
the manifest, in the code and in the documentation at the same time. It is also
the first release published to PyPI (`pip install markerfinder`), shipped from
the tagged GitHub release through an OIDC trusted-publisher workflow rather
than a long-lived API token.

### Added

- **PyPI distribution.** `.github/workflows/publish.yml` builds the sdist and
  wheel on a tagged release and uploads them to PyPI with
  `pypa/gh-action-pypi-publish` under `permissions: id-token: write` — no
  API token is stored in the repository. The `license` field moved to the SPDX
  expression form (`license = "MIT"`) and the redundant license classifier was
  dropped, clearing the setuptools deprecation from the build log.
- `scripts/make_review_package.py` builds a genuinely anonymous peer-review
  package: identifiers (name, e-mail, ORCID iD, own-repository URLs, BibTeX
  fields) are derived from the repository itself and scrubbed, the build fails
  on any residual identifying-string hit, and the script's own source is
  written so it cannot re-leak what it removes. `tests/unit/test_make_review_package.py`
  covers detection and removal agreement.
- `worked_example/provenance_manifest.json` ships the provenance record the
  manuscript cites, and the worked-example recipe states the correct severity
  for a manifest hash mismatch (integrity signal, not a reproducibility
  blocker).

- **Continuous integration.** `.github/workflows/ci.yml` runs the full fast
  suite on CPython 3.10, 3.11, 3.12, 3.13 and 3.14, with `permissions:
  contents: read`, per-job `timeout-minutes`, and every action pinned to a
  commit SHA. A single `ci-required` job aggregates the matrix so branch
  protection has one stable check name to require. `--check` is now a gate, not
  just a manual command.
- `tests/benchmark/` is covered by CI for the first time. Its 51 contract tests
  counted toward the README badge but had no regression protection.
- `.github/workflows/validation-nightly.yml` runs the acceptance layer on a
  schedule. It was previously a job in the PR gate where it could never have
  passed: that workflow installed only the pip package, while
  `validation/run_validation.py` exits 2 unless `hmmsearch`, `mafft`, `trimal`,
  `FastTree` and `astral` are all on `PATH`.
- `.github/dependabot.yml` for the pip and github-actions ecosystems.
- `CITATION.cff`, including the retraction status of the methodological source
  paper.
- `six` is now declared as a runtime dependency. See Fixed below.
- `pytest-xdist` added to the `dev` extra; the documented
  `run_validation.py -n 8` is xdist's flag and could not run before.
- The nightly acceptance workflow installs `iqtree` (`>=3.0.0`, matching
  `environment.yml`). The pipeline treats `iqtree3` as an optional builder,
  but the v09 support-scale case exercises `--gene-tree-builder iqtree`
  directly and cannot measure either builder's support units without it.

### Fixed

- **`--marker-preset` now dispatches the preset's selection strategy.** The preset
  wrote its strategy to `SelectionConfig.selection_strategy`, a name that is not a
  declared field of the dataclass (the field is `strategy`, and `modules/marker_selection.py`
  reads `self.config.strategy`), so the assignment landed on an ad-hoc instance attribute:
  the log line advertised e.g. `strategy=info_max` while marker selection still ran
  `GREEDY`. `--marker-preset` therefore set only the occupancy floor and the marker budget,
  not the triple the README documents; the archived 0.1.0 behaved that way, and the
  shipped tree now matches the documented semantics. Regression test:
  `tests/unit/test_marker_preset_strategy_regression.py` (`markerfinder/cli/config_build.py:334-341`).
- **A clean install produced an unimportable `ete3`.** ete3 3.1.3 declares no
  requirements at all, yet `ete3/webplugin/webapp.py` does
  `import six.moves.cPickle` and `ete3/__init__.py` star-imports that module.
  Nothing in the declared dependency closure pulled `six` in, so every
  MAD / monophyly measurement either failed or silently degraded to
  "unmeasured". This was invisible on any machine where some unrelated package
  happened to have `six` installed. `six>=1.16` is now declared.
- `markerfinder/utils/gtdb_tk_markers.py` referenced `ft_last_err` inside an
  f-string, but that name was never assigned on any code path. Whenever
  FastTree was missing — the ordinary case — the "FastTree unavailable"
  warning raised `NameError` instead of warning, turning a recoverable
  degradation into a crash.
- `--check` reported the interpreter as out of range on Python 3.13 and later
  and exited non-zero, citing a `requires-python` value the manifest no longer
  contained. The check label is now derived from the constants rather than
  restating a literal.
- The documentation-number guard skipped itself on Python 3.13+ on the premise
  that `ete3` could not import there. That premise stopped being true once the
  `cgi` stand-in landed, so the three tests that keep the README numbers honest
  were silently disabled on every supported interpreter. It now asks whether
  `ete3` is importable.
- `cli/self_test.py` and `cli/validation.py` used `Optional` and `argparse` in
  annotations without importing them, so `typing.get_type_hints()` failed on
  those functions.
- The `mag_named` validation fixture's four genome symlinks dangle in a fresh
  checkout, and the assertion only checked that the directory existed, so a
  missing link surfaced much later as an unparseable-FASTA error from inside the
  pipeline. It now fails with the `--prepare` command to run.
- `test_ete3_is_usable_on_every_supported_interpreter` swallowed the underlying
  exception and reported only "ete3 is installed but unimportable". The actual
  cause is now in the assertion message.
- **`MANIFEST.sha256` recorded stale digests for two shipped must-pass
  tables.** The committed `mustpass/mustpass.yaml` and
  `mustpass/mustpass_domain.yaml` do not match the hashes the manifest listed,
  so the manifest checksum case failed on every clean checkout and could only
  pass where a working tree still held the bytes the manifest was generated
  from. The two entries now record the digests of the files as committed.
- **The acceptance layer now reports an unfetchable dependency as skipped,
  not failed.** The hmm-mode cases and the `database_versions` assertion can
  only hold when the third-party TIGRFAM/Pfam profile library exists
  (`db/gtdb_markers`, `validation/data/hmms/core` — both git-ignored by
  design, see `db/README.md`). On a checkout that never fetched it they now
  skip with that stated reason — a `hmm_library` fixture, a guard in the
  `hmm_dir` fixture, and a `iqtree3`-presence guard on the v09 support-scale
  case — instead of erroring in a way that was indistinguishable from a
  regression. The capability matrix is unaffected: it is built from
  collection-time declarations, not from run outcomes.

### Changed

- **Supported interpreters: CPython 3.10–3.12 → 3.10 and newer, with no upper
  bound.** The ceiling existed because ete3 imported the standard-library `cgi`
  module, which CPython removed in 3.13 (PEP 594).
  `markerfinder/_cgi_compat.py` installs a minimal stand-in when, and only when,
  the real module is absent, so the constraint no longer applies. The absence of
  an upper bound is itself under test.
- `environment.yml` no longer pins `python=3.10`.
- Eight test skip reasons that said "must run in a 3.10-3.12 environment" now
  describe the actual condition, which is a broken `ete3` installation.
- The version is now read from `markerfinder/_version.py` by
  `[tool.setuptools.dynamic]`; it used to be restated as a literal in both
  `pyproject.toml` and `_version.py`.
- README, docs and manuals were reconciled against the manifest, in both
  languages. The English and Chinese documents had drifted apart — the
  bilingual test checked structure, not content, which is why the drift could
  survive it. (Superseded: see `[Unreleased]`.)

## [0.1.0] - 2026-09-30

Initial public release. See the
[release notes](https://github.com/ZengZichao/MarkerFinder/releases/tag/v0.1.0).
