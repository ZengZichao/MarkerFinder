## Summary

<!-- What changes and why, in two or three sentences. 中文描述也可以。 -->

## Checklist

- [ ] `ci-required` is green (unit + integration + benchmark on CPython 3.10–3.14)
- [ ] Documentation updated in both languages (`*.EN.md` and `*.CN.md` twins stay in sync — enforced by `tests/unit/test_docs_are_bilingual.py`)
- [ ] Numbers quoted in the docs still match the measured run (`tests/unit/test_docs_numbers_are_current.py`)
- [ ] For pipeline-behaviour changes: `validation/run_validation.py` re-run on the changed tree
