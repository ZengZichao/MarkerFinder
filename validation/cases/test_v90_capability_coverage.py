"""V-90 — the coverage ratchet: every advertised capability has a case.

This file is what makes "the tests cover everything the software claims"
checkable instead of aspirational. The surface is read from the code — the
argparse parser, the subcommand table, the exit-code constants and the product
inventory — and the coverage is read from the ``@pytest.mark.capability``
markers the cases themselves carry, dumped at collection time by
``conftest.pytest_collection_modifyitems``.

Two directions are enforced:

* a capability in the surface with no case fails the suite (a claim nobody
  tested — the defect class this whole rework exists to remove);
* a capability claimed by a case but absent from the surface also fails (an
  option that was renamed or deleted while its test marker stayed behind).

The matrix is exported to ``results/capability_matrix.tsv`` for the test report.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import capabilities  # Noqa: E402 (validation/ on sys.path by this file)


@pytest.mark.capability("workflow:failure-loudness")
def test_the_snapshot_covers_the_whole_collected_suite(request, results_dir):
    """Guard against a narrowed run producing a narrowed matrix.

    ``pytest -k something`` shrinks the collection, and a matrix built from
    that would silently claim that the missing cases were never needed. The
    snapshot records how many items were collected; this case only judges the
    matrix when the suite was collected whole.
    """
    collected_now = len(request.session.items)
    snapshot_size = capabilities.collected_case_count()
    assert snapshot_size, (
        "no collection snapshot found under validation/.work/capability_snapshot"
    )
    if collected_now < snapshot_size:
        pytest.skip(
            f"only {collected_now} of {snapshot_size} collected cases are "
            "running, so the matrix would be incomplete by construction"
        )


@pytest.mark.capability("workflow:provenance-recorded")
def test_every_advertised_capability_has_a_case(request, results_dir,
                                                record_metric):
    collected_now = len(request.session.items)
    if collected_now < capabilities.collected_case_count():
        pytest.skip("partial collection: the matrix is not authoritative")

    missing = capabilities.export_matrix(results_dir / "capability_matrix.tsv")
    covered = capabilities.covered_capabilities()
    required = capabilities.required_capabilities()
    record_metric("v90_coverage", "capabilities_required", len(required))
    record_metric("v90_coverage", "capabilities_covered", len(covered))
    record_metric("v90_coverage", "capabilities_uncovered", sorted(missing))
    assert not missing, (
        "these advertised capabilities have no case:\n  " + "\n  ".join(missing) +
        "\nAdd a case with @pytest.mark.capability(\"<id>\"), or remove the "
        "capability from the surface if the software no longer offers it."
    )


@pytest.mark.capability("workflow:failure-loudness")
def test_no_case_claims_a_capability_that_does_not_exist(request, results_dir):
    collected_now = len(request.session.items)
    if collected_now < capabilities.collected_case_count():
        pytest.skip("partial collection: the matrix is not authoritative")

    required = set(capabilities.required_capabilities())
    covered = capabilities.covered_capabilities()
    undeclared = sorted(set(covered) - required)
    assert not undeclared, (
        "these capability ids are claimed by cases but are not part of the "
        f"advertised surface: {undeclared}\n"
        "Either the code gained an entry point that is not documented as a "
        "capability, or a case still points at an option that was removed."
    )


@pytest.mark.capability("workflow:provenance-recorded")
def test_the_option_surface_itself_is_what_we_think(request):
    """The parser is the documentation of the command line. A change here is a
    change to what the suite must cover, so the surface is pinned: an option
    added without a role, or an option that quietly disappears, fails here
    rather than in a user's run.
    """
    options = capabilities.cli_options()
    dests = [dest for _flags, dest in options]

    # Several dests are deliberately reachable under two names: the documented
    # Deprecated aliases (--tree/--species-tree). Anything outside that
    # Whitelist means two different options write the same field, which is a
    # Silent precedence question.
    known_alias_groups = {"species_tree"}
    by_dest: dict = {}
    for _flags, dest in options:
        by_dest.setdefault(dest, []).append(_flags)
    collisions = {d: f for d, f in by_dest.items()
                  if len(f) > 1 and d not in known_alias_groups}
    assert not collisions, (
        f"two options share one dest without being a documented alias pair: "
        f"{collisions}"
    )
    for _flags, dest in options:
        assert dest in capabilities.OPTION_ROLES, (
            f"option {dest!r} is exposed by the parser but has no recorded role "
            "in validation/capabilities.py: either it is a new capability that "
            "needs a case, or its role (deprecated alias / no-op / terminal "
            "flag) needs writing down"
        )
    assert len(options) >= 60, (
        f"only {len(options)} CLI options found; the documented surface is "
        "larger, so something is missing from the parser or from this check"
    )


@pytest.mark.capability("subcommand:aggregate-run")
def test_subcommand_and_product_inventories_are_non_empty():
    subs = capabilities.subcommands()
    assert len(subs) == 4, subs
    assert capabilities.PRODUCTS, "the product inventory is empty"
    codes = capabilities.exit_codes()
    assert len(codes) >= 6, codes
