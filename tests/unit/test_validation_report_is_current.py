"""The published test report must agree with the archived evidence it cites.

``docs/TEST-REPORT.{EN,CN}.md`` is shipped with the release and states measured
numbers. Prose drifts: a case added, a metric recomputed, a gate that stopped
firing. Rather than ask a human to re-read the document, this test re-derives
every headline figure from the artefacts the suite itself wrote —
``validation/results/{junit.xml,metrics.json,capability_matrix.tsv}``, the data
manifest and the capability surface — and requires the report to carry those
exact figures, in BOTH language versions.

A missing artefact is a failure, not a skip: these files are part of the release
bundle, so their absence means the bundle was built without its evidence.
"""

from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "validation" / "results"
DATA = REPO / "validation" / "data"
REPORTS = [REPO / "docs" / "TEST-REPORT.EN.md", REPO / "docs" / "TEST-REPORT.CN.md"]

sys.path.insert(0, str(REPO / "validation"))


@pytest.fixture(scope="module")
def junit() -> ET.Element:
    path = RESULTS / "junit.xml"
    assert path.exists(), (
        f"{path} missing: run `python validation/run_validation.py --all` so the "
        "published evidence exists before it is quoted"
    )
    root = ET.parse(path).getroot()
    return root if root.tag == "testsuite" else root.find("testsuite")


@pytest.fixture(scope="module")
def metrics() -> dict:
    path = RESULTS / "metrics.json"
    assert path.exists(), f"{path} missing: see run_validation.py"
    return json.loads(path.read_text(encoding="utf-8"))


def _count_xfails(suite: ET.Element) -> int:
    """Pytest records an xfail as a ``<skipped>`` whose text says expected."""
    n = 0
    for case in suite.iter("testcase"):
        skip = case.find("skipped")
        if skip is None:
            continue
        blob = " ".join(filter(None, [skip.get("type") or "",
                                      skip.get("message") or "",
                                      skip.text or ""]))
        if "xfail" in blob.lower() or "expected" in blob.lower():
            n += 1
    return n


class TestReportHeadlinesMatchTheArtefacts:
    def test_validation_suite_tally(self, junit):
        tests = int(junit.get("tests"))
        failures = int(junit.get("failures"))
        errors = int(junit.get("errors"))
        skipped = int(junit.get("skipped") or 0)
        xfailed = _count_xfails(junit)
        passed = tests - failures - errors - skipped
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert str(tests) in text, f"{path.name} never mentions {tests} cases"
            assert f"{passed} passed" in text, (
                f"{path.name} does not state {passed} passed (junit says "
                f"{passed} passed, {failures} failed, {errors} errors, "
                f"{xfailed} xfailed out of {tests})"
            )
            assert text.count("0 failed") >= 1, (
                f"{path.name} must state the failure count; junit says {failures}"
            )
            if failures or errors:
                pytest.fail(
                    f"the archived junit records {failures} failures and {errors} "
                    f"errors; the report may not claim a clean run — re-run "
                    "validation and re-write the report from the evidence"
                )

    def test_capability_surface_is_fully_covered(self, metrics):
        from capabilities import required_capabilities

        required = required_capabilities()
        coverage = metrics["v90_coverage"]
        assert coverage["capabilities_uncovered"] == [], coverage
        assert coverage["capabilities_required"] == len(required), (
            f"metrics recorded {coverage['capabilities_required']} capabilities "
            f"but the parser-derived surface now has {len(required)}: the report "
            "and the archived matrix are stale"
        )
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert str(len(required)) in text, (
                f"{path.name} does not mention the {len(required)}-entry surface"
            )

    def test_coverage_by_origin_matches_the_report(self, metrics):
        """The per-source rows of the report's coverage table."""
        from capabilities import (
            PRODUCTS, SUBCOMMAND_CAPABILITIES, WORKFLOW_CAPABILITIES, cli_options,
            exit_codes,
        )

        expected = {
            "options": len({dest for _f, dest in cli_options()}),
            "subcommands": len(SUBCOMMAND_CAPABILITIES),
            "exit codes": len(exit_codes()),
            "products": len(PRODUCTS),
            "workflow properties": len(WORKFLOW_CAPABILITIES),
        }
        total = metrics["v90_coverage"]["capabilities_required"]
        # The row label per language, so the check is on the row's meaning and not
        # On a number that happens to appear somewhere else in the document.
        kinds = {
            REPORTS[0]: {"options": expected["options"],
                         "subcommand": expected["subcommands"],
                         "exit code": expected["exit codes"],
                         "product": expected["products"],
                         "workflow": expected["workflow properties"]},
            REPORTS[1]: {"选项": expected["options"],
                         "子命令": expected["subcommands"],
                         "退出码": expected["exit codes"],
                         "产物": expected["products"],
                         "流程性质": expected["workflow properties"]},
        }
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            for kind, count in kinds[path].items():
                # The coverage table's row for this source must read
                # "| <count> | <count> |" (covered == required).
                pattern = (r"\|[^|\n]*" + re.escape(kind) + r"[^|\n]*\|\s*"
                           + str(count) + r"\s*\|\s*" + str(count) + r"\s*\|")
                assert re.search(pattern, text, flags=re.I), (
                    f"{path.name} lacks a '{kind}' coverage row reading "
                    f"'{count} | {count}' (surface derived from the code)"
                )
            assert str(total) in text

    def test_determinism_and_replay_figures(self, metrics):
        det = metrics["v15_determinism"]
        replay = metrics["v02_replay"]
        assert replay["differing_parameters"] == [], replay
        assert det["files_differing"] == [], det
        # The phrase differs per language; the figure may not.
        zero_word = {REPORTS[0]: "0 differing", REPORTS[1]: "0 个不同"}
        replay_word = {REPORTS[0]: "differing parameters: **[]**",
                       REPORTS[1]: "不一致参数：**[]**"}
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert str(det["files_compared"]) in text, (
                f"{path.name} does not quote the {det['files_compared']} compared "
                "files from v15_determinism"
            )
            assert zero_word[path] in text, (
                f"{path.name} must state the differing-file count as "
                f"{len(det['files_differing'])}"
            )
            assert replay_word[path] in text, (
                f"{path.name} must state the replay's differing-parameter count "
                f"as empty ({replay['differing_parameters']})"
            )

    def test_gene_tree_cache_and_gate_figures(self, metrics):
        cache = metrics["v14_tree_cache"]
        gate = metrics["v11_gate_pass"]["gate_summary"]
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert str(cache["cached_trees"]) in text, (
                f"{path.name} does not quote {cache['cached_trees']} cached trees"
            )
            assert cache["cached_trees"] == cache["stamped_trees"], cache
            assert "groups_tested=8" in text, (
                f"{path.name} does not quote the must-pass gate's tested-group "
                f"count (archived: {gate})"
            )


class TestDetectionMetrics:
    """The planted-transfer benchmark numbers, from the archived metrics."""

    def _metrics(self, metrics):
        path = RESULTS / "chimera_metrics.json"
        assert path.exists(), f"{path} missing: run the validation suite with --all"
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert "directionality" in payload, sorted(payload)
        assert "quality_at_exclusion" in payload, sorted(payload)
        return payload

    def test_directionality_matches(self, metrics):
        payload = self._metrics(metrics)
        direction = payload["directionality"]
        quality = payload["quality_at_exclusion"]
        planted = direction["planted"]
        maxed = [m for m, r in
                 direction["planted_risk_in_chimera_run"].items() if r >= 1.0]
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert f"{len(maxed)} / {len(planted)}" in text, (
                f"{path.name} must state how many planted markers were excluded "
                f"({len(maxed)} of {len(planted)})"
            )
            assert str(direction["n_markers_scored_in_both"]) in text, (
                f"{path.name} omits the scored-marker count "
                f"{direction['n_markers_scored_in_both']}"
            )
            # The bystander claim: no unplanted marker may move.
            assert len(direction["unplanted_markers_that_moved"]) == 0, (
                "the paired runs disagree on unplanted markers; the report's "
                "'0 bystanders' sentence is no longer true"
            )
            assert "Background exclusion rate" in text or "背景剔除率" in text
            rate = quality["background_exclusion_rate"]
            assert str(round(rate, 4)) in text or f"{rate}" in text, (
                f"{path.name} does not quote the background rate {rate}"
            )

    def test_precision_and_recall_are_quoted(self, metrics):
        quality = self._metrics(metrics)["quality_at_exclusion"]
        precision = round(quality["precision"], 3)
        recall = round(quality["recall"], 3)
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert f"{precision}" in text, (
                f"{path.name} must quote precision {precision} (archived "
                f"{quality['precision']})"
            )
            assert (f"{recall}" in text or f"{round(recall, 2)}" in text
                    or "0.667" in text or "0.67" in text), (
                f"{path.name} must quote recall {quality['recall']}"
            )


class TestDataProvenanceFigures:
    def test_manifest_and_provenance_counts(self):
        manifest = (DATA / "MANIFEST.sha256").read_text(encoding="utf-8")
        n_files = len([ln for ln in manifest.splitlines() if ln.strip()])
        header, *rows = (DATA / "PROVENANCE.tsv").read_text(
            encoding="utf-8").splitlines()
        columns, genomes = len(header.split("\t")), len([r for r in rows if r.strip()])
        assert genomes == 12, genomes
        for path in REPORTS:
            text = path.read_text(encoding="utf-8")
            assert f"{columns} " in text or f"{columns}\u00d7" in text \
                or f"{columns} columns" in text or f"{columns} 列" in text, (
                f"{path.name} must state the provenance column count ({columns})"
            )
            assert str(n_files) in text, (
                f"{path.name} does not quote the {n_files}-file manifest"
            )

    def test_report_files_exist_and_are_bilingual(self):
        for path in REPORTS:
            assert path.exists(), path
            assert path.stat().st_size > 4000, path
        en = REPORTS[0].read_text(encoding="utf-8")
        cn = REPORTS[1].read_text(encoding="utf-8")
        # Same section count, so neither version summarises away a chapter.
        assert en.count("\n## ") == cn.count("\n## "), (
            "the English and Chinese reports no longer have the same sections"
        )
