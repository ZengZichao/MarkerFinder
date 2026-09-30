#!/usr/bin/env python3
"""Benchmark provenance is enforced, not just recommended.

 asks for "数据不入库；来源 DOI / 许可证 / 获取日期 / 校验和记录在案". Before
``provenance.py`` there was no field anywhere in the shipped skeleton to hold
any of it, so a user who filled in accessions and ran the benchmark ended up
with a number nobody else could re-derive.

These tests are two-sided on purpose: they require the shipped file to be
reported INCOMPLETE (it genuinely is — the order choice is the user's),
and they require a filled-in record to pass, so the gate cannot degenerate into
"always fail" or "always pass".
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1] / "benchmark"
if str(BENCH) not in sys.path:
    sys.path.insert(0, str(BENCH))

import provenance  # Noqa: E402
from provenance import (  # Noqa: E402
    REQUIRED_FIELDS,
    load_accessions,
    missing_provenance,
    require_provenance,
    sha256_file,
)


def _complete():
    data = {"gtdb_release": "r220"}
    data["provenance"] = {
        field: "filled" for field in REQUIRED_FIELDS if field != "gtdb_release"
    }
    data["provenance"]["doi"] = "10.5281/zenodo.1234567"
    data["provenance"]["retrieved_utc"] = "2026-09-20T10:00:00Z"
    data["provenance"]["checksum_sha256"] = "a" * 64
    return data


class TestValidator:
    def test_complete_record_has_nothing_missing(self):
        assert missing_provenance(_complete()) == []

    @pytest.mark.parametrize("field", REQUIRED_FIELDS)
    def test_each_required_field_is_individually_required(self, field):
        data = _complete()
        if field == "gtdb_release":
            data["gtdb_release"] = ""
        else:
            data["provenance"][field] = ""
        assert missing_provenance(data) == [field]

    @pytest.mark.parametrize("junk", ["", "   ", None, "TBD", "TODO",
                                      "N/A", "placeholder",
                                      "REPLACE_WITH_ORDER_NAME", "<fill me>"])
    def test_placeholder_values_count_as_missing(self, junk):
        data = _complete()
        data["provenance"]["doi"] = junk
        assert "doi" in missing_provenance(data)

    def test_missing_provenance_block_reports_every_field(self):
        assert missing_provenance({"gtdb_release": "r53"}) == [
            field for field in REQUIRED_FIELDS if field != "gtdb_release"
        ]

    def test_require_raises_exit_2_and_names_the_field(self, capsys):
        data = _complete()
        data["provenance"]["license"] = ""
        with pytest.raises(SystemExit) as exc:
            require_provenance(data, context="run_benchmark")
        assert exc.value.code == 2
        printed = capsys.readouterr().err
        assert "NOT EXECUTED" in printed and "license" in printed


class TestShippedState:
    def test_shipped_skeleton_is_reported_incomplete(self):
        """Honest status of: the fields exist, the values do not yet."""
        missing = missing_provenance(load_accessions())
        assert missing == [
            "source_registry", "doi", "license", "retrieved_utc",
            "checksum_sha256",
        ], missing
        assert "gtdb_release" not in missing  # Already pinned to r53

    def test_accessions_file_still_declares_the_provenance_block(self):
        data = load_accessions()
        assert isinstance(data.get("provenance"), dict)


class TestChecksumHelper:
    def test_sha256_of_a_known_content(self, tmp_path):
        f = tmp_path / "genome.faa"
        f.write_bytes(b"abc")
        assert sha256_file(f) == hashlib.sha256(b"abc").hexdigest()

    def test_streaming_matches_for_a_large_payload(self, tmp_path):
        payload = b">x\n" + (b"ACGT" * 200000)
        f = tmp_path / "big.faa"
        f.write_bytes(payload)
        assert sha256_file(f) == hashlib.sha256(payload).hexdigest()


class TestRunnerGate:
    def _run(self, monkeypatch, tmp_path, capsys, accessions):
        import run_benchmark

        for tool in run_benchmark.REQUIRED_TOOLS:
            assert tool  # Documents the gate we are bypassing below
        # Pretend every external tool is on PATH: these tests are about the
        # Provenance gate, not about tool discovery.
        monkeypatch.setattr(
            "shutil.which", lambda name: "/usr/bin/" + str(name)
        )
        downloads = tmp_path / "downloads"
        downloads.mkdir()
        (downloads / "manifest.sha256").write_text("x\n", encoding="utf-8")
        monkeypatch.setattr(provenance, "load_accessions", lambda *a, **k: accessions)
        monkeypatch.setattr(
            sys, "argv", ["run_benchmark.py", "--downloads", str(downloads)]
        )
        code = run_benchmark.main()
        return code, capsys.readouterr()

    def test_incomplete_provenance_blocks_the_run(self, monkeypatch, tmp_path,
                                                  capsys):
        code, captured = self._run(
            monkeypatch, tmp_path, capsys,
            {"gtdb_release": "r53", "provenance": {}},
        )
        combined = captured.out + captured.err
        assert code == 2
        assert "provenance incomplete" in combined
        assert "doi" in combined

    def test_complete_provenance_is_echoed_into_the_run_record(self, monkeypatch,
                                                               tmp_path, capsys):
        code, captured = self._run(monkeypatch, tmp_path, capsys, _complete())
        combined = captured.out + captured.err
        assert "provenance incomplete" not in combined
        # The run record must carry what the set came from.
        assert "Provenance recorded" in combined
        assert "10.5281/zenodo.1234567" in combined
        # Still refuses to fabricate metrics while accessions are placeholders.
        assert code == 2
        assert "placeholders" in combined
