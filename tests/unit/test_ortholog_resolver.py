"""Unit tests for the reserved OrthologResolver interface."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from markerfinder.config import OrthologConfig
from markerfinder.modules.ortholog_resolver import (
    OrthologResolver,
    _extract_candidate_id,
    _extract_candidate_info,
    _run_diamond_blastp_simple,
)


class TestCandidateHelpers:
    def test_extract_candidate_info_dict(self):
        assert _extract_candidate_info({"id": "a", "seq": "ACGT"}) == ("a", "ACGT")

    def test_extract_candidate_info_object(self):
        class Obj:
            id = "b"
            seq = "TGCA"
        assert _extract_candidate_info(Obj()) == ("b", "TGCA")

    def test_extract_candidate_id_defaults_to_object_id(self):
        class Obj:
            id = "c"
        assert _extract_candidate_id(Obj()) == "c"


class TestDiamondBlastpSimple:
    def _mock_run(self, out_lines):
        def runner(cmd, **kwargs):
            # Diamond makedb
            if "makedb" in cmd:
                return MagicMock(returncode=0)
            # Diamond blastp -- write TSV output
            out_path = cmd[cmd.index("--out") + 1]
            Path(out_path).resolve().write_text(
                "".join(line + "\n" for line in out_lines), encoding="utf-8"
            )
            return MagicMock(returncode=0)
        return runner

    def test_returns_scores(self, tmp_path):
        lines = ["query\ts1\t100"]
        with patch("subprocess.run", side_effect=self._mock_run(lines)):
            scores = _run_diamond_blastp_simple("ACGT", [{"id": "s1", "seq": "ACGT"}], tmp_dir=str(tmp_path))
        assert scores == {"s1": 100.0}

    def test_missing_executable_returns_empty(self, tmp_path):
        with patch("subprocess.run", side_effect=FileNotFoundError("diamond")):
            scores = _run_diamond_blastp_simple("ACGT", [{"id": "s1", "seq": "ACGT"}], tmp_dir=str(tmp_path))
        assert scores == {}

    def test_bad_output_lines_ignored(self, tmp_path):
        lines = ["query\ts1", "query\ts2\t90\textra"]
        with patch("subprocess.run", side_effect=self._mock_run(lines)):
            scores = _run_diamond_blastp_simple("ACGT", [{"id": "s2", "seq": "ACGT"}], tmp_dir=str(tmp_path))
        assert scores == {"s2": 90.0}


class TestOrthologResolver:
    def test_resolve_single_candidate(self):
        r = OrthologResolver(OrthologConfig())
        assert r.resolve_ortholog("cog", "g1", [{"id": "a", "seq": "ACGT"}], [], {}) == {"id": "a", "seq": "ACGT"}

    def test_resolve_empty_candidates(self):
        r = OrthologResolver(OrthologConfig())
        assert r.resolve_ortholog("cog", "g1", [], [], {}) is None

    def test_length_fallback(self):
        r = OrthologResolver(OrthologConfig())
        cands = [{"id": "a", "seq": "A" * 100}, {"id": "b", "seq": "A" * 115}]
        refs = {"g2": [{"id": "r", "seq": "A" * 110}]}
        result = r._length_consistency_fallback(cands, refs)
        assert result["id"] == "b"

    def test_graph_clustering_no_refs(self):
        r = OrthologResolver(OrthologConfig())
        cands = [{"id": "a", "seq": "ACGT"}, {"id": "b", "seq": "TGCA"}]
        result = r._graph_clustering_selection("cog", cands, {})
        assert result["id"] == "a"

    def test_bbh_no_refs(self):
        r = OrthologResolver(OrthologConfig())
        cands = [{"id": "a", "seq": "ACGT"}]
        assert r._bbh_verification("cog", cands, [], {}) is None
