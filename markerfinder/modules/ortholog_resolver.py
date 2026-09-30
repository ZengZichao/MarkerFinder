"""Ortholog resolution helpers and the ``OrthologResolver`` class.

This module is a **reserved interface** for future multi-copy marker
resolution. The current main pipeline (``gtdb_tk`` and ``hmm`` modes) keeps
only the best-scoring hit per (genome, marker) and does **not** invoke
``OrthologResolver``; the class and its helper functions are retained so that
future versions can add DIAMOND BBH-based ortholog selection without breaking
the public API or existing tests.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

from markerfinder.config import OrthologConfig
import logging

logger = logging.getLogger(__name__)


def _extract_candidate_info(candidate: Any) -> Tuple[str, str]:
    """从候选序列中提取 (id, seq)，同时支持 dict 和对象。"""
    if isinstance(candidate, dict):
        return str(candidate.get("id", id(candidate))), str(candidate.get("seq", ""))
    return getattr(candidate, "id", str(id(candidate))), str(getattr(candidate, "seq", ""))


def _extract_candidate_id(candidate: Any) -> str:
    """从候选序列中提取 id，同时支持 dict 和对象。"""
    if isinstance(candidate, dict):
        return str(candidate.get("id", id(candidate)))
    return getattr(candidate, "id", str(id(candidate)))


def _run_diamond_blastp_simple(
    query_seq: str,
    subject_seqs: List[Any],
    tmp_dir: str = "",
) -> Dict[str, float]:
    """运行 DIAMOND blastp 进行序列比对，返回 {subject_id: bitscore}。"""
    scores: Dict[str, float] = {}

    if not tmp_dir:
        import tempfile as _tf
        tmp_dir = _tf.gettempdir()

    query_fd, query_path = tempfile.mkstemp(suffix=".fasta", dir=tmp_dir)
    os.close(query_fd)
    os.chmod(query_path, 0o600)
    subject_fd, subject_path = tempfile.mkstemp(suffix=".fasta", dir=tmp_dir)
    os.close(subject_fd)
    os.chmod(subject_path, 0o600)
    db_path = os.path.join(tmp_dir, f"ortholog_db_{os.path.basename(query_path)}")
    out_fd, out_path = tempfile.mkstemp(suffix=".tsv", dir=tmp_dir)
    os.close(out_fd)
    os.chmod(out_path, 0o600)
    logger.debug(f"Created temp files: {query_path}, {subject_path}, {out_path}")

    try:
        query_target = Path(query_path).resolve()
        query_target.write_text(f">query\n{query_seq}\n", encoding="utf-8", newline="\n")

        subject_lines: List[str] = []
        for s in subject_seqs:
            if isinstance(s, dict):
                sid = str(s.get("id", id(s)))
                seq = str(s.get("seq", s))
            else:
                sid = getattr(s, "id", str(id(s)))
                seq = str(getattr(s, "seq", s))
            subject_lines.append(f">{sid}\n{seq}\n")
        subject_target = Path(subject_path).resolve()
        subject_target.write_text("".join(subject_lines), encoding="utf-8", newline="\n")

        subprocess.run(
            ["diamond", "makedb", "--in", subject_path, "-d", db_path, "--quiet"],
            capture_output=True, check=True, timeout=60,
        )

        subprocess.run(
            [
                "diamond", "blastp",
                "--query", query_path,
                "-d", db_path,
                "--outfmt", "6", "qseqid", "sseqid", "bitscore",
                "--out", out_path,
                "--quiet",
            ],
            capture_output=True, check=True, timeout=120,
        )

        if Path(out_path).exists():
            with open(out_path, encoding="utf-8", newline="") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        scores[parts[1]] = float(parts[2])

    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass
    finally:
        for p in [query_path, subject_path, out_path]:
            Path(p).unlink(missing_ok=True)
        for suffix in [".dmnd"]:
            Path(db_path + suffix).unlink(missing_ok=True)

    return scores


class OrthologResolver:
    """Reserved interface for multi-copy ortholog resolution.

    The current main pipeline does not use this class; it is retained so that
    future versions can implement DIAMOND BBH / graph-clustering / length-based
    ortholog selection without changing the public API.
    """

    def __init__(self, config: OrthologConfig):
        self.config = config
        self.bbh_evalue_threshold = config.bbh_evalue
        self.bbh_identity_threshold = config.bbh_identity
        self.bbh_coverage_threshold = config.bbh_coverage

    def resolve_ortholog(
        self,
        cog_id: str,
        query_genome_id: str,
        candidate_sequences: List[Any],
        reference_genome_ids: List[str],
        reference_sequences: Dict[str, List[Any]],
    ) -> Optional[Any]:
        if len(candidate_sequences) <= 1:
            return candidate_sequences[0] if candidate_sequences else None

        bbh_result = self._bbh_verification(
            cog_id, candidate_sequences, reference_genome_ids, reference_sequences
        )
        if bbh_result is not None:
            return bbh_result

        graph_result = self._graph_clustering_selection(
            cog_id, candidate_sequences, reference_sequences
        )
        if graph_result is not None:
            return graph_result

        return self._length_consistency_fallback(candidate_sequences, reference_sequences)

    def _bbh_verification(
        self,
        cog_id: str,
        candidates: List[Any],
        ref_genome_ids: List[str],
        ref_sequences: Dict[str, List[Any]],
    ) -> Optional[Any]:
        """BBH验证: 用DIAMOND做双向最佳比对。"""
        bbh_scores: Dict[str, float] = {}

        for candidate in candidates:
            cand_id, cand_seq = _extract_candidate_info(candidate)
            if not cand_seq:
                continue

            all_ref_seqs = []
            for ref_id in ref_genome_ids:
                all_ref_seqs.extend(ref_sequences.get(ref_id, []))

            if not all_ref_seqs:
                continue

            forward_hits = _run_diamond_blastp_simple(cand_seq, all_ref_seqs, self.config.tmp_dir)
            if not forward_hits:
                continue

            best_ref_id = max(forward_hits, key=lambda k: forward_hits[k])
            best_ref = next((s for s in all_ref_seqs if _extract_candidate_id(s) == best_ref_id), None)
            if best_ref is None:
                continue

            _, best_ref_seq = _extract_candidate_info(best_ref)
            reverse_hits = _run_diamond_blastp_simple(best_ref_seq, candidates, self.config.tmp_dir)

            if reverse_hits:
                reverse_best = max(reverse_hits, key=lambda k: reverse_hits[k])
                if reverse_best == cand_id:
                    bbh_scores[cand_id] = bbh_scores.get(cand_id, 0) + 1

        if not bbh_scores:
            return None

        best_id = max(bbh_scores, key=lambda k: bbh_scores[k])
        for c in candidates:
            if _extract_candidate_id(c) == best_id:
                return c

        return None

    def _graph_clustering_selection(
        self,
        cog_id: str,
        candidates: List[Any],
        ref_sequences: Dict[str, List[Any]],
    ) -> Optional[Any]:
        """图聚类选择: 选择与参考序列最相似的候选。"""
        if not candidates:
            return None

        all_ref_seqs = []
        for seqs in ref_sequences.values():
            all_ref_seqs.extend(seqs)

        if not all_ref_seqs:
            return candidates[0]

        best_candidate = candidates[0]
        best_score = 0.0

        for candidate in candidates:
            _, cand_seq = _extract_candidate_info(candidate)
            if not cand_seq:
                continue
            hits = _run_diamond_blastp_simple(cand_seq, all_ref_seqs, self.config.tmp_dir)
            max_score = max(hits.values()) if hits else 0
            if max_score > best_score:
                best_score = max_score
                best_candidate = candidate

        return best_candidate

    def _length_consistency_fallback(
        self,
        candidates: List[Any],
        ref_sequences: Dict[str, List[Any]],
    ) -> Optional[Any]:
        if not candidates:
            return None

        ref_lengths = []
        for seqs in ref_sequences.values():
            for s in seqs:
                _, seq_str = _extract_candidate_info(s)
                ref_lengths.append(len(seq_str))

        if not ref_lengths:
            return max(candidates, key=lambda c: len(_extract_candidate_info(c)[1]))

        ref_lengths.sort()
        median_length = ref_lengths[len(ref_lengths) // 2]

        return min(candidates, key=lambda c: abs(len(_extract_candidate_info(c)[1]) - median_length))
