"""Adaptive marker selection: HMMER-based scan of input genomes.

Two operating modes, selected by ``SelectionConfig.marker_mode``:

- ``gtdb_tk`` (see ``utils.gtdb_tk_markers``): consume per-marker aligned
  FASTA files produced by an external ``gtdb-tk align`` run. No HMM/DIAMOND
  searching happens here.
- ``hmm``: given a directory of Pfam/TIGRFAM ``.HMM``/``.hmm`` profiles, run
  ``hmmsearch`` on each input genome, build a per-marker occupancy matrix
  (SINGLE / MULTI / ABSENT), select a high-occupancy marker subset, and
  extract the matched protein subsequence from each genome FASTA.

This module implements both modes. For ``hmm`` the marker set is derived from
the files present in ``marker_hmm_dir`` (one profile per marker, filename stem=id),
not from a hard-coded COG pool.
"""

from __future__ import annotations

import os
import subprocess
import uuid
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from markerfinder.phases import canonical_log_tag
from markerfinder.config import SelectionConfig
from markerfinder.models.genome import GeneState, Genome, OccupancyMatrix
from markerfinder.models.marker import (
    MarkerQualityScore,
    MarkerLevel,
    SelectedMarkerSet,
    SelectionStrategy,
)
from markerfinder.models.pipeline_types import (
    MarkerSelectionResult,
    NoMarkerAvailableError,
)
from markerfinder.exceptions import ExternalToolError
from markerfinder.utils.io import parse_fasta
from markerfinder.utils.concurrency import parallel_map
import logging

logger = logging.getLogger(__name__)


# 每趟 MarkerFinder run 唯一的 8-hex 标记。所有 domtblout / per-run tmp 文件带上此后缀,
# 避免多趟并行跑共享默认 /tmp 时文件名冲突(ENOENT 或覆盖).
_RUN_UUID = uuid.uuid4().hex[:8]


_MERGED_LIB_STEMS = {"Pfam-A", "tigrfam", "Pfam-B"}

def gather_marker_ids(marker_hmm_dir: str) -> List[str]:
    """Derive sorted marker ids from individual.HMM/.hmm files in ``marker_hmm_dir``.

    Only per-marker profiles are considered; top-level merged libraries
    (``Pfam-A.hmm``, ``tigrfam.hmm``) are excluded.
    """
    d = Path(marker_hmm_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"HMM marker directory not found: {marker_hmm_dir}")
    profiles = sorted({
        p for p in d.rglob("*.hmm") if p.stem not in _MERGED_LIB_STEMS
    } | {
        p for p in d.rglob("*.HMM") if p.stem not in _MERGED_LIB_STEMS
    })
    if not profiles:
        raise FileNotFoundError(
            f"No .HMM/.hmm per-marker profile files found under {marker_hmm_dir}. "
            "TIGRFAM/Pfam profiles (one file per marker, e.g. TIGR00006.HMM under individual_hmms/) are required."
        )
    return sorted(set(p.stem for p in profiles))


def hmm_path_for(marker_hmm_dir: str, marker_id: str) -> Optional[str]:
    """Locate the.HMM/.hmm file whose stem equals marker_id."""
    d = Path(marker_hmm_dir)
    for ext in (".HMM", ".hmm", ".HMM.gz", ".hmm.gz"):
        p = d / f"{marker_id}{ext}"
        if p.exists():
            return str(p)
        cand = d / "individual_hmms" / f"{marker_id}{ext}"
        if cand.exists():
            return str(cand)
    # Recursive fallback
    for p in d.rglob("*"):
        if p.is_file() and p.stem == marker_id and p.suffix in (".HMM", ".hmm"):
            return str(p)
    return None


def _run_hmmsearch(
    fasta_path: str,
    hmm_path: str,
    cpus: int = 1,
    tmp_dir: str = "",
) -> Dict[str, List[Dict]]:
    """Run hmmsearch of a single genome against a single HMM profile."""
    result: Dict[str, List[Dict]] = defaultdict(list)

    if not tmp_dir:
        import tempfile
        tmp_dir = tempfile.gettempdir()

    if not Path(hmm_path).exists():
        logger.warning(f"  HMM profile not found: {hmm_path}")
        return result

    stem_f = Path(fasta_path).stem
    stem_h = Path(hmm_path).stem
    # Create the scratch directory before pointing hmmsearch at it. A missing
    # Directory makes ``subprocess`` raise FileNotFoundError too, and the handler
    # Below would then report "hmmsearch not found on PATH" for a tool that is
    # Installed — a wrong diagnosis is worse than a raw error.
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    domtblout = os.path.join(
        tmp_dir, f"hmmsearch_{stem_f}_{stem_h}_{_RUN_UUID}.domtblout")
    cmd = [
        "hmmsearch", "--noali", "--cpu", str(cpus),
        "--domtblout", domtblout, hmm_path, fasta_path,
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=600)
    except FileNotFoundError as e:
        # 工具缺失/失败必须**硬失败**, 绝不能用空结果集冒充"该 marker 在所有基因组
        # 中 ABSENT" —— 否则下游会把"工具挂了"误判为"标记确实缺失", 静默丢失信息.
        raise ExternalToolError(
            "hmmsearch not found on PATH. Install HMMER "
            "(e.g. `conda install -c bioconda hmmer`) before running HMM-mode "
            "marker selection. A missing tool must not be silently treated as "
            "an empty hit set (which would falsely mark every marker ABSENT)."
        ) from e
    except subprocess.CalledProcessError as e:
        raise ExternalToolError(
            f"hmmsearch failed for {fasta_path} vs {stem_h}: "
            f"{getattr(e, 'stderr', '')[:200]}"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            f"hmmsearch timed out for {fasta_path} vs {stem_h}"
        ) from e

    if not Path(domtblout).exists():
        return result

    with open(domtblout, encoding="utf-8", newline="") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            target_name = parts[0]
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            domain_evalue = float(parts[12])
            domain_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            result[target_name].append({
                "target": target_name,
                "score": full_score,
                "evalue": full_evalue,
                "domain_score": domain_score,
                "domain_evalue": domain_evalue,
                "length": ali_to - ali_from + 1,
                "hmm_start": hmm_from, "hmm_end": hmm_to,
                "ali_start": ali_from, "ali_end": ali_to,
                "env_start": env_from, "env_end": env_to,
            })

    Path(domtblout).unlink(missing_ok=True)
    return result


def _parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse a FASTA file into {header_first_token: sequence} (shared impl)."""
    return parse_fasta(fasta_path)


def _scan_one_candidate(
    genome_id: str,
    marker_id: str,
    fasta_path: str,
    hmm_path: Optional[str],
    cpus: int,
    tmp_dir: str,
    min_hmm_score: float,
) -> Tuple[str, str, Dict[str, List[Dict]], int]:
    """Scan a single (genome, marker) pair.

    Pure, side-effect-free worker: it runs ``_run_hmmsearch`` and returns the
    raw hits plus the qualifying-hit count. The caller is responsible for
    merging results into the occupancy matrix and the hmm cache, so this stays
    safe to run concurrently (threads or processes) without shared-state races.
    The hmmsearch domtblout filename embeds the genome/marker stems plus the
    module-level ``_RUN_UUID``, so concurrent scans never collide on temp files.
    """
    if hmm_path is None:
        return genome_id, marker_id, {}, 0
    hits = _run_hmmsearch(fasta_path, hmm_path, cpus, tmp_dir)
    n_hits = sum(
        1 for t, th in hits.items()
        if any(h["score"] >= min_hmm_score for h in th)
    )
    return genome_id, marker_id, hits, n_hits


class AdaptiveMarkerSelector:
    def __init__(self, config: SelectionConfig, marker_ids: Optional[List[str]] = None):
        self.config = config
        self.marker_ids: List[str] = marker_ids if marker_ids is not None else []
        self.min_occupancy_threshold = config.min_occupancy
        self.max_markers = config.max_markers
        self.min_hmm_score = config.min_hmm_score
        self._hmm_cache: Optional[Dict[str, Dict[str, Dict[str, List[Dict]]]]] = None

    def scan_all_candidates(self, genomes: List[Genome]) -> OccupancyMatrix:
        matrix = OccupancyMatrix(genomes=[g.id for g in genomes], cogs=list(self.marker_ids))

        # Build the full (genome, marker) task list, then run the hmmsearch
        # Scans concurrently. hmmsearch is an external subprocess that spends
        # Most of its wall-clock time blocked on I/O, so running the scans in a
        # Thread pool overlaps them with no change to the per-pair result. Each
        # Task writes a uniquely-named domtblout
        # (``hmmsearch_<genome>_<marker>_<_RUN_UUID>.domtblout``) and returns
        # Its own hits, so there is no shared mutable state to race on; the
        # Caller merges results serially afterwards.
        tasks = [
            (g.id, marker_id, g.effective_protein_path,
             hmm_path_for(self.config.marker_hmm_dir, marker_id))
            for g in genomes
            for marker_id in self.marker_ids
        ]
        logger.info(
            f"    Scanning {len(genomes)} genomes x {len(self.marker_ids)} "
            f"HMM profiles ({len(tasks)} pairs)..."
        )

        def _worker(task):
            gid, mid, fasta, hmm = task
            return _scan_one_candidate(
                gid, mid, fasta, hmm,
                self.config.cpus, self.config.tmp_dir, self.min_hmm_score,
            )

        results = parallel_map(_worker, tasks, max_workers=max(1, min(self.config.cpus, len(tasks))))

        for genome_id, marker_id, hits, n_hits in results:
            if self._hmm_cache is not None:
                self._hmm_cache.setdefault(genome_id, {})[marker_id] = hits
            if n_hits == 0:
                matrix.set(genome_id, marker_id, GeneState.ABSENT)
            elif n_hits == 1:
                matrix.set(genome_id, marker_id, GeneState.SINGLE_COPY)
            else:
                matrix.set(genome_id, marker_id, GeneState.MULTI_COPY)
        return matrix

    def scan_target(self, fasta_path: str, marker_id: str) -> List[Dict]:
        hmm_file = hmm_path_for(self.config.marker_hmm_dir, marker_id)
        if hmm_file is None:
            return []
        hits = _run_hmmsearch(fasta_path, hmm_file, 1, self.config.tmp_dir)
        return [h for th in hits.values() for h in th if h["score"] >= self.min_hmm_score]

    def calculate_occupancy_scores(self, matrix: OccupancyMatrix) -> Dict[str, float]:
        scores: Dict[str, float] = {}
        for cog in matrix.cogs:
            states = matrix.get_column(cog)
            single = sum(1 for s in states if s == GeneState.SINGLE_COPY)
            if self.config.quality_weighted and matrix.genome_qualities:
                tw = sum(matrix.get_quality(g) for g in matrix.genomes)
                wc = sum(matrix.get_quality(g) for g in matrix.genomes if matrix.get(g, cog) == GeneState.SINGLE_COPY)
                scores[cog] = wc / tw if tw > 0 else 0
            else:
                scores[cog] = single / len(matrix.genomes) if matrix.genomes else 0
        return scores

    def select_optimal_marker_set(self, matrix: OccupancyMatrix, occupancy_scores: Dict[str, float]) -> SelectedMarkerSet:
        strat = self.config.strategy
        if strat == SelectionStrategy.GREEDY:
            return self._greedy_selection(matrix, occupancy_scores)
        elif strat == SelectionStrategy.INFO_MAX:
            return self._information_maximization(matrix, occupancy_scores)
        elif strat == SelectionStrategy.RATE_BALANCED:
            return self._rate_balanced_selection(matrix, occupancy_scores)
        elif strat == SelectionStrategy.SPARSE_OPTIMIZED:
            return self._sparse_optimized_selection(matrix, occupancy_scores)
        raise ValueError(f"Unknown strategy: {strat}")

    def _greedy_selection(self, matrix, scores):
        eligible = {c: s for c, s in scores.items() if s >= self.min_occupancy_threshold}
        sorted_cogs = sorted(eligible.items(), key=lambda x: x[1], reverse=True)
        selected = [c for c, _ in sorted_cogs[:self.max_markers]]
        final = {c: eligible[c] for c in selected}
        if not selected:
            raise NoMarkerAvailableError(
                f"No markers meet min_occupancy {self.min_occupancy_threshold}. Max found: "
                f"{max(scores.values()) if scores else 0:.4f}. Check {self.config.marker_hmm_dir}."
            )
        return SelectedMarkerSet(markers=selected, occupancy_scores=final, quality_scores={},
                                 strategy=SelectionStrategy.GREEDY, total_candidates=len(self.marker_ids),
                                 mean_occupancy=sum(final.values())/len(final) if final else 0)

    def _information_maximization(self, matrix, scores):
        eligible = {c: s for c, s in scores.items() if s >= self.min_occupancy_threshold}
        if not eligible:
            raise NoMarkerAvailableError("No markers meet min_occupancy threshold")
        # Precompute, once, the set of single-copy genomes per marker. The
        # Redundancy penalty scored below queries pairwise overlaps of these
        # Sets; computing them up-front turns every overlap query from O(G)
        # (a full genome scan) into an O(1) set operation, removing the
        # O(K^2 * G) blow-up of the inner selection loop. The numeric result is
        # Identical to scanning the matrix each time.
        sc_sets = self._single_copy_sets(matrix)
        selected, seen, rw = [], set(), 0.5
        for _ in range(min(self.max_markers, len(eligible))):
            best_cog, best_gain = None, -1.0
            for cog, sc in eligible.items():
                if cog in seen:
                    continue
                if not selected:
                    gain = sc
                else:
                    # ``selected`` holds marker ids (str), which is exactly what
                    # ``_compute_occupancy_overlap`` expects for its second
                    # Argument. Passing the occupancy *score* (a float) here was a
                    # Bug: it never matched a column in the matrix, so the
                    # Redundancy penalty was always 0 and InfoMax collapsed into
                    # Plain greedy.
                    red = sum(self._compute_occupancy_overlap(matrix, cog, c2, sc_sets) for c2 in selected) / len(selected)
                    gain = sc - rw * red
                if gain > best_gain:
                    best_gain, best_cog = gain, cog
            if best_cog is None or best_gain <= 0:
                break
            selected.append(best_cog); seen.add(best_cog)
        final = {c: eligible[c] for c in selected}
        return SelectedMarkerSet(markers=selected, occupancy_scores=final, quality_scores={},
                                 strategy=SelectionStrategy.INFO_MAX, total_candidates=len(scores),
                                 mean_occupancy=sum(final.values())/len(final) if final else 0)

    def _single_copy_sets(self, matrix) -> Dict[str, set]:
        """Precompute ``{marker_id: {genome_id that is SINGLE_COPY}}`` once.

        Used by ``_information_maximization`` so that the pairwise redundancy
        penalty can be evaluated with O(1) set operations instead of re-scanning
        every genome for every pair.
        """
        out: Dict[str, set] = {}
        for c in matrix.cogs:
            s = set()
            for g in matrix.genomes:
                if matrix.get(g, c) == GeneState.SINGLE_COPY:
                    s.add(g)
            out[c] = s
        return out

    def _compute_occupancy_overlap(self, matrix, c1, c2, _sc_sets=None):
        if _sc_sets is not None:
            s1 = _sc_sets.get(c1, set())
            s2 = _sc_sets.get(c2, set())
        else:
            s1, s2 = set(), set()
            for g in matrix.genomes:
                if matrix.get(g, c1) == GeneState.SINGLE_COPY: s1.add(g)
                if matrix.get(g, c2) == GeneState.SINGLE_COPY: s2.add(g)
        union = len(s1 | s2)
        if union == 0:
            return 0.0
        return len(s1 & s2) / union

    def _rate_balanced_selection(self, matrix, scores):
        eligible = {c: s for c, s in scores.items() if s >= self.min_occupancy_threshold}
        if not eligible: raise NoMarkerAvailableError("No markers meet min_occupancy threshold")
        composite = {c: 0.6*eligible[c] + 0.4*self._estimate_rate_balance(matrix, c) for c in eligible}
        sel = [c for c,_ in sorted(composite.items(), key=lambda x:x[1], reverse=True)[:self.max_markers]]
        final = {c: eligible[c] for c in sel}
        return SelectedMarkerSet(markers=sel, occupancy_scores=final, quality_scores={},
                                 strategy=SelectionStrategy.RATE_BALANCED, total_candidates=len(scores),
                                 mean_occupancy=sum(final.values())/len(final) if final else 0)

    def _estimate_rate_balance(self, matrix, cog):
        states = matrix.get_column(cog)
        t = len(states) if states else 1
        sr = sum(1 for s in states if s == GeneState.SINGLE_COPY)/t
        mr = sum(1 for s in states if s == GeneState.MULTI_COPY)/t
        dev = abs(mr - 0.05)/(0.05 + 0.01)
        return sr * max(0.0, 1.0 - dev)

    def _sparse_optimized_selection(self, matrix, scores):
        eligible = {c: s for c, s in scores.items() if s > 0}
        if not eligible: raise NoMarkerAvailableError("No markers with any hits found")
        sel = [c for c,_ in sorted(eligible.items(), key=lambda x:x[1], reverse=True)[:self.max_markers]]
        final = {c: eligible[c] for c in sel}
        return SelectedMarkerSet(markers=sel, occupancy_scores=final, quality_scores={},
                                 strategy=SelectionStrategy.SPARSE_OPTIMIZED, total_candidates=len(matrix.cogs),
                                 mean_occupancy=sum(final.values())/len(final) if final else 0)


class AdaptiveMarkerSelectionModule:
    def __init__(self, config: SelectionConfig):
        self.config = config
        self._hmm_cache: Optional[Dict[str, Dict[str, Dict[str, List[Dict]]]]] = None
        self._marker_ids: Optional[List[str]] = None
        md = config.marker_hmm_dir
        if config.marker_mode == "hmm" and md and Path(md).is_dir():
            try:
                self._marker_ids = gather_marker_ids(md)
                logger.info(f"  [{config.marker_mode}] discovered {len(self._marker_ids)} HMM profiles in {md}")
            except FileNotFoundError as e:
                logger.warning(f"  [{config.marker_mode}] cant read marker_hmm_dir: {e}")
                self._marker_ids = []
            self.selector = AdaptiveMarkerSelector(config, marker_ids=self._marker_ids)
        else:
            self.selector = AdaptiveMarkerSelector(config, marker_ids=[])

    def run(self, genomes: List[Genome], phase_context: object = None) -> MarkerSelectionResult:
        logger.info(f"{canonical_log_tag('1')} Marker Selection (mode={self.config.marker_mode}; {len(self.selector.marker_ids)} profiles)")
        if phase_context is not None and hasattr(phase_context, "adaptive_params") and phase_context.adaptive_params:
            p = phase_context.adaptive_params
            self.selector.min_occupancy_threshold = p.min_occupancy
            self.selector.max_markers = p.max_markers
            self.selector.min_hmm_score = p.min_hmm_score
            logger.info(f"  Applied adaptive params: occ>={p.min_occupancy}, max={p.max_markers}, hmm_score>={p.min_hmm_score}")

        self._hmm_cache = {}
        # Share the cache dict with the selector: scan_all_candidates stores each
        # Genome's hmmsearch hits into ``self._hmm_cache`` (seen via the selector),
        # Which ``_compute_quality_scores`` later reads to derive real hmm_score and
        # Alignment-length coherence. Without this aliasing the selector's own cache
        # Attr stays ``None`` and the per-bitscore signal is lost.
        self.selector._hmm_cache = self._hmm_cache
        matrix = self.selector.scan_all_candidates(genomes)
        occ = self.selector.calculate_occupancy_scores(matrix)
        ms = self.selector.select_optimal_marker_set(matrix, occ)
        # Evaluate quality from the real per-bitscore + domtblout alignment-length
        # Coherence (sequences have not been extracted yet at this point; the HMM
        # Hit lengths are the authoritative length-coherence signal).
        qs = self._compute_quality_scores(ms, genomes, matrix, hmm_cache=self._hmm_cache)
        logger.info(f"  Selected {len(ms.markers)} markers (mean occ {ms.mean_occupancy:.4f})")
        return MarkerSelectionResult(marker_set=ms, quality_scores=qs, occupancy_matrix=matrix, resolved_sequences={})

    def run_gtdb_tk(self, genomes, markers_dir, tmp_dir="/tmp", adaptive_params=None,
                    gene_trees_dir: Optional[str] = None):
        from markerfinder.utils.gtdb_tk_markers import load_gtdb_markers, build_gene_trees
        logger.info(f"{canonical_log_tag('1')} GTDB-TK marker loading (ar53/bac120)")
        gids = [g.id for g in genomes]
        matrix, mseq, occ = load_gtdb_markers(markers_dir, genome_ids=gids)

        # Apply MAG-aware adaptive parameters when provided by Phase 0.
        min_occupancy = self.config.min_occupancy
        max_markers = self.config.max_markers
        min_hmm_score = self.config.min_hmm_score
        if adaptive_params is not None:
            min_occupancy = adaptive_params.min_occupancy
            max_markers = adaptive_params.max_markers
            min_hmm_score = adaptive_params.min_hmm_score
            logger.info(
                f"  Applied adaptive params: occ>={min_occupancy}, "
                f"max={max_markers}, hmm_score>={min_hmm_score}"
            )

        elig = {m: s for m, s in occ.items() if s >= min_occupancy}
        if not elig:
            raise NoMarkerAvailableError(
                f"No GTDB markers meet min_occupancy {min_occupancy}. "
                f"Max occ: {max(occ.values()) if occ else 0:.4f}"
            )
        # Respect max_markers even in gtdb_tk mode for consistency with hmm mode.
        if len(elig) > max_markers:
            elig = dict(sorted(elig.items(), key=lambda x: x[1], reverse=True)[:max_markers])
            logger.info(f"  Trimmed to top {max_markers} markers by occupancy")

        ms = SelectedMarkerSet(
            markers=list(elig.keys()),
            occupancy_scores=elig,
            quality_scores={},
            strategy=SelectionStrategy.GREEDY,
            total_candidates=len(occ),
            mean_occupancy=sum(elig.values()) / len(elig) if elig else 0,
        )
        qs = self._compute_quality_scores(ms, genomes, matrix)
        logger.info(f"  Loaded {len(ms.markers)} GTDB markers (mean occ {ms.mean_occupancy:.4f})")
        # 基因树缓存目录优先使用调用方给的 <output>/Phase4_trees/gene_trees,
        # 使 filter/infer 能复用同一批树; 未提供时退回 tmp 子目录(旧行为).
        gtd = gene_trees_dir or os.path.join(tmp_dir, "gene_trees")
        gts = build_gene_trees(mseq, gtd, method="fasttree", cpus=self.config.cpus)
        return (
            MarkerSelectionResult(
                marker_set=ms, quality_scores=qs, occupancy_matrix=matrix, resolved_sequences={}
            ),
            mseq,
            gts,
        )

    def extract_marker_sequences(self, genomes, marker_set, matrix, tmp_dir="/tmp"):
        mg = {m: [] for m in marker_set.markers}
        hd = self.config.marker_hmm_dir
        if not hd or not Path(hd).is_dir():
            logger.warning("  extract_marker_sequences: marker_hmm_dir not configured")
            return mg
        for gn in genomes:
            gseqs = _parse_fasta(gn.effective_protein_path)
            for mid in marker_set.markers:
                if matrix.get(gn.id, mid) == GeneState.ABSENT:
                    continue
                hf = hmm_path_for(hd, mid)
                if hf is None:
                    continue
                hits = (self._hmm_cache.get(gn.id, {}).get(mid) if self._hmm_cache is not None else None)
                if hits is None:
                    hits = _run_hmmsearch(gn.effective_protein_path, hf, 1, tmp_dir)
                allh = [h for th in hits.values() for h in th if h["score"] >= self.selector.min_hmm_score]
                if not allh:
                    continue
                # One representative sequence per (genome, marker). GTDB's
                # Conserved-protein set is single-copy by design, so when a
                # Genome has several qualifying hits (typically the same marker
                # On chromosome + plasmid, overlapping domains, etc.) we keep
                # Only the highest-bitscore hit and drop the rest. Writing all
                # Of them under the same genome-level id would collide later in
                # The supermatrix concatenation (SeqIO.to_dict raises
                # "Duplicate key"). The occupancy matrix still counts all
                # Qualifying hits upstream (see scan_all_candidates), so the
                # MULTI_COPY / SINGLE_COPY assessment used by adaptive
                # Thresholds is unaffected.
                #
                # Deduplicate by target first: hmmsearch domtblout may report
                # Multiple domain rows for the same protein with identical or
                # Near-identical scores. We keep the best score per target, then
                # Pick the globally best target if several distinct proteins
                # Qualify.
                best_per_target: Dict[str, Dict] = {}
                for h in allh:
                    tgt = h["target"]
                    if tgt not in best_per_target or h["score"] > best_per_target[tgt]["score"]:
                        best_per_target[tgt] = h
                allh = list(best_per_target.values())

                # 同源多拷贝(paralogy)单例选择策略: 同一基因组对某个标记可能命中
                # 多个副本(paralogs)——domtblout 中表现为多个 target 行(已按 target
                # 去重保留每个 target 的最佳分值, 见上方 best_per_target)。这里再从这些
                # 不同 target 中取 *全局最高分值* 的单个副本作为该基因组的代表序列
                # ("单例")。这是刻意的: 我们不为同一基因组保留多个副本(那会膨胀超级
                # 矩阵并歪曲占用矩阵), 也不随机丢弃——永远保留最可信(hit 分值最高)的
                # 那个, 保证结果可复现且偏向高质量比对。不同基因组中同名但序列不同的
                # 合法情况不受影响(各基因组独立选其自身最佳副本)。
                best = max(allh, key=lambda h: h["score"])
                if len(allh) > 1:
                    dropped = [(h["target"], h["score"]) for h in allh if h is not best]
                    dropped_str = ", ".join(f"{t} (score={s:.1f})" for t, s in dropped)
                    logger.debug(
                        f"  [hmm-dedup] {mid} / {gn.id}: {len(allh)} distinct targets passing "
                        f"score>={self.selector.min_hmm_score}; keeping best "
                        f"target={best['target']} score={best['score']:.1f}, "
                        f"dropping {len(allh)-1} lower-scoring: {dropped_str}"
                    )
                self._append_marker_seq(mg, mid, gn.id, best, gseqs)

        # Defensive stance: every marker's extracted sequence list must have
        # Unique genome ids. If this trips, extraction produced the exact
        # Collision the dedup above is meant to prevent — fail loudly with the
        # Offending marker/genomes rather than letting it silently reach Phase 4.
        for mid, seqs in mg.items():
            ids = [s["id"] for s in seqs]
            if len(ids) != len(set(ids)):
                from collections import Counter
                dup = sorted(k for k, c in Counter(ids).items() if c > 1)
                raise RuntimeError(
                    f"Marker extraction produced duplicate genome ids for {mid}: {dup}"
                )
        return mg

    @staticmethod
    def _append_marker_seq(mg, mid, gid, hit, gseqs):
        tgt = hit["target"]
        a0, a1 = hit.get("ali_start", 0), hit.get("ali_end", 0)
        seq = ""
        if tgt in gseqs:
            full = gseqs[tgt]
            seq = full[a0-1:a1] if a0 > 0 and a1 > 0 and a1 <= len(full) else full
        # NOTE: keep ``id`` equal to the bare genome id (no mid/tgt suffix).
        # Downstream supermatrix concatenation in phylogenetic_inference.py groups
        # Rows by id; a per-marker suffix would make each (genome, marker)
        # Pair look like a separate taxon (N_markers * N_genomes taxa instead of
        # N_genomes). ``target``/``score``/``evalue`` keep the per-hit detail.
        entry = {"id": gid, "genome_id": gid, "target": tgt,
                 "seq": seq, "score": hit["score"], "evalue": hit["evalue"]}
        mg[mid].append(entry)

    def _compute_quality_scores(self, marker_set, genomes, matrix,
                                hmm_cache: Optional[Dict[str, Dict[str, Dict[str, List[Dict]]]]] = None) -> Dict[str, MarkerQualityScore]:
        """Per-marker quality score.

        Uses only quantities that are actually computed by this phase:
        - hmm_score: mean over genomes of the best per-bitscore (full_score from
          hmmsearch). Stored as the *raw* mean bitscore so the quality formula's
          ``hmm_norm = min(hmm_score / 100.0, 1.0)`` normalizes on the intended
          0-100+ scale (bitscores of ~100 saturate the component). Real, per-marker.
        - occupancy: from the occupancy score. Real, per-marker.
        - length_ratio: alignment-length coherence across genomes (1 - CV of the
          per-genome HMM hit lengths, in [0,1]). Real, per-marker.
        - copy_number_cv: from the occupancy matrix (MULTI_COPY fraction). Real.

        Two components remain placeholders because they require downstream phases
        not yet executed at this point:
        - phylogenetic_informativeness: needs a gene-tree / site-composition
          analysis (Phase 3, ``run_infer``). Centered at 0.5 until then.
        - hgt_risk_score: needs the HGT filter (Phase 2). This is back-filled by
          ``MarkerFinderPipeline.run_filter`` from the HGT evaluation's
          ``overall_risk`` once available (see ``HGTFilterModule`` output).
        """
        hmm_cache = hmm_cache or self._hmm_cache or {}
        genome_ids = [g.id for g in genomes]

        out: Dict[str, MarkerQualityScore] = {}
        for cog in marker_set.markers:
            states = matrix.get_column(cog)
            s = sum(1 for st in states if st == GeneState.SINGLE_COPY)
            m = sum(1 for st in states if st == GeneState.MULTI_COPY)
            t = len(states) if states else 1

            # --- hmm_score: mean of the best per-genome bitscore, normalized ---
            best_scores: List[float] = []
            ali_lengths: List[float] = []
            for gid in genome_ids:
                hits_by_target = (hmm_cache.get(gid) or {}).get(cog)
                all_hits = [h for th in (hits_by_target or {}).values() for h in th]
                if all_hits:
                    best = max(all_hits, key=lambda h: h.get("score", 0.0))
                    best_scores.append(float(best.get("score", 0.0)))  # Non-negative bitscore
                    ali_lengths.append(float(best.get("length", 0)))
            if best_scores:
                # Store the RAW mean bitscore. The quality formula normalizes via
                # ``hmm_norm = min(hmm_score / 100.0, 1.0)``, so a score around/below
                # 100 gives a proportional contribution and stronger hits saturate it.
                hmm_score = sum(best_scores) / len(best_scores)
            else:
                hmm_score = 0.0

            # --- length_ratio: alignment-length coherence across genomes (1 - CV) ---
            if len(ali_lengths) >= 2:
                mean_len = sum(ali_lengths) / len(ali_lengths)
                var_len = sum((x - mean_len) ** 2 for x in ali_lengths) / len(ali_lengths)
                sd_len = var_len ** 0.5
                cv = sd_len / mean_len if mean_len > 0 else 0.0
                length_ratio = max(0.0, min(1.0, 1.0 - cv))
            elif len(ali_lengths) == 1:
                length_ratio = 1.0
            else:
                length_ratio = 0.0

            out[cog] = MarkerQualityScore(
                marker_id=cog,
                hmm_score=hmm_score,
                occupancy=marker_set.occupancy_scores.get(cog, 0.0),
                length_ratio=length_ratio,
                copy_number_cv=m / t if t else 0,
                phylogenetic_informativeness=None,  # Unmeasured until a gene tree exists (Phase 3 back-fill in run_infer)
                hgt_risk_score=0.0,                 # Placeholder: back-filled from HGT filter (Phase 2 back-fill in run_filter)
                functional_category="",
            )
        return out
