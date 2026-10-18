"""GTDB-TK marker loading and per-marker gene-tree building.

MarkerFinder consumes the conserved-protein sequences that GTDB-TK produces for
a genome set instead of discovering markers via HMM search or DIAMOND
self-search:

* ``ar53`` (53 archaeal) / ``bac120`` (120 bacterial) marker proteins, one
  FASTA file per marker, containing every genome's (raw, UNALIGNED) sequence.
* Optionally the GTDB-TK concatenated species tree (Newick), used as the
  reference species tree for the HGT phylogenetic-step screening.

No HMM/DIAMOND marker discovery is performed — the marker set is exactly the
GTDB conserved-protein set. IMPORTANT: GTDB-TK's ar53/bac120 output FASTA files
are NOT aligned or trimmed, so before building a gene tree each marker must be
(1) aligned with MAFFT and (2) gap-trimmed with trimal. These two steps are
required preparation for a valid tree, not "extra" homology search.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from markerfinder.models.genome import GeneState, OccupancyMatrix
from markerfinder.utils.io import parse_fasta
from markerfinder.utils.concurrency import parallel_map
import json
import logging

logger = logging.getLogger(__name__)


def _parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse a FASTA file into {header_first_token: sequence} (shared impl)."""
    return parse_fasta(fasta_path)


def load_gtdb_markers(
    markers_dir: str,
    genome_ids: Optional[List[str]] = None,
) -> Tuple[OccupancyMatrix, Dict[str, List[dict]], Dict[str, float]]:
    """Load per-marker raw/unaligned FASTA files produced from GTDB-TK output.

    The directory should contain one raw/unaligned FASTA per conserved protein
    (filename stem = marker id), each holding every genome's sequence.

    Returns:
        (occupancy_matrix, marker_sequences, occupancy_scores)
        - marker_sequences: {marker_id: [{"id", "genome_id", "seq"},...]}
        - occupancy_scores: {marker_id: fraction of genomes present}
    """
    markers_dir = Path(markers_dir)
    if not markers_dir.is_dir():
        raise FileNotFoundError(f"GTDB markers directory not found: {markers_dir}")

    fasta_files = sorted(markers_dir.glob("*.faa")) + sorted(markers_dir.glob("*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(f"No per-marker FASTA files (*.faa/*.fasta) in {markers_dir}")

    raw: Dict[str, Dict[str, str]] = {}
    all_genomes: List[str] = list(genome_ids) if genome_ids else []
    for f in fasta_files:
        marker_id = f.stem
        seqs = _parse_fasta(str(f))
        raw[marker_id] = seqs
        for gid in seqs:
            if gid not in all_genomes:
                all_genomes.append(gid)

    matrix = OccupancyMatrix(genomes=all_genomes, cogs=list(raw.keys()))
    marker_sequences: Dict[str, List[dict]] = {}
    occupancy_scores: Dict[str, float] = {}

    for marker_id, seqs in raw.items():
        present = [gid for gid in all_genomes if gid in seqs]
        for gid in all_genomes:
            if gid in seqs:
                matrix.set(gid, marker_id, GeneState.SINGLE_COPY)
            else:
                matrix.set(gid, marker_id, GeneState.ABSENT)
        # 关键正确性修复: FASTA 表头/tip 标签使用**裸 genome id**,
        # 而不是 "{gid}_{marker_id}".
        # 下游 _write_fasta_for_seqs 以 seq["id"] 作为 FASTA 表头 ->
        # _concatenate_fasta_alignments 按表头对齐各行形成超级矩阵; 基因树 tip
        # 也来自同一字段(HGT 筛选)或 genome_id(build_gene_trees). 若用
        # "{gid}_{marker_id}" 作为表头, 每个 (基因组, 标记) 会成为独立的分类
        # 单元, 超级矩阵按列拼接错位, ASTRAL 合并基因树时同一基因组的不同 tip
        # 也无法对应. 统一为裸 gid 后, 同一基因组在不同 marker 的对齐行 / 基因
        # 树 tip 才能被正确按列拼接 / 对应. marker_id 另存为独立字段供需要时使用.
        marker_sequences[marker_id] = [
            {"id": gid, "genome_id": gid, "marker_id": marker_id, "seq": seqs[gid]}
            for gid in present
        ]
        occupancy_scores[marker_id] = (len(present) / len(all_genomes)) if all_genomes else 0.0

    logger.info(
        f"  Loaded {len(raw)} GTDB markers across {len(all_genomes)} genomes "
        f"(mean occupancy: {sum(occupancy_scores.values()) / len(occupancy_scores):.3f})"
    )
    return matrix, marker_sequences, occupancy_scores


def split_user_msa(
    msa_fasta: str,
    boundaries: List[Tuple[str, int]],
    out_dir: str,
) -> None:
    """Split a GTDB-TK concatenated per-genome MSA into per-marker FASTA files.

    GTDB-TK's ``*.user_msa.fasta`` concatenates every marker's alignment into a
    single sequence per genome (markers back-to-back, in a fixed order). Given
    the marker boundaries ``[(name, aligned_length),...]`` in that same order,
    this writes one aligned FASTA per marker into ``out_dir``.

    ``boundaries`` can be generated from GTDB-TK's packaged marker set, e.g.::

        from gtdb_toolkit import MARKERS # pseudo; see scripts/split_gtdbtk_msa.py
        bounds = [(m.name, m.aligned_len) for m in MARKERS if m.domain == "bac120"]

    Args:
        msa_fasta: path to GTDB-TK ``*.user_msa.fasta``
        boundaries: ordered list of (marker_name, aligned_length)
        out_dir: directory to write per-marker ``<name>.faa`` files
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_marker: Dict[str, List[str]] = {name: [] for name, _ in boundaries}
    try:
        with open(msa_fasta, encoding="utf-8", newline="") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    header = line[1:].split()[0]
                    for name, _ in boundaries:
                        per_marker[name].append(f">{header}")
                else:
                    pos = 0
                    for name, length in boundaries:
                        chunk = line[pos:pos + length]
                        per_marker[name].append(chunk)
                        pos += length
    except FileNotFoundError:
        logger.error(f"  GTDB-TK user MSA not found: {msa_fasta}")
        return

    for name, lines in per_marker.items():
        marker_target = (out_dir / f"{name}.faa").resolve()
        marker_target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    logger.info(f"  Split {msa_fasta} into {len(boundaries)} per-marker FASTA files in {out_dir}")


def _gene_tree_input_id(seqs: List[dict]) -> str:
    """Sha256 over the (genome id, sequence) pairs a gene tree is built from.

    Order-independent: the same marker file re-emitted in another order is the
    same input, so the checkpoint is not invalidated by a rewrite.
    """
    digest = hashlib.sha256()
    for s in sorted(seqs, key=lambda x: str(x.get("genome_id", ""))):
        digest.update(
            f"{s.get('genome_id', '')}\0{s.get('seq', '')}\n".encode("utf-8")
        )
    return digest.hexdigest()


def _build_one_gene_tree(
    marker_id: str,
    seqs: List[dict],
    out_dir: str,
    method: str = "fasttree",
    min_tips: int = 4,
    trim_threshold: float = 0.2,
    cpus: int = 1,
) -> Tuple[str, Optional[str]]:
    """Build the gene tree for a single marker (worker for ``build_gene_trees``).

    Contains the exact per-marker logic previously inlined in ``build_gene_trees``
    (write unaligned FASTA, MAFFT align, trimal trim, FastTree/IQ-TREE build,
    checkpoint reuse). Returns ``(marker_id, nwk_path_or_None)`` so the caller
    can merge results with no shared mutable state. Each marker uses its own
    ``{marker_id}.faa/.aln/.aln.trim/.nwk`` files, so concurrent workers never
    touch the same file.
    """
    out_dir = Path(out_dir)
    if len(seqs) < min_tips:
        logger.info(f"  Skipping gene tree for {marker_id} ({len(seqs)} genomes < {min_tips})")
        return marker_id, None

    unaligned = out_dir / f"{marker_id}.faa"
    aln = out_dir / f"{marker_id}.aln.faa"
    trimmed = out_dir / f"{marker_id}.aln.trim.faa"
    nwk_path = out_dir / f"{marker_id}.nwk"
    stamp_path = out_dir / f"{marker_id}.nwk.input_sha256"
    input_id = _gene_tree_input_id(seqs)

    # Checkpoint: skip rebuild only if the existing tree was built from THIS
    # Input. A checkpoint keyed by marker id alone is handed back to the next run
    # That reuses the directory (``--force`` over an existing output directory is
    # Documented), so a tree inferred from different sequences would silently
    # Become this run's gene tree — and the HGT screen would grade against it.
    if nwk_path.exists() and nwk_path.stat().st_size > 0:
        stamped = (
            stamp_path.read_text(encoding="utf-8").strip()
            if stamp_path.exists() else ""
        )
        if stamped == input_id:
            logger.info(f"  Reusing existing gene tree for {marker_id}")
            return marker_id, str(nwk_path)
        logger.info(
            f"  Existing gene tree for {marker_id} was built from different "
            f"sequences; rebuilding"
        )

    unaligned_lines: List[str] = []
    for s in seqs:
        # Use the genome id (not the marker-specific id) as the tip
        # Label so gene-tree tips match the taxonomy table and stay
        # Consistent across markers for downstream concatenation.
        unaligned_lines.append(f">{s['genome_id']}\n{s['seq']}\n")
    unaligned_target = unaligned.resolve()
    unaligned_target.write_text("".join(unaligned_lines), encoding="utf-8", newline="\n")

    try:
        # 1) MAFFT alignment (writes alignment to stdout)
        # MAFFT's multi-threaded mode (--thread N, N>1) produces non-deterministic
        # Alignments: thread scheduling affects gap placement in edge cases, which
        # Propagates through trimal → PIS → marker_quality_score. Pin to a single
        # Thread here so repeated runs are byte-identical (reproducibility, test_10).
        # Per-marker parallelism is already provided by parallel_map at the caller.
        aln_target = aln.resolve()
        with aln_target.open("w", encoding="utf-8") as af:
            subprocess.run(
                ["mafft", "--thread", "1", "--quiet", "--auto", str(unaligned)],
                stdout=af, check=True, timeout=900,
            )
        if not aln.exists() or aln.stat().st_size == 0:
            raise subprocess.CalledProcessError(1, "mafft")

        # 2) trimal trimming
        subprocess.run(
            ["trimal", "-in", str(aln), "-out", str(trimmed), "-gt", str(trim_threshold)],
            check=True, timeout=300,
        )
        if not trimmed.exists() or trimmed.stat().st_size == 0:
            logger.warning(
                f"  trimal produced an empty alignment for {marker_id}; "
                f"falling back to untrimmed MAFFT alignment"
            )
            trimmed = aln

        # Persist the per-marker informativeness as a sidecar
        # Next to the tree so the pipeline can back-fill marker_summary
        # Without re-reading alignments.
        try:
            from markerfinder.utils.informative_sites import (
                effective_columns,
                gap_fraction,
                parsimony_informative_sites,
            )

            aln_map: Dict[str, str] = {}
            _sid = ""
            with open(trimmed, encoding="utf-8") as _fh:
                for _line in _fh:
                    _line = _line.strip()
                    if _line.startswith(">"):
                        _sid = _line[1:].split()[0]
                        aln_map[_sid] = ""
                    elif _sid:
                        aln_map[_sid] += _line
            if aln_map:
                sidecar = (out_dir / f"{marker_id}.pis").resolve()
                sidecar.write_text(
                    json.dumps({
                        "pis": parsimony_informative_sites(aln_map),
                        "effective_columns": effective_columns(aln_map),
                        "gap_fraction": round(gap_fraction(aln_map), 6),
                    }),
                    encoding="utf-8", newline="\n",
                )
        except Exception as e:  # Noqa: BLE001 - PIS is auxiliary, never fatal
            logger.debug(f"  PIS sidecar for {marker_id} not written: {e}")

        # Composition / GC diagnostics on the SAME trimmed
        # Alignment, appended to the sidecar. Kept in its own guard so a
        # Composition failure cannot destroy the PIS evidence beside it.
        try:
            from markerfinder.modules.composition import (
                marker_composition_metrics,
            )

            if aln_map:
                _metrics = marker_composition_metrics(aln_map, protein=True)
                _sidecar = (out_dir / f"{marker_id}.pis").resolve()
                if _sidecar.exists():
                    _payload = json.loads(
                        _sidecar.read_text(encoding="utf-8")
                    )
                else:
                    _payload = {}
                _payload.update({
                    "rcv": _metrics["rcv"],
                    "gc_bias": _metrics["gc_bias"],
                    "n_sequences": _metrics["n_sequences"],
                })
                _sidecar.write_text(
                    json.dumps(_payload), encoding="utf-8", newline="\n",
                )
        except Exception as e:  # Noqa: BLE001 - auxiliary, never fatal
            logger.debug(f"  Composition metrics for {marker_id} not written: {e}")

        # 3) tree building on the trimmed alignment
        if method == "iqtree":
            subprocess.run(
                [
                    "iqtree3", "-quiet", "-nt", str(cpus), "-s", str(trimmed),
                    "-pre", str(out_dir / marker_id),
                    "-seed", "1",
                ],
                capture_output=True, text=True, check=True, timeout=600,
            )
            final = out_dir / f"{marker_id}.treefile"
            if final.exists():
                final.rename(nwk_path)
        else:
            # 关键正确性修复: GTDB 比对是**蛋白序列**, 必须用 WAG 模型
            # (-wag). 原代码使用的 "-nt N" 在 FastTree 中是 *核苷酸 JC69 模型*
            # 开关(并非线程数), 用在蛋白比对上会选错模型(或因子进程名解析错误
            # 而建树失败), 导致所有 GTDB 基因树拓扑/支长不可信.
            # - 移除 "-nt N" 位置参数与核苷酸模型开关;
            # - 使用 "-wag" 蛋白模型(与本项目 hgt_filter.py /
            # Phylogenetic_inference.py 的其它 FastTree 调用一致);
            # - FastTree 本身不支持多线程, 故不传递线程参数;
            # - 依次尝试 fasttree / FastTree 两种大小写, 首个命中即使用
            # (与项目其它处 FastTree 调用风格统一).
            ft_candidates = ["fasttree", "FastTree"]
            ft_ok = False
            # Why the last attempt failed. The summary warning below has to be
            # able to name a cause, and it used to interpolate `ft_last_err` --
            # a name that was never assigned on any path. So the "FastTree
            # unavailable" branch, which is the normal branch when fasttree is
            # not installed, raised NameError instead of warning: a degradation
            # that was supposed to be quiet and recoverable turned into a crash
            # inside the logger call that was supposed to describe it.
            ft_last_err = "no FastTree candidate ran"
            for ft_exe in ft_candidates:
                try:
                    if ft_exe == "fasttree":
                        subprocess.run(
                            ["fasttree", "-wag", "-quiet", "-seed", "1", "-out", str(nwk_path), str(trimmed)],
                            capture_output=True, text=True, check=True, timeout=600,
                        )
                    else:
                        subprocess.run(
                            ["FastTree", "-wag", "-quiet", "-seed", "1", "-out", str(nwk_path), str(trimmed)],
                            capture_output=True, text=True, check=True, timeout=600,
                        )
                    ft_ok = True
                    break
                except FileNotFoundError as e:
                    ft_last_err = f"{ft_exe} is not on PATH ({e})"
                    continue
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                    ft_last_err = f"{ft_exe} failed ({e})"
                    logger.warning(
                        f"  FastTree({ft_exe}) gene-tree build failed for "
                        f"{marker_id}: {e}"
                    )
                    break
            if not ft_ok:
                logger.warning(
                    f"  FastTree unavailable/failed for {marker_id}; "
                    f"gene tree will be skipped. Last error: {ft_last_err}"
                )
        path = str(nwk_path) if nwk_path.exists() and nwk_path.stat().st_size > 0 else None
        if path:
            stamp_path.write_text(input_id + "\n", encoding="utf-8")
        return marker_id, path
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.warning(f"  Gene tree build failed for {marker_id}: {e}")
        return marker_id, None
    finally:
        for tmp in (unaligned, aln):
            tmp.unlink(missing_ok=True)


def build_gene_trees(
    marker_sequences: Dict[str, List[dict]],
    out_dir: str,
    method: str = "fasttree",
    min_tips: int = 4,
    trim_threshold: float = 0.2,
    cpus: int = 1,
) -> Dict[str, Optional[str]]:
    """Build a gene tree for every marker from GTDB-TK's (unaligned) sequences.

    GTDB-TK's ar53/bac120 output FASTA files are raw, UNALIGNED protein
    sequences, so before tree building we must:
      1. align with MAFFT (``--auto``)
      2. trim gaps with trimal (``-gt <trim_threshold>``)
    Only then do we run the tree builder (FastTree or IQ-TREE). This is
    required preparation for a valid gene tree, not "extra" homology search —
    no HMM/DIAMOND marker discovery is performed (markers come straight from
    GTDB-TK).

    Args:
        marker_sequences: {marker_id: [{"id", "genome_id", "seq"},...]}
        out_dir: directory for intermediate alignments + ``<marker_id>.nwk``
        method: "fasttree" (default) or "iqtree"
        min_tips: skip markers with fewer than this many genomes
        trim_threshold: trimal gap threshold (fraction of gaps per column;
            columns with more gaps than this are removed)

    Returns:
        {marker_id: newick_path or None}
    """
    # MAFFT + trimal are required to turn GTDB-TK's raw sequences into a
    # Trimmable alignment. Fail loudly if either is missing.
    for bin_name in ("mafft", "trimal"):
        if not shutil.which(bin_name):
            raise FileNotFoundError(
                f"gtdb_tk mode requires '{bin_name}' on PATH to align/trim "
                f"GTDB-TK ar53/bac120 sequences before tree building. "
                f"Install it (e.g. `conda install -c bioconda {bin_name}`)."
            )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Each marker's MAFFT+trimal+FastTree/IQ-TREE chain is an independent,
    # I/O-bound external-tool job writing its own {marker_id}.faa/.aln/.nwk
    # Files, so the markers can be built concurrently. Results are identical to
    # The previous serial loop; the per-marker checkpoint reuse (skip when a
    # Non-empty {marker_id}.nwk already exists) is preserved inside the worker.
    tasks = [(marker_id, seqs) for marker_id, seqs in marker_sequences.items()]

    def _worker(task):
        marker_id, seqs = task
        return _build_one_gene_tree(
            marker_id, seqs, str(out_dir), method, min_tips, trim_threshold, cpus,
        )

    built = parallel_map(_worker, tasks, max_workers=max(1, min(cpus, len(tasks))))

    trees: Dict[str, Optional[str]] = {marker_id: path for marker_id, path in built}

    n_built = sum(1 for v in trees.values() if v)
    logger.info(f"  Built {n_built}/{len(marker_sequences)} per-marker gene trees (MAFFT+trimal prep)")
    return trees
