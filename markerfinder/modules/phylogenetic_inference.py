from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from markerfinder.phases import canonical_log_tag
from markerfinder.config import PhylogeneticConfig
from markerfinder.utils.io import write_fasta
from markerfinder.models.alignment import PartitionEntry, PartitionFile
from markerfinder.models.pipeline_types import (
    CoalescentResult,
    ConflictReport,
    PhylogeneticResult,
    PhylogeneticInferenceError,
    SupermatrixResult,
)
from markerfinder.models.tree import Tree, TreeRecommendation
from markerfinder.utils.concurrency import parallel_map
import logging

logger = logging.getLogger(__name__)


# 标记比对长度的回退默认值 (单一事实源).
# 当比对文件不可读 / 无有效比对时, 占位分区 Nexus 与串联长度回退使用的标记长度.
# 与旧逻辑 "AlignIO 读取失败时回退 300" 的取值保持一致.
DEFAULT_MARKER_LENGTH = 300


def _write_fasta_for_seqs(seqs: List[dict], output_path: str) -> None:
    """将序列字典列表写入FASTA文件 (复用 utils.io.write_fasta)。"""
    sequences = {
        str(s.get("id", "unknown")): str(s.get("seq", s.get("target", "")))
        for s in seqs
    }
    write_fasta(sequences, output_path)


def _run_mafft(input_fasta: str, output_fasta: str, cpus: int = 1) -> bool:
    """运行 MAFFT 多序列比对。成功返回 True；工具缺失/失败则返回 False（上层据此
    跳过该标记或得到 0 个基因树），**不**抛异常——MAFFT 缺失意味着该标记无法比对，
    属“无数据”而非“错误结论”，逐标记跳过是预期行为（Coalescent 返回 0 基因树、
    Supermatrix 跳过该标记）。失败时会输出明确 warning 日志，避免静默跳过。
    （注：HGT/占用矩阵/质量等“伪装成正常结果”的静默降级已在各自模块修复，
    此处遵循“可优雅降级为空”的设计契约，与测试 ``test_run_with_seqs`` 一致。）"""
    try:
        out_target = Path(output_fasta).resolve()
        with out_target.open("w", encoding="utf-8", newline="\n") as out:
            # MAFFT multi-threaded mode is non-deterministic (thread scheduling
            # Affects gap placement). Pin to 1 thread for reproducibility;
            # Per-marker parallelism is already provided by parallel_map.
            subprocess.run(
                ["mafft", "--auto", "--quiet", "--thread", "1", input_fasta],
                stdout=out, stderr=subprocess.DEVNULL, check=True, timeout=600,
            )
        return Path(output_fasta).exists()
    except FileNotFoundError as e:
        logger.warning(
            "MAFFT not found on PATH; skipping alignment for this marker. "
            "Install it (e.g. `conda install -c bioconda mafft`) to enable "
            "phylogenetic inference."
        )
        return False
    except subprocess.CalledProcessError as e:
        logger.warning(f"MAFFT failed for {input_fasta}: {e}")
        return False
    except subprocess.TimeoutExpired as e:
        logger.warning(f"MAFFT timed out for {input_fasta}: {e}")
        return False


def _run_trimal(input_aln: str, output_aln: str) -> bool:
    """运行 trimAl 比对修剪。"""
    try:
        subprocess.run(
            ["trimal", "-automated1", "-in", input_aln, "-out", output_aln],
            capture_output=True, check=True, timeout=300,
        )
        return Path(output_aln).exists()
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.warning(f"  trimAl failed: {e}")
        return False


def _run_iqtree3(
    aln_path: str,
    prefix: str,
    cpus: int = 1,
    ufboot: int = 1000,
    fast_mode: bool = False,
    partition_path: Optional[str] = None,
    timeout: int = 3600,
) -> Optional[Tree]:
    """运行 IQ-TREE3 ML建树。"""
    cmd: List[str] = [
        "iqtree3",
        "-s", aln_path,
        "-m", "LG+F" if fast_mode else "MFP",
        "-nt", str(cpus),
        "-pre", prefix,
        "--quiet",
        "-seed", "1",
    ]
    if fast_mode:
        # Ultrafast bootstrap (-B/-bnni) 与 -fast 不兼容:
        # "Ultrafast bootstrap does not work with -fast... option"
        # -fast 是类似 FastTree 的快速近似搜索,本身不提供 bootstrap,故省略 -B。
        cmd.insert(1, "-fast")
    else:
        # Iqtree3 要求 -B 的重复数 >= 1000；低于此值会直接拒绝。
        ufboot = max(int(ufboot), 1000)
        cmd.extend(["-B", str(ufboot), "-bnni"])
    if partition_path and Path(partition_path).exists():
        cmd.extend(["-p", partition_path])

    try:
        subprocess.run(
            ["iqtree3"] + cmd[1:],
            capture_output=True, text=True, check=True, timeout=timeout,
        )
    except FileNotFoundError:
        logger.error("  iqtree3 not found. Please install IQ-TREE3.")
        return None
    except subprocess.CalledProcessError as e:
        logger.warning(f"  IQ-TREE3 failed: {e.stderr[:300]}")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("  IQ-TREE3 timed out")
        return None

    tree_path = f"{prefix}.treefile"
    if Path(tree_path).exists():
        return Tree.read(tree_path)
    return None


def _run_fasttree(input_aln: str, output_tree: str) -> Optional[Tree]:
    """运行 FastTree2 快速建树 (蛋白序列使用 WAG 模型).

    conda/各平台的 FastTree 二进制名大小写不统一(fasttree / FastTree),
    故依次尝试两种名称,首个命中即使用.
    """
    for exe in ("fasttree", "FastTree"):
        tree = _run_fasttree_named(exe, input_aln, output_tree)
        if tree is not None:
            return tree
    logger.warning("  FastTree2 not found (tried fasttree / FastTree). Please install FastTree2.")
    return None


def _run_fasttree_named(exe: str, input_aln: str, output_tree: str) -> Optional[Tree]:
    """Run one concrete FastTree binary name; None means try the next name.

    ``-seed 1`` pins FastTree's internal RNG so that repeated runs on the
    same alignment produce identical trees (required for reproducibility).
    """
    try:
        if exe == "fasttree":
            subprocess.run(
                ["fasttree", "-wag", "-quiet", "-seed", "1", "-out", output_tree, input_aln],
                capture_output=True, check=True, timeout=600,
            )
        else:
            subprocess.run(
                ["FastTree", "-wag", "-quiet", "-seed", "1", "-out", output_tree, input_aln],
                capture_output=True, check=True, timeout=600,
            )
        if Path(output_tree).exists():
            return Tree.read(output_tree)
        return None
    except FileNotFoundError:
        return None
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.warning(f"  FastTree(%s) failed: {e}", exe)
        return None


def _run_astral(gene_trees_file: str, output_tree: str, timeout: int = 1800) -> Optional[Tree]:
    """运行 ASTRAL-III 合并物种树。"""
    try:
        subprocess.run(
            ["astral", "-i", gene_trees_file, "-o", output_tree, "-t", "1"],
            capture_output=True, text=True, check=True, timeout=timeout,
        )
        if Path(output_tree).exists():
            # ASTRAL 以 -t 1 运行: LPP 后验概率作为常规 Newick 节点标签
            # (分支支持度) 输出, 直接保留. strip_nhx_annotations 作为防御性
            # 后路保留, 以防上游被改为 -t 2 产出含 NHX 注释的树.
            try:
                raw = Path(output_tree).read_text(encoding="utf-8").strip()
                from markerfinder.utils.tree_utils import strip_nhx_annotations
                clean = strip_nhx_annotations(raw)
                if clean:
                    return Tree(newick=clean)
            except Exception:
                pass
            return Tree.read(output_tree)
        return None
    except FileNotFoundError:
        logger.error("  astral not found. Please install ASTRAL-III.")
        return None
    except subprocess.CalledProcessError as e:
        logger.warning(f"  ASTRAL failed: {e.stderr[:300]}")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("  ASTRAL timed out")
        return None


def _concatenate_fasta_alignments(
    alignments: Dict[str, str],
    output_path: str,
    placeholder: str = "-",
) -> Tuple[List[str], int]:
    """串联多个FASTA比对。

    Returns:
        (species_order, total_sites)
    """
    from Bio import SeqIO

    all_species: List[str] = []
    aln_seqs: Dict[str, Dict[str, str]] = {}
    aln_lengths: Dict[str, int] = {}

    for marker_id, aln_path in alignments.items():
        if not Path(aln_path).exists():
            continue
        records = SeqIO.to_dict(SeqIO.parse(aln_path, "fasta"))
        aln_seqs[marker_id] = {k: str(v.seq) for k, v in records.items()}
        aln_lengths[marker_id] = len(next(iter(records.values())).seq) if records else 0
        for sid in records:
            if sid not in all_species:
                all_species.append(sid)

    all_species.sort()
    marker_order = list(alignments.keys())

    supermatrix_lines: List[str] = []
    for species in all_species:
        parts: List[str] = []
        for mid in marker_order:
            seq = aln_seqs.get(mid, {}).get(species, placeholder * aln_lengths.get(mid, 0))
            parts.append(seq)
        full_seq = "".join(parts)
        supermatrix_lines.append(f">{species}")
        for i in range(0, len(full_seq), 80):
            supermatrix_lines.append(full_seq[i:i+80])

    total_sites = sum(aln_lengths.values())
    matrix_target = Path(output_path).resolve()
    matrix_target.write_text("\n".join(supermatrix_lines) + "\n", encoding="utf-8", newline="\n")
    return all_species, total_sites


def _write_partition_nexus(
    marker_lengths: Dict[str, int],
    output_path: str,
) -> None:
    """写入 Nexus 分区文件。"""
    lines: List[str] = ["#nexus", "begin sets;"]
    pos = 1
    for mid, length in marker_lengths.items():
        end = pos + length - 1
        lines.append(f"  charset {mid} = {pos}-{end};")
        pos = end + 1
    lines.append("end;")
    partition_target = Path(output_path).resolve()
    partition_target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


class SupermatrixInference:
    """串联法系统发育推断。"""

    def __init__(self, config: PhylogeneticConfig):
        self.config = config

    def run(
        self, marker_genes: Dict[str, List[dict]], genomes: List[object]
    ) -> SupermatrixResult:
        logger.info("Mode A: Supermatrix Phylogenetic Inference")
        tmp_dir = self.config.tmp_dir
        # Create the scratch directory here rather than relying on somebody
        # Higher up having made it: this module is also driven directly (per-step
        # Runs, library use), and a missing directory otherwise surfaces as a
        # FileNotFoundError on a path the caller never typed.
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        prefix = self.config.output_prefix

        alignments: Dict[str, str] = {}
        marker_lengths: Dict[str, int] = {}

        for marker_id, seqs in marker_genes.items():
            if not seqs or len(seqs) < 4:
                continue

            fasta_path = os.path.join(tmp_dir, f"{marker_id}.faa")
            aln_path = os.path.join(tmp_dir, f"{marker_id}.aln")
            trim_path = os.path.join(tmp_dir, f"{marker_id}.aln.trim")

            _write_fasta_for_seqs(seqs, fasta_path)

            if not _run_mafft(fasta_path, aln_path, self.config.cpus):
                continue

            trimmed = _run_trimal(aln_path, trim_path)
            final_aln = trim_path if trimmed else aln_path

            # 比对长度必须与 _concatenate_fasta_alignments 实际串联所用的长度严格
            # 一致: 二者都从同一个 final_aln 文件读取, 否则分区 Nexus 坐标会与
            # 串联 FASTA 错配(旧逻辑在 AlignIO 读取失败时回退 300, 但串联按真实
            # 长度拼接, 导致 IQ-TREE 分区错位). 读取失败意味着比对文件损坏 ->
            # 该标记不应以伪造长度进入分区, 直接排除.
            try:
                from Bio import SeqIO
                _records = list(SeqIO.parse(final_aln, "fasta"))
                if not _records:
                    raise ValueError("empty alignment")
                marker_lengths[marker_id] = len(_records[0].seq)
            except Exception as e:
                logger.warning(
                    f"  Marker {marker_id} alignment unreadable ({e}); "
                    f"excluding it from the supermatrix."
                )
                continue

            alignments[marker_id] = final_aln

        if not alignments:
            # 无任何可用比对时绝不能返回占位树 ";": dataclass 实例恒为
            # 真值, 会被下游可用性判定当成有效串联树写盘并判定整趟运行"成功".
            # Tree=None 让 report 跳过文件写出、pipeline 判定 species_tree_source
            # =none, CLI 以 EXIT_DATA_ERROR 明确告知无可用物种树.
            logger.warning(
                "  No markers could be aligned. No supermatrix tree will be produced."
            )
            partition = PartitionFile(entries=[
                PartitionEntry(name=mid, start=1, end=DEFAULT_MARKER_LENGTH, model="AUTO")
                for mid in marker_genes
            ])
            return SupermatrixResult(
                alignment=None, partition=partition,
                tree=None, marker_order=list(marker_genes.keys()),
            )

        concat_path = os.path.join(tmp_dir, f"{prefix}.concat.fasta")
        species_order, total_sites = _concatenate_fasta_alignments(
            alignments, concat_path
        )

        partition_entries = []
        pos = 1
        for mid in alignments:
            length = marker_lengths.get(mid, DEFAULT_MARKER_LENGTH)
            partition_entries.append(PartitionEntry(name=mid, start=pos, end=pos + length - 1))
            pos += length
        partition = PartitionFile(entries=partition_entries)

        partition_path = os.path.join(tmp_dir, f"{prefix}.partition.nex")
        partition.write_nexus(partition_path)

        tree_prefix = os.path.join(tmp_dir, f"{prefix}_concat")
        tree = _run_iqtree3(
            concat_path, tree_prefix,
            cpus=self.config.cpus,
            ufboot=self.config.ufboot_replicates,
            fast_mode=self.config.fast_mode,
            partition_path=partition_path,
            timeout=self.config.iqtree_timeout,
        )

        if tree is None:
            logger.warning("  IQ-TREE3 failed. Trying FastTree2...")
            fasttree_path = os.path.join(tmp_dir, f"{prefix}.fasttree")
            tree = _run_fasttree(concat_path, fasttree_path)

        if tree is None:
            # 外部建树工具(IQ-TREE3 / FastTree)缺失或失败, 绝不能静默返回空树
            #; 假装物种树构建成功 —— 那是把工具失败伪装成生物学结论. 显式硬
            # 失败, 让调用方(Pipeline)以可定位的异常终止.
            raise PhylogeneticInferenceError(
                "Species-tree inference failed: both IQ-TREE3 and FastTree "
                "produced no usable tree (tools missing or crashed)."
            )

        return SupermatrixResult(
            alignment=None,
            partition=partition,
            tree=tree,
            marker_order=list(alignments.keys()),
            best_model="MFP",
            avg_ufboot=None,
        )


def _coalescent_build_one(
    marker_id: str,
    seqs: List[dict],
    tmp_dir: str,
    use_fasttree: bool,
    gene_trees_dir: str,
    cpus: int,
    ufboot_replicates: int,
    fast_mode: bool,
    iqtree_timeout: int,
) -> Tuple[str, Optional[Tree]]:
    """Build the gene tree for a single marker (worker for ``CoalescentInference``).

    Faithful port of the per-marker build branch of ``CoalescentInference.run``.
    Returns ``(marker_id, Tree_or_None)``. The caller handles cache reuse (read
    of an already-built ``{marker_id}.nwk``) before dispatching builds, so this
    worker only ever builds (and writes its own per-marker ``.nwk``/``.fasttree``
    file — distinct filenames across markers, hence safe to run concurrently).
    """
    fasta_path = os.path.join(tmp_dir, f"{marker_id}.faa")
    aln_path = os.path.join(tmp_dir, f"{marker_id}.aln")
    trim_path = os.path.join(tmp_dir, f"{marker_id}.aln.trim")

    _write_fasta_for_seqs(seqs, fasta_path)

    if not _run_mafft(fasta_path, aln_path, cpus):
        return marker_id, None

    trimmed = _run_trimal(aln_path, trim_path)
    final_aln = trim_path if trimmed else aln_path

    if use_fasttree:
        # 写入缓存目录供后续复用
        if gene_trees_dir:
            Path(gene_trees_dir).mkdir(parents=True, exist_ok=True)
            tree_path = str(Path(gene_trees_dir) / f"{marker_id}.nwk")
        else:
            tree_path = os.path.join(tmp_dir, f"{marker_id}.fasttree")
        tree = _run_fasttree(final_aln, tree_path)
    else:
        tree_prefix = os.path.join(tmp_dir, f"{marker_id}")
        tree = _run_iqtree3(
            final_aln, tree_prefix,
            cpus=1,
            ufboot=ufboot_replicates,
            fast_mode=fast_mode,
            timeout=iqtree_timeout,
        )
    return marker_id, tree


class CoalescentInference:
    """合并法系统发育推断。"""

    def __init__(self, config: PhylogeneticConfig, gene_trees_dir: str = ""):
        self.config = config
        self.gene_trees_dir = gene_trees_dir

    def run(
        self, marker_genes: Dict[str, List[dict]], genomes: List[object]
    ) -> CoalescentResult:
        logger.info("Mode B: Coalescent Phylogenetic Inference")
        tmp_dir = self.config.tmp_dir
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        prefix = self.config.output_prefix
        # Explicit two-option gene-tree builder selection. `gene_tree_builder`
        # ("fasttree" / "iqtree") takes precedence over the legacy `use_fasttree`
        # Boolean for coalescent inference.
        gene_tree_builder = getattr(self.config, "gene_tree_builder", None)
        if gene_tree_builder is None:
            gene_tree_builder = "fasttree" if self.config.use_fasttree else "iqtree"
        use_fasttree = gene_tree_builder == "fasttree"
        logger.info(f"  Gene-tree builder for coalescent inference: {gene_tree_builder}")

        gene_trees: Dict[str, Tree] = {}
        # Per-marker provenance — did the coalescent leg reuse
        # The exact FastTree the HGT/concat leg produced? Shared upstream
        # Caps recommendation confidence (R2 / Error-A era lesson).
        gene_tree_provenance: Dict[str, str] = {}

        # Cache reuse (read of an already-built {marker_id}.nwk from a prior
        # Phase, e.g. HGT) is done serially here — it is cheap and must happen
        # Before any rebuild so we never clobber a tree another phase produced
        # In this shared directory.
        build_tasks: List[Tuple[str, List[dict]]] = []
        for marker_id, seqs in marker_genes.items():
            if not seqs or len(seqs) < 4:
                continue

            # 优先复用系统发育步骤缓存的 FastTree 基因树（仅 FastTree 路径可复用）
            if use_fasttree and self.gene_trees_dir:
                cache_path = Path(self.gene_trees_dir) / f"{marker_id}.nwk"
                if cache_path.exists() and cache_path.stat().st_size > 0:
                    try:
                        gene_trees[marker_id] = Tree.read(str(cache_path))
                        gene_tree_provenance[marker_id] = "cached"
                        continue
                    except Exception:
                        pass  # 缓存损坏，重建

            gene_tree_provenance[marker_id] = "rebuilt"
            build_tasks.append((marker_id, seqs))

        # Build the remaining (uncached) markers concurrently. Each build writes
        # Its own per-marker.nwk/.fasttree file, so concurrent workers never
        # Touch the same file. Results are identical to the old serial loop.
        def _worker(task):
            mid, seqs = task
            return _coalescent_build_one(
                mid, seqs, tmp_dir, use_fasttree, self.gene_trees_dir,
                self.config.cpus, self.config.ufboot_replicates,
                self.config.fast_mode, self.config.iqtree_timeout,
            )

        built = parallel_map(_worker, build_tasks, max_workers=max(1, min(self.config.cpus, len(build_tasks))))
        for mid, tree in built:
            if tree is not None:
                gene_trees[mid] = tree

        if not gene_trees:
            # 与串联法同理: 占位树 ";" 会被下游真值判定当成有效物种树,
            # 改为 None + source="none", 让 CLI 明确报"无可用物种树".
            logger.warning("  No gene trees could be built. No coalescent tree will be produced.")
            return CoalescentResult(
                species_tree=None,
                gene_trees={}, filtered_gene_trees={},
                n_total_genes=0, n_passed_filter=0,
                gene_tree_provenance=gene_tree_provenance,
            )

        logger.info(f"  Built {len(gene_trees)} gene trees")

        filtered_trees = self._filter_gene_trees(gene_trees)
        logger.info(f"  {len(filtered_trees)} gene trees passed quality filter")

        species_tree = self._run_astral(filtered_trees, prefix, tmp_dir)
        species_tree_source = "astral" if species_tree is not None else None

        if species_tree is None:
            # ASTRAL-III 失败时不再退化为单棵基因树; 明确记录失败并返回空结果,
            # 由报告阶段决定是否输出物种树文件.
            logger.error(
                "  ASTRAL-III failed to produce a coalescent species tree. "
                "No coalescent species tree will be output. "
                "Check the ASTRAL-III installation/log for details."
            )
            species_tree = None
            species_tree_source = "none"

        return CoalescentResult(
            species_tree=species_tree,
            gene_trees=gene_trees,
            filtered_gene_trees=filtered_trees,
            gene_tree_provenance=gene_tree_provenance,
            n_total_genes=len(gene_trees),
            n_passed_filter=len(filtered_trees),
            avg_quartet_support=None,
            species_tree_source=species_tree_source or "none",
        )

    def _filter_gene_trees(self, gene_trees: Dict[str, Tree]) -> Dict[str, Tree]:
        filtered: Dict[str, Tree] = {}
        # Optional, gentle support gate. ``None`` keeps the historical behaviour
        # (only the ``n_tips >= 4`` filter). When set, gene trees whose average
        # Local support falls below the threshold are dropped so they do not
        # Degrade the coalescent (ASTRAL-III) species tree.
        min_support = getattr(self.config, "min_gene_tree_support", None)
        unmeasurable: List[str] = []
        for mid, tree in gene_trees.items():
            if tree.n_tips < 4:
                continue
            if min_support is not None:
                avg_support = tree.get_average_support()
                if avg_support is None:
                    # "support could not be read" is neither "support 0"
                    # Nor a number at all. Comparing it raised
                    # TypeError: '<' not supported between instances of
                    # 'NoneType' and 'float' and killed this filter (canonical
                    # Phase 3). A quality gate
                    # Fails closed, but it must say which trees it could not
                    # Judge instead of dropping them silently.
                    unmeasurable.append(mid)
                    continue
                if avg_support < min_support:
                    logger.debug(
                        f"  Dropping gene tree {mid}: avg support "
                        f"{avg_support:.2f} < {min_support:.2f}"
                    )
                    continue
            filtered[mid] = tree
        if unmeasurable:
            shown = ", ".join(sorted(unmeasurable)[:10])
            more = f" (+{len(unmeasurable) - 10} more)" if len(unmeasurable) > 10 else ""
            logger.warning(
                f"  Gene-tree support gate: {len(unmeasurable)} gene tree(s) "
                f"dropped because average local support was NOT MEASURABLE "
                f"(unparseable Newick, or ete3 unavailable in this "
                f"interpreter): {shown}{more}"
            )
        return filtered

    def _run_astral(self, gene_trees: Dict[str, Tree], prefix: str, tmp_dir: str) -> Optional[Tree]:
        if not gene_trees:
            return None

        gene_trees_file = os.path.join(tmp_dir, f"{prefix}.genetrees")
        tree_lines: List[str] = [tree.to_newick() + "\n" for tree in gene_trees.values()]
        trees_target = Path(gene_trees_file).resolve()
        trees_target.write_text("".join(tree_lines), encoding="utf-8", newline="\n")

        output_tree = os.path.join(tmp_dir, f"{prefix}.astral.tree")
        return _run_astral(gene_trees_file, output_tree, self.config.astral_timeout)


class ConflictDetector:
    """Conflict detection between concatenation and coalescent species trees.

    The current implementation reports dataset-level topology disagreement
    (normalized RF distance and quartet agreement) and a human-readable
    recommendation. Branch-level conflict classification (HGT vs. ILS vs.
    method bias) is **not implemented**; the corresponding fields in
    ``ConflictReport`` are kept as reserved placeholders for backward
    compatibility.
    """

    # The band thresholds used to be constructor parameters that
    # Nothing could reach -- the module always built ConflictDetector with the
    # Defaults, and quartet_threshold had zero readers anywhere in the package
    # (a knob that turns nothing is worse than no knob, because it is documented
    # As if it worked). The moderate bound now comes from the same config field
    # The recommendation uses, so the two cannot drift apart, and
    # Quartet_threshold is gone.
    def __init__(self, config: object = None, rf_threshold: Optional[float] = None):
        self.config = config
        self.rf_threshold = (
            float(rf_threshold)
            if rf_threshold is not None
            else float(getattr(config, "recommend_rf_moderate", 0.3))
        )

    @property
    def low_bound(self) -> float:
        """Lower edge of the "moderate" band: one config field, one source."""
        return float(getattr(self.config, "recommend_rf_low", 0.1))

    def detect_conflicts(
        self,
        concat_tree: Tree,
        astral_tree: Tree,
        gene_trees: Dict[str, Tree],
    ) -> ConflictReport:
        logger.info(f"{canonical_log_tag('3')} Conflict Detection & Assessment")
        from markerfinder.utils.tree_utils import calculate_rf_distance, get_quartet_topology

        # Concat and ASTRAL trees are mutual references —
        # Both must be fully resolved before either is used to judge the other.
        from markerfinder.utils.reference import validate_reference_tree

        v1 = validate_reference_tree(concat_tree.newick, astral_tree.newick)
        v2 = validate_reference_tree(astral_tree.newick, concat_tree.newick)
        if not (v1.ok and v2.ok):
            bad = v1 if not v1.ok else v2
            logger.warning(
                f"  Conflict detection skipped: reference illegal "
                f"[{bad.code}] {bad.detail}"
            )
            return ConflictReport(
                rf_distance=None,
                normalized_rf=None,
                quartet_agreement=None,
                topology_note=(
                    f"Conflict detection not performed: reference illegal "
                    f"[{bad.code}] {bad.detail}"
                ),
            )

        raw_rf, norm_rf = calculate_rf_distance(concat_tree.newick, astral_tree.newick)
        if norm_rf is not None:
            logger.info(f"  Normalized RF distance: {norm_rf:.3f}")
        else:
            logger.info("  Normalized RF distance: N/A (computation failed)")

        quartet_agreement = self._quartet_agreement(
            concat_tree.newick, astral_tree.newick
        )
        if quartet_agreement is not None:
            logger.info(f"  Quartet agreement: {quartet_agreement:.3f}")

        topology_note = self._topology_note(
            concat_tree.newick, astral_tree.newick, norm_rf, quartet_agreement
        )
        if topology_note:
            logger.info(f"  {topology_note}")

        # Branch-level conflict classification fields (n_conflicting_branches,
        # Conflicting_branches, gene_tree_agreement, conflict_summary) are
        # Intentionally left at their default placeholder values because the
        # Current implementation only computes dataset-level RF/quartet metrics.
        return ConflictReport(
            rf_distance=raw_rf,
            normalized_rf=norm_rf,
            quartet_agreement=quartet_agreement,
            topology_note=topology_note,
        )

    def _quartet_agreement(
        self,
        newick1: str,
        newick2: str,
        max_quartets: int = 200,
    ) -> Optional[float]:
        """计算两棵树在所有 quartet 上的拓扑一致比例。

        使用无根 quartet 子树的 RF 距离判断一致性, 避免因 rooting/分支顺序
        导致的字符串差异. 当 RF=0 时该值应为 1.0.
        """
        try:
            from itertools import combinations
            # Single sanctioned entry point.
            from markerfinder.utils.etree import require_ete3

            EteTree = require_ete3().Tree
            have_ete3 = True
        except Exception:
            # Ete3 缺失时改用纯 Python split-set quartet 拓扑
            # 比较（utils.etree.quartet_topology，规范化串相等即一致）。
            have_ete3 = False

        if not have_ete3:
            return self._quartet_agreement_pure_python(
                newick1, newick2, max_quartets=max_quartets
            )

        try:
            from itertools import combinations

            t1 = EteTree(newick1)
            t2 = EteTree(newick2)
            tips1 = set(n.name for n in t1.get_leaves() if n.name)
            tips2 = set(n.name for n in t2.get_leaves() if n.name)
            common = sorted(tips1 & tips2)
            if len(common) < 4:
                return None

            all_quartets = list(combinations(common, 4))
            n_total = len(all_quartets)
            if n_total > max_quartets:
                # 系统性采样, 避免大数据集下组合爆炸.
                step = n_total / max_quartets
                selected = [all_quartets[int(i * step)] for i in range(max_quartets)]
            else:
                selected = all_quartets

            consistent = 0
            total = 0
            for quartet in selected:
                try:
                    q1 = t1.copy()
                    q1.prune(list(quartet), preserve_branch_length=False)
                    q2 = t2.copy()
                    q2.prune(list(quartet), preserve_branch_length=False)
                    rf_result = q1.robinson_foulds(q2, unrooted_trees=True)
                    rf = rf_result[0] if rf_result else -1
                    total += 1
                    if rf == 0:
                        consistent += 1
                except Exception:
                    continue

            if total == 0:
                return None
            return consistent / total
        except Exception as e:
            logger.debug(f"  Quartet agreement calculation failed: {e}")
            return None

    @staticmethod
    def _quartet_agreement_pure_python(
        newick1: str,
        newick2: str,
        max_quartets: int = 200,
    ) -> Optional[float]:
        """Ete3 缺失时的 quartet agreement 回退。

        与 ete3 路径同语义：两棵树在 quartet 上的诱导拓扑（规范化 2-2 划分串）
        相等即一致；未解析（星状）诱导拓扑两两相等也计为一致，与 ete3
        prune 后字符串相等的行为一致。>64 共享尖端显式不可测。
        """
        from itertools import combinations

        from markerfinder.utils import etree as _etree

        tips1 = _etree.tip_set(newick1)
        tips2 = _etree.tip_set(newick2)
        common = sorted(tips1 & tips2)
        if len(common) < 4:
            return None
        if len(common) > _etree.MAX_PURE_PYTHON_TIPS:
            logger.debug(
                "  Quartet agreement not measurable: "
                f"{len(common)} shared tips exceed pure-Python cap "
                f"({_etree.MAX_PURE_PYTHON_TIPS})"
            )
            return None

        all_quartets = list(combinations(common, 4))
        n_total = len(all_quartets)
        if n_total > max_quartets:
            # 系统性采样, 与 ete3 路径相同的确定性采样规则.
            step = n_total / max_quartets
            selected = [all_quartets[int(i * step)] for i in range(max_quartets)]
        else:
            selected = all_quartets

        consistent = 0
        total = 0
        for quartet in selected:
            q1 = _etree.quartet_topology(newick1, quartet)
            q2 = _etree.quartet_topology(newick2, quartet)
            if q1 is None or q2 is None:
                continue
            total += 1
            if q1 == q2:
                consistent += 1

        if total == 0:
            return None
        return consistent / total

    def _topology_note(
        self,
        newick1: str,
        newick2: str,
        norm_rf: float,
        quartet_agreement: Optional[float],
    ) -> str:
        """根据 RF 距离与 quartet agreement 生成人类可读的拓扑说明。"""
        if norm_rf < 0:
            return "RF distance could not be computed; conflict assessment unavailable."
        if norm_rf == 0.0:
            note = (
                "Concatenation and coalescent trees are topologically identical "
                "(RF distance = 0)."
            )
            if newick1.strip() != newick2.strip():
                note += (
                    " Differences in the Newick strings are due to rooting or branch "
                    "ordering, not topology."
                )
            if quartet_agreement is not None and quartet_agreement < 1.0:
                note += (
                    f" Quartet agreement is {quartet_agreement:.3f}, which may reflect "
                    "very short internal branches or numeric noise."
                )
            return note
        # Both bounds come from configuration, and the low bound is
        # The SAME field recommend_tree reads -- no threshold/3 restatement, which
        # Would let the prose and the verdict disagree.
        if norm_rf < self.low_bound:
            return (
                "Low topological conflict between concatenation and coalescent trees."
            )
        if norm_rf < self.rf_threshold:
            return (
                "Moderate topological conflict; coalescent tree is preferred for "
                "species tree inference."
            )
        return (
            "Strong topological conflict; potential HGT, ILS, or method bias signal."
        )


class PhylogeneticInferenceModule:
    """路径三: 双策略系统发育推断主模块。"""

    def __init__(self, config: PhylogeneticConfig, gene_trees_dir: str = ""):
        self.config = config
        self.supermatrix = SupermatrixInference(config)
        self.coalescent = CoalescentInference(config, gene_trees_dir=gene_trees_dir)
        # Hand the detector the config so its bands are the run's bands.
        self.conflict_detector = ConflictDetector(config=config)

    def run(
        self,
        marker_genes: Dict[str, List[dict]],
        genomes: List[object],
        phase_context: object = None,
    ) -> PhylogeneticResult:
        result = PhylogeneticResult()

        if not marker_genes:
            logger.warning("  No marker genes provided. Skipping phylogenetic inference.")
            return result

        result.supermatrix = self.supermatrix.run(marker_genes, genomes)

        if self.config.coalescent_mode == "off":
            logger.info("  Coalescent mode: OFF - skipping coalescent inference")
            result.coalescent = None
        elif self.config.coalescent_mode in ("post-filter", "always"):
            logger.info(f"  Coalescent mode: {self.config.coalescent_mode.upper()}")
            result.coalescent = self.coalescent.run(marker_genes, genomes)

        if result.supermatrix and result.coalescent:
            if result.supermatrix.tree and result.coalescent.species_tree:
                concat_tree = result.supermatrix.tree
                astral_tree = result.coalescent.species_tree
                if concat_tree.n_tips >= 4 and astral_tree.n_tips >= 4:
                    result.conflict_report = self.conflict_detector.detect_conflicts(
                        concat_tree,
                        astral_tree,
                        result.coalescent.gene_trees,
                    )
                    # Measure how independent the two legs are.
                    from markerfinder.models.pipeline_types import IndependenceReport

                    provenance = getattr(result.coalescent, "gene_tree_provenance", {}) or {}
                    cached = sum(1 for v in provenance.values() if v == "cached")
                    shared_ratio = (cached / len(provenance)) if provenance else 0.0
                    result.conflict_report.independence = IndependenceReport(
                        marker_set_jaccard=1.0 if result.coalescent.gene_trees else 0.0,
                        shared_gene_tree_ratio=shared_ratio,
                        same_trimming_regime=False,
                        ref_built_from_tested_markers=True,
                    )
                else:
                    logger.info("  Skipping conflict detection: one or both trees have <4 tips")

        # The dataset-level recommendation is computed INSIDE
        # The run and stored on the result, not left as an off-run helper. The
        # Pipeline summary asserts properties of this verdict ("confidence
        # Capped at medium when the legs share upstream"), so the verdict itself
        # Must exist -- the first product-level pass found that the
        # Summary made a claim about an object nothing ever computed or showed.
        result.tree_recommendation = self.recommend_tree(result)

        return result

    def recommend_tree(self, result: PhylogeneticResult) -> TreeRecommendation:
        if not result.coalescent:
            return TreeRecommendation(
                recommended_tree=result.supermatrix.tree if result.supermatrix else None,
                confidence="medium",
                reason="Coalescent method not run (mode=off). Using concatenation tree.",
                warnings=["Coalescent inference was disabled"],
            )

        if not result.conflict_report:
            return TreeRecommendation(
                recommended_tree=result.supermatrix.tree or result.coalescent.species_tree,
                confidence="unknown",
                reason="Only one method was run",
            )

        rf = result.conflict_report.normalized_rf
        qa = getattr(result.conflict_report, "quartet_agreement", None)
        # Shared upstream caps recommendation confidence at
        # "medium" — "both frameworks agree" is weaker evidence when the
        # Both frameworks were fed the same gene trees.
        independence = getattr(result.conflict_report, "independence", None)
        shared_upstream = bool(
            independence is not None
            and (independence.shared_gene_tree_ratio > 0.0
                 or independence.ref_built_from_tested_markers)
        )

        def _metrics_str() -> str:
            # Decision basis shared in one place (consumed by both Interactive and
            # PlainText report generators via this recommendation): the normalized
            # RF distance and the quartet consistency, so the recommendation
            # Reason always states the conflict evidence, not just availability.
            parts = []
            if rf is not None:
                parts.append(f"normalized RF={rf:.3f}")
            else:
                parts.append("normalized RF=N/A")
            if qa is not None:
                parts.append(f"quartet agreement={qa:.3f}")
            return "; ".join(parts)

        if rf is None:
            # When conflict metrics are unmeasurable the
            # Pipeline must NOT hand out a tree as if nothing were wrong.
            return TreeRecommendation(
                recommended_tree=None,
                confidence="inconclusive",
                reason=(
                    f"Conflict metrics unavailable ({_metrics_str()}) — "
                    f"recommendation is inconclusive; re-run with a usable "
                    f"tree-comparison backend (ete3) or inspect both trees manually."
                ),
                warnings=["Conflict metrics unmeasurable: no tree recommended"],
            )
        # Recommendation bands come from configuration.
        rf_low = float(getattr(self.config, "recommend_rf_low", 0.1))
        rf_moderate = float(getattr(self.config, "recommend_rf_moderate", 0.3))
        if rf < rf_low:
            reason = f"Low conflict ({_metrics_str()}). Both methods agree."
            confidence = "high"
            if shared_upstream:
                confidence = "medium"
                reason += (
                    " Confidence capped to medium: the two legs share upstream "
                    "gene trees (cached reuse / reference built from the tested "
                    "markers), so their agreement is not fully independent evidence."
                )
            return TreeRecommendation(
                recommended_tree=result.supermatrix.tree,
                confidence=confidence,
                reason=reason,
            )
        elif rf < rf_moderate:
            return TreeRecommendation(
                recommended_tree=result.coalescent.species_tree,
                confidence="medium",
                reason=f"Moderate conflict ({_metrics_str()}). Coalescent method preferred.",
                warnings=["Some branches may be affected by ILS or HGT"],
            )
        else:
            return TreeRecommendation(
                recommended_tree=result.coalescent.species_tree,
                confidence="low",
                reason=f"Strong conflict ({_metrics_str()}). Significant HGT or method bias detected.",
                warnings=["High topological conflict between methods"],
            )
