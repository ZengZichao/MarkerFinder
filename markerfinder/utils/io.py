"""File I/O utilities."""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Optional

from markerfinder.models.genome import Genome, GenomeSource

logger = logging.getLogger(__name__)


def _collect_files(input_path: Path, file_ext: Optional[str] = None) -> List[Path]:
    """Collect candidate genome FASTA files, ordered by stem.

    When ``file_ext`` is None the preferred order is:
      1. ``*.faa`` (protein)
      2. ``*.fna`` (nucleotide)
      3. ``*.fasta`` / ``*.fa`` fallback
    """
    if file_ext:
        suffix = file_ext if file_ext.startswith(".") else f".{file_ext}"
        files = sorted(input_path.glob(f"*{suffix}"))
        if not files:
            files = sorted(input_path.glob(f"*{file_ext.replace('.', '')}*"))
            if not files:
                files = sorted(input_path.glob("*.fasta")) + sorted(input_path.glob("*.fa"))
        return files

    faa_files = sorted(input_path.glob("*.faa"))
    fna_files = sorted(input_path.glob("*.fna"))
    if faa_files or fna_files:
        # De-duplicate by stem while preserving faa-first order.
        seen = set()
        files = []
        for p in faa_files + fna_files:
            if p.stem not in seen:
                files.append(p)
                seen.add(p.stem)
        files = _drop_empty_files(files)
        _report_ignored_files(input_path, files)
        return files

    files = sorted(input_path.glob("*.fasta")) + sorted(input_path.glob("*.fa"))
    files = _drop_empty_files(files)
    _report_ignored_files(input_path, files)
    return files


def _drop_empty_files(files: List[Path]) -> List[Path]:
    """A zero-length FASTA is not a genome.

    Left in, it enters the tip set with no sequence, contributes nothing to any
    marker, and dilutes every per-marker occupancy fraction computed against the
    genome count — so a user's typo (``touch bin_12.faa``) quietly changes the
    statistics. It is dropped, and the drop is reported.
    """
    keep, empty = [], []
    for p in files:
        try:
            (keep if p.stat().st_size else empty).append(p)
        except OSError:      # Unreadable: let the downstream reader report it
            keep.append(p)
    if empty:
        logger.warning(
            f"  Input directory: skipping {len(empty)} empty FASTA file(s): "
            f"{', '.join(p.name for p in empty)}"
        )
    return keep


def _report_ignored_files(input_path: Path, used: List[Path]) -> None:
    """Name every file in the input directory that is not being read.

    A directory that also holds ``README.md`` or a spreadsheet used to be
    processed as if it contained only genomes. The names are what a reader needs
    to notice a missing genome — an input that is silently not an input is the
    beginning of a wrong result, and the genome count printed in the summary
    otherwise disagrees with the directory.
    """
    used_names = {p.name for p in used}
    ignored = sorted(
        p.name for p in input_path.iterdir()
        if p.is_file() and p.name not in used_names
    )
    if ignored:
        logger.warning(
            f"  Input directory {input_path}: {len(ignored)} file(s) not read "
            f"(no recognised FASTA extension): {', '.join(ignored)}"
        )


def _mol_type_for_file(path: Path) -> str:
    if path.suffix.lower() == ".fna":
        return "nucleotide"
    if path.suffix.lower() in (".faa", ".fasta", ".fa"):
        return "protein"
    return "protein"


def iter_genome_ids(input_dir: str, file_ext: Optional[str] = None) -> List[str]:
    """List genome IDs (stem of the filename) from a directory, without loading sequences.

    Mirrors the ordering and format rules of:func:`load_genomes_from_directory`
    so that downstream ID lookups see the same set of genomes in the same order.
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_path.is_dir():
        raise ValueError(f"Expected a directory, got file: {input_dir}")

    files = _collect_files(input_path, file_ext=file_ext)
    return [p.stem for p in files]


def load_genomes_from_directory(input_dir: str, file_ext: Optional[str] = None) -> List[Genome]:
    """Load genome sequence files from a directory.

    Supports three input layouts:

    1. Paired protein + nucleotide (recommended):
       Directory contains ``*.faa`` and matching ``*.fna`` files with the same
       stem. Protein sequences drive HMM scans and marker extraction;
       nucleotide sequences drive CheckM and composition-based HGT filtering.

    2. Protein-only:
       Directory contains only ``*.faa`` files. HMM scans work normally;
       CheckM and DNA composition steps fall back to defaults / skip.

    3. Nucleotide-only:
       Directory contains only ``*.fna`` files. This is allowed as a fallback,
       but HMM scans (protein HMMs against DNA) may fail or produce no hits;
       a warning is logged.

    Args:
        input_dir: path to directory containing genome FASTA files
        file_ext: explicit file extension filter; if None, auto-detect

    Returns:
        List of Genome objects with mol_type and path fields set accordingly
    """
    genomes: List[Genome] = []
    input_path = Path(input_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    if input_path.is_file():
        raise ValueError(f"Expected a directory, got file: {input_dir}")

    files = _collect_files(input_path, file_ext=file_ext)

    if file_ext is None:
        faa_by_stem = {p.stem: p for p in input_path.glob("*.faa")}
        fna_by_stem = {p.stem: p for p in input_path.glob("*.fna")}
    else:
        faa_by_stem: dict = {}
        fna_by_stem: dict = {}

    for fa_file in files:
        genome_id = fa_file.stem
        source_type = _infer_source_type(genome_id)

        if file_ext is not None:
            mol_type = _mol_type_for_file(fa_file)
            genome = Genome(
                id=genome_id,
                fasta_path=str(fa_file),
                source_type=source_type,
                mol_type=mol_type,
            )
        else:
            faa = faa_by_stem.get(genome_id)
            fna = fna_by_stem.get(genome_id)

            if faa and fna:
                mol_type = "protein"
                genome = Genome(
                    id=genome_id,
                    fasta_path=str(faa),
                    source_type=source_type,
                    mol_type=mol_type,
                    protein_fasta_path=str(faa),
                    nucleotide_fasta_path=str(fna),
                )
            elif faa:
                mol_type = "protein"
                genome = Genome(
                    id=genome_id,
                    fasta_path=str(faa),
                    source_type=source_type,
                    mol_type=mol_type,
                )
            elif fna:
                mol_type = "nucleotide"
                genome = Genome(
                    id=genome_id,
                    fasta_path=str(fna),
                    source_type=source_type,
                    mol_type=mol_type,
                )
                logger.warning(
                    f"  Nucleotide-only input for {genome_id}: HMM scans may fail "
                    f"because protein HMMs are searched against DNA sequences."
                )
            else:
                # Fallback for.fasta/.fa files
                mol_type = "protein"
                genome = Genome(
                    id=genome_id,
                    fasta_path=str(fa_file),
                    source_type=source_type,
                    mol_type=mol_type,
                )

        genomes.append(genome)

    logger.info(f"Loaded {len(genomes)} genomes from {input_dir}")
    return genomes


# Delimiter-aware match so substrings such as 'mag' inside 'magellan' or
# 'bin' inside 'mybin' do not trigger a false MAG/SAG classification. The
# Keyword must be bounded by a non-letter character (or a string edge), which
# Matches real layouts like 'MAG_001', 'genome_bin.15', or 'SAG_xyz' while
# Rejecting 'magellan'.
_MAG_BIN_RE = re.compile(r"(?:^|[^a-z])(mag|bin)(?:[^a-z]|$)", re.IGNORECASE)
_SAG_RE = re.compile(r"(?:^|[^a-z])sag(?:[^a-z]|$)", re.IGNORECASE)


def _infer_source_type(genome_id: str) -> GenomeSource:
    """Infer genome source type from filename.

    Uses a delimiter-aware match so a substring such as 'mag' inside 'magellan'
    or 'bin' inside 'mybin' does not produce a false MAG/SAG classification.
    MAG/BIN takes precedence over SAG, matching the previous substring priority.
    """
    if _MAG_BIN_RE.search(genome_id):
        return GenomeSource.MAG
    if _SAG_RE.search(genome_id):
        return GenomeSource.SAG
    return GenomeSource.ISOLATE


def write_fasta(sequences: dict, output_path: str) -> str:
    """Write sequences in FASTA format.

    Args:
        sequences: {id: sequence_string} dictionary
        output_path: output file path
    """
    lines: List[str] = []
    for seq_id, seq_str in sequences.items():
        lines.append(f">{seq_id}")
        for i in range(0, len(seq_str), 80):
            lines.append(seq_str[i : i + 80])
    target = Path(output_path).resolve()
    target.write_text("\n".join(lines) + "\n" if lines else "", encoding="utf-8", newline="\n")
    return output_path


def ensure_dir(path: str) -> str:
    """Ensure directory exists."""
    os.makedirs(path, exist_ok=True)
    return path


def parse_fasta(path: str) -> Dict[str, str]:
    """Parse a FASTA file into ``{header_first_token: sequence}``.

    The header's first whitespace-delimited token is used as the sequence id.
    This matches the legacy ``_parse_fasta`` behaviour previously duplicated in
    ``marker_selection`` and ``gtdb_tk_markers``. Returns an empty dict when the
    file is missing (a warning is logged).

    重复 ID 不再静默覆盖: 每次覆盖都会打 WARNING(保留后出现的序列,
    与既有 last-wins 行为一致), 避免上游文件混入重复行时无声丢数据.
    """
    seqs: Dict[str, str] = {}
    cur_id = ""
    cur_seq: List[str] = []
    try:
        with open(path, encoding="utf-8", newline="") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_id:
                        if cur_id in seqs:
                            logger.warning(
                                f"Duplicate FASTA id '{cur_id}' in {path}; "
                                f"keeping the later occurrence"
                            )
                        seqs[cur_id] = "".join(cur_seq)
                    header = line[1:]
                    cur_id = header.split()[0] if header.split() else header
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_id:
                if cur_id in seqs:
                    logger.warning(
                        f"Duplicate FASTA id '{cur_id}' in {path}; "
                        f"keeping the later occurrence"
                    )
                seqs[cur_id] = "".join(cur_seq)
    except FileNotFoundError:
        logger.warning(f"FASTA file not found: {path}")
    return seqs
