from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple


class GenomeSource(Enum):
    ISOLATE = "isolate"
    MAG = "mag"
    SAG = "sag"


class GeneState(Enum):
    ABSENT = "absent"
    SINGLE_COPY = "single_copy"
    MULTI_COPY = "multi_copy"


@dataclass
class GenomeQuality:
    completeness: float = 0.0
    contamination: float = 0.0
    quality_score: float = 0.0
    strain_heterogeneity: Optional[float] = None
    n_contigs: Optional[int] = None
    genome_size: Optional[int] = None
    gc_content: Optional[float] = None
    coding_density: Optional[float] = None
    genome_id: Optional[str] = None


@dataclass
class Genome:
    id: str
    fasta_path: str
    source_type: GenomeSource = GenomeSource.ISOLATE
    mol_type: str = "protein"
    taxonomic_group: Optional[str] = None
    gc_content: Optional[float] = None
    gc_std: Optional[float] = None
    quality: Optional[GenomeQuality] = None
    # Optional paired sequence files. When both are provided, protein file
    # Drives HMM scans / marker extraction and nucleotide file drives
    # CheckM and composition-based HGT filtering.
    protein_fasta_path: Optional[str] = None
    nucleotide_fasta_path: Optional[str] = None
    _codon_cache: Optional[Dict[str, float]] = None
    _tetra_cache: Optional[List[float]] = None

    @property
    def effective_protein_path(self) -> str:
        """Return the protein FASTA path to use for HMM scans."""
        return self.protein_fasta_path or self.fasta_path

    @property
    def effective_nucleotide_path(self) -> Optional[str]:
        """Return the nucleotide FASTA path, if one is available."""
        return self.nucleotide_fasta_path

    def _sequence_path_for_composition(self) -> Optional[str]:
        """Choose the FASTA path used for GC/codon/tetranucleotide stats.

        Backward-compatible fallback chain:
        1. Explicit nucleotide file if paired input is available.
        2. Legacy ``fasta_path`` for nucleotide-only genomes or old callers.
        3. ``fasta_path`` as last resort to avoid breaking existing tests/code
           that pass DNA sequences through ``fasta_path``.
        """
        if self.nucleotide_fasta_path:
            return self.nucleotide_fasta_path
        return self.fasta_path

    def get_background_codon_usage(self) -> Dict[str, float]:
        """Compute codon-usage frequencies from this genome's sequence file.

.. deprecated:: 0.1.0
            The composition-based HGT step has been removed; this method is
            now in a **half-retired** state and has no in-pipeline consumer.
            It is kept only for backward-compatible ad-hoc analysis.

        NOTE: the underlying sequences MUST be aligned, equal-length protein
        (or coding) sequences. Codons are read frame-by-frame (step 3) across
        the concatenated sequence; un-aligned sequences of differing length
        shift the reading frame and yield meaningless codon statistics. No
        length check is performed here — callers must guarantee alignment.
        """
        if self._codon_cache is not None:
            return self._codon_cache
        codon_counts: Dict[str, int] = {}
        total = 0
        seq_path = self._sequence_path_for_composition()
        if seq_path is None:
            self._codon_cache = {}
            return {}
        try:
            with open(seq_path, encoding="utf-8", newline="") as f:
                seq_parts: List[str] = []
                for line in f:
                    if line.startswith(">"):
                        if seq_parts:
                            seq = "".join(seq_parts).upper()
                            for i in range(0, len(seq) - 2, 3):
                                codon = seq[i:i+3]
                                if len(codon) == 3 and "N" not in codon:
                                    codon_counts[codon] = codon_counts.get(codon, 0) + 1
                                    total += 1
                        seq_parts = []
                    else:
                        seq_parts.append(line.strip())
                if seq_parts:
                    seq = "".join(seq_parts).upper()
                    for i in range(0, len(seq) - 2, 3):
                        codon = seq[i:i+3]
                        if len(codon) == 3 and "N" not in codon:
                            codon_counts[codon] = codon_counts.get(codon, 0) + 1
                            total += 1
        except (FileNotFoundError, OSError):
            self._codon_cache = {}
            return {}
        if total == 0:
            self._codon_cache = {}
            return {}
        self._codon_cache = {k: v / total for k, v in codon_counts.items()}
        return self._codon_cache

    def get_tetranucleotide_signature(self) -> List[float]:
        """Compute the tetranucleotide (4-mer) composition signature.

.. deprecated:: 0.1.0
            The composition-based HGT step has been removed; this method is
            now in a **half-retired** state and has no in-pipeline consumer.
            It is kept only for backward-compatible ad-hoc analysis.

        NOTE: the input sequences MUST be aligned, equal-length nucleotide
        sequences and must be in the SAME coordinate system as any other genome
        being compared (the signature is unit-normalized per genome but the
        k-mer windows are positional). Computing this on raw, un-aligned
        sequences produces inconsistent signatures and corrupts MAG
        heterogeneity / strain assessment. No length check is performed here —
        callers must guarantee alignment before calling.
        """
        if self._tetra_cache is not None:
            return self._tetra_cache
        tetra_counts: Dict[str, int] = {}
        total = 0
        seq_path = self._sequence_path_for_composition()
        if seq_path is None:
            self._tetra_cache = []
            return []
        try:
            with open(seq_path, encoding="utf-8", newline="") as f:
                seq_parts: List[str] = []
                for line in f:
                    if line.startswith(">"):
                        if seq_parts:
                            seq = "".join(seq_parts).upper()
                            for i in range(len(seq) - 3):
                                tetra = seq[i:i+4]
                                if "N" not in tetra and len(tetra) == 4:
                                    tetra_counts[tetra] = tetra_counts.get(tetra, 0) + 1
                                    total += 1
                        seq_parts = []
                    else:
                        seq_parts.append(line.strip())
                if seq_parts:
                    seq = "".join(seq_parts).upper()
                    for i in range(len(seq) - 3):
                        tetra = seq[i:i+4]
                        if "N" not in tetra and len(tetra) == 4:
                            tetra_counts[tetra] = tetra_counts.get(tetra, 0) + 1
                            total += 1
        except (FileNotFoundError, OSError):
            self._tetra_cache = []
            return []
        if total == 0:
            self._tetra_cache = []
            return []
        sorted_keys = sorted(tetra_counts.keys())
        self._tetra_cache = [tetra_counts.get(k, 0) / total for k in sorted_keys]
        return self._tetra_cache


@dataclass
class OccupancyMatrix:
    """Per-genome × per-marker (COG) presence/absence state matrix.

    IMPORTANT: the sequences used to derive these states MUST be aligned,
    equal-length sequences. Occupancy (single/multi-copy vs absent) is computed
    per aligned column / per marker across genomes; feeding raw, un-aligned
    protein sequences of differing lengths distorts copy-number and occupancy
    estimates and biases downstream marker-quality scores. Callers are
    responsible for passing aligned input — this class performs no length
    reconciliation.
    """

    genomes: List[str] = field(default_factory=list)
    cogs: List[str] = field(default_factory=list)
    _data: Dict[Tuple[str, str], GeneState] = field(default_factory=dict)
    genome_qualities: Dict[str, float] = field(default_factory=dict)

    def get(self, genome: str, cog: str) -> GeneState:
        return self._data.get((genome, cog), GeneState.ABSENT)

    def set(self, genome: str, cog: str, state: GeneState) -> None:
        self._data[(genome, cog)] = state

    def get_column(self, cog: str) -> List[GeneState]:
        return [self.get(g, cog) for g in self.genomes]

    def get_quality(self, genome: str) -> float:
        return self.genome_qualities.get(genome, 1.0)

    def get_all_states(self) -> List[GeneState]:
        return list(self._data.values())

    @property
    def n_genomes(self) -> int:
        return len(self.genomes)

    @property
    def n_cogs(self) -> int:
        return len(self.cogs)

    def subset(self, selected_cogs: List[str], selected_genomes: List[str]) -> OccupancyMatrix:
        new_matrix = OccupancyMatrix(genomes=selected_genomes, cogs=selected_cogs)
        for g in selected_genomes:
            for c in selected_cogs:
                new_matrix.set(g, c, self.get(g, c))
        new_matrix.genome_qualities = {
            g: self.genome_qualities.get(g, 1.0) for g in selected_genomes
        }
        return new_matrix
