from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Tuple

from markerfinder.exceptions import PhyloToolError
from markerfinder.models.genome import GenomeQuality
from markerfinder.models.evidence import MeasureState, UnknownReason
from markerfinder.models.marker import MarkerLevel


class ConflictType(Enum):
    HGT_SIGNAL = "hgt"
    AMBIGUOUS = "ambiguous"
    # ILS_SIGNAL / METHOD_BIAS removed by: declared with no
    # Assignment point anywhere in the package. Keeping unimplemented conflict
    # Categories is the same shape of unmet promise as the paper's "reserved
    # Interface". HGT_SIGNAL stays because assigns it; AMBIGUOUS stays as
    # The dataclass default.


class NoMarkerAvailableError(PhyloToolError):
    pass


class PipelineError(PhyloToolError):
    pass


class PhylogeneticInferenceError(PhyloToolError):
    pass


@dataclass
class PipelineParameters:
    pass


@dataclass
class BlastHit:
    bitscore: float = 0.0
    evalue: float = 1.0
    species: str = ""
    is_ingroup: bool = False
    taxonomic_distance: float = 0.0


@dataclass
class PhyloBranch:
    branch_id: str = ""
    left_clade: Set[str] = field(default_factory=set)
    right_clade: Set[str] = field(default_factory=set)
    support: float = 0.0


@dataclass
class PhyloStepResult:
    """Output of the phylogenetic (phylogenetic/protein-tree) HGT screening step.

    Carries the per-marker phylogenetic HGT risk from the phylogenetic
    (phylogenetic-per-marker) screen. This is the risk that a marker's
    protein tree disagrees with the species tree or taxonomic hierarchy,
    a hallmark of horizontal transfer.

    Two available sub-strategies:
      - MAD-rooted monophyly-proportion screen (taxonomy-table path): fraction of
        the informative taxa at the measured rank whose tips the gene tree
        groups as one side of a split. Root-invariant by construction: a taxon
        counts when the taxon OR its complement is a clade, so the verdict does
        not depend on where the rooting rule placed the root. ``None`` when no
        taxon at the rank has an informative split (>=2 representatives on both
        sides of the split).
      - RF + quartet-consistency screen (species-tree path): RF distance and
        quartet consistency between the marker's gene tree and the provided
        species tree.

    Under ``monophyly_rank='auto'`` (the default) the measured rank is derived
    from the tree: one below the taxonomic scope of its tips. Naming a rank
    measures at that rank instead, and the scope is only reported.
    """
    gene_id: str = ""
    # A component that was not measured is ``None`` —
    # Never a numeric placeholder (0.0/0.5/1.0).
    rf_distance: Optional[int] = None
    normalized_rf: Optional[float] = None
    quartet_consistency: Optional[float] = None
    n_conflicting_branches: int = 0
    conflicting_branches: List[ClassifiedConflict] = field(default_factory=list)
    # Sh_test_pvalue removed by: declared but never assigned —
    # An unfulfilled promise rather than a feature.
    monophyly_proportion: Optional[float] = None
    overall_risk: float = 0.0
    is_suspicious: bool = False
    # Per-component evidence states. Default MEASURED keeps
    # Constructing sites that already fill real values fully compatible.
    rf_state: MeasureState = MeasureState.MEASURED
    quartet_state: MeasureState = MeasureState.MEASURED
    monophyly_state: MeasureState = MeasureState.NOT_APPLICABLE
    rf_reason: str = ""
    quartet_reason: str = ""
    # "two_signal_weighted" | "monophyly_only" | "unscreened".
    risk_basis: str = ""
    # The monophyly path records which rank it actually
    # Tested and the taxon counts behind the proportion.
    test_level_used: Optional[str] = None
    n_total: Optional[int] = None
    n_mono: Optional[int] = None
    # The measured proportion came from a different rank than
    # The configured one — cross-rank proportions must not share a threshold.
    cross_rank_comparison: bool = False

    def unmeasured_components(self) -> List["UnknownReason"]:
        """Critical components of this result that were not measured.

        A REJECTED state (reference illegality) dominates and maps to
        REFERENCE_ILLEGAL. The monophyly path (monophyly measured) is complete
        by itself; the species-tree path requires BOTH rf and quartet.
        """
        if MeasureState.REJECTED in (self.rf_state, self.quartet_state):
            return [UnknownReason.REFERENCE_ILLEGAL]
        if self.monophyly_state is MeasureState.MEASURED:
            return []
        unmet: List[UnknownReason] = []
        if self.rf_state is not MeasureState.MEASURED:
            unmet.append(UnknownReason.RF_UNMEASURABLE)
        if self.quartet_state is not MeasureState.MEASURED:
            unmet.append(UnknownReason.QUARTET_UNMEASURABLE)
        return unmet

    def reason_for(self, reason: "UnknownReason") -> str:
        """Human-readable detail for an unmeasured component."""
        if reason is UnknownReason.RF_UNMEASURABLE:
            return self.rf_reason or "RF distance was not measurable"
        if reason is UnknownReason.QUARTET_UNMEASURABLE:
            return self.quartet_reason or "quartet consistency was not measurable"
        return ""

    def measured_scores(self) -> Dict[str, float]:
        """Real signal values available for weighted synthesis."""
        scores: Dict[str, float] = {}
        if self.rf_state is MeasureState.MEASURED and self.normalized_rf is not None:
            scores["rf"] = self.normalized_rf
        if self.quartet_state is MeasureState.MEASURED and self.quartet_consistency is not None:
            scores["quartet"] = 1.0 - self.quartet_consistency
        return scores


@dataclass
class ConflictingBranch:
    branch: PhyloBranch = field(default_factory=PhyloBranch)
    conflict_type: ConflictType = ConflictType.AMBIGUOUS
    concat_support: float = 0.0
    astral_support: float = 0.0
    split_support: float = 0.0


@dataclass
class QuartetConflict:
    species: Tuple[str, ...] = ()
    topology_tree1: str = ""
    topology_tree2: str = ""
    conflict_type: ConflictType = ConflictType.AMBIGUOUS


@dataclass
class GeneTreeAgreement:
    support_for_concat: Dict[str, float] = field(default_factory=dict)
    support_for_astral: Dict[str, float] = field(default_factory=dict)
    split_support: Dict[str, float] = field(default_factory=dict)


@dataclass
class ClassifiedConflict:
    branch: PhyloBranch = field(default_factory=PhyloBranch)
    conflict_type: ConflictType = ConflictType.AMBIGUOUS
    concat_support: float = 0.0
    astral_support: float = 0.0
    split_support: float = 0.0


@dataclass
class HGTEvaluation:
    marker_id: str = ""
    overall_risk: float = 0.0
    level: MarkerLevel = MarkerLevel.LEVEL_1
    confidence: str = "low"
    step_details: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    notes: str = ""
    # How the risk was constituted —
    # "two_signal_weighted" | "monophyly_only" | "unscreened".
    risk_basis: str = ""
    # The originating PhyloStepResult, kept for the
    # Assertion layer so it need not re-read logs. None for UNKNOWN.
    phylo_step: Optional[PhyloStepResult] = None
    # Per-marker provenance card. Machine-readable
    # Mirror of the TSV columns; schema_version 2.
    decision_card: Dict[str, object] = field(default_factory=dict)


@dataclass
class HGTReport:
    total_markers: int = 0
    level1_count: int = 0
    level2_count: int = 0
    level3_count: int = 0
    # HGT 筛查未能运行(无参照树/无 taxonomy 表/基因树不可得)时的 UNKNOWN 计数.
    # UNKNOWN 标记被保留用于建树, 但其 HGT 状态未经过筛查, 报告需单独呈现.
    unknown_count: int = 0
    mean_risk: float = 0.0
    high_confidence_count: int = 0
    marker_evaluations: List[HGTEvaluation] = field(default_factory=list)
    phylogenetic_enabled: bool = True
    # Evidence coverage — markers whose HGT state was truly
    # Measured (graded, not UNKNOWN). Rendered at the top of
    # Pipeline_summary.txt; below ``require_evidence_coverage`` triggers a
    # Prominent warning.
    n_actually_measured: int = 0
    evidence_coverage: float = 0.0
    unknown_reason_counts: Dict[str, int] = field(default_factory=dict)
    # Far-distance mode must be visible at report top level.
    far_active: bool = False


@dataclass
class SupermatrixResult:
    alignment: Optional[object] = None
    partition: Optional[object] = None
    tree: Optional[object] = None
    marker_order: List[str] = field(default_factory=list)
    best_model: Optional[str] = None
    avg_ufboot: Optional[float] = None


@dataclass
class CoalescentResult:
    species_tree: Optional[object] = None
    gene_trees: Dict[str, object] = field(default_factory=dict)
    filtered_gene_trees: Dict[str, object] = field(default_factory=dict)
    n_total_genes: int = 0
    n_passed_filter: int = 0
    avg_quartet_support: Optional[float] = None
    # {marker_id: "cached"|"rebuilt"} gene-tree provenance.
    gene_tree_provenance: Dict[str, str] = field(default_factory=dict)
    # 物种树来源标记: 当 ASTRAL-III 缺失 / 失败时回退到 consensus / first_gene_tree.
    # 取值: 'astral' | 'consensus' | 'first_gene_tree' | 'none'.
    species_tree_source: str = "none"


@dataclass
class IndependenceReport:
    """How independent the two legs' evidence actually is.

    ``shared_gene_tree_ratio``: fraction of coalescent gene trees reused from
    the concatenation/HGT cache. ``same_trimming_regime``: whether both legs
    trimmed alignments identically (today: they never do — documents three
    regimes). ``ref_built_from_tested_markers``: the GTDB-TK default builds the
    reference from the very markers under test — almost always True, which is
    an honest fact, not a bug.
    """
    marker_set_jaccard: float = 0.0
    shared_gene_tree_ratio: float = 0.0
    same_trimming_regime: bool = False
    ref_built_from_tested_markers: bool = True


@dataclass
class ConflictReport:
    """Topology disagreement between the concatenation and coalescent trees.

    Currently computed fields: ``rf_distance``, ``normalized_rf``,
    ``quartet_agreement``, and ``topology_note``.

    The following fields are reserved placeholders for future branch-level
    conflict classification and are **not calculated** by the current
    ``ConflictDetector`` implementation:
      - ``n_quartet_conflicts`` / ``quartet_conflicts``
      - ``n_conflicting_branches`` / ``conflicting_branches``
      - ``gene_tree_agreement``
      - ``conflict_summary``

    They are retained to keep the public dataclass stable for downstream
    reports and existing tests.
    """
    rf_distance: Optional[int] = None
    normalized_rf: Optional[float] = None
    n_quartet_conflicts: int = 0
    quartet_conflicts: List[QuartetConflict] = field(default_factory=list)
    n_conflicting_branches: int = 0
    conflicting_branches: List[ClassifiedConflict] = field(default_factory=list)
    gene_tree_agreement: Optional[GeneTreeAgreement] = None
    conflict_summary: Dict[str, int] = field(default_factory=dict)
    quartet_agreement: Optional[float] = None
    topology_note: str = ""
    #
    independence: Optional[IndependenceReport] = None


@dataclass
class PhylogeneticResult:
    supermatrix: Optional[SupermatrixResult] = None
    coalescent: Optional[CoalescentResult] = None
    conflict_report: Optional[ConflictReport] = None
    # The dataset-level tree verdict computed by
    # ``recommend_tree`` inside the run. Typed loosely (``object``) to keep
    # Models free of an import cycle, same as ``marker_set`` below -- the
    # Report renders its confidence/reason, and the independence block's
    # "capped at medium" claim must agree with THIS object, not dangle.
    tree_recommendation: Optional[object] = None
    # {marker_id: ConsistencyGrade.value}. Filled by the run
    # Only when --hgt-mode consistency; empty under the default risk mode, in
    # Which case the report renders the column as NA rather than inventing a
    # Grade. Storing plain strings keeps models free of a module cycle.
    consistency_grades: Dict[str, str] = field(default_factory=dict)


@dataclass
class MarkerSelectionResult:
    marker_set: Optional[object] = None
    quality_scores: Dict[str, object] = field(default_factory=dict)
    occupancy_matrix: Optional[object] = None
    resolved_sequences: Dict[str, Dict[str, object]] = field(default_factory=dict)


@dataclass
class QualityData:
    genome_qualities: List[GenomeQuality] = field(default_factory=list)
    reference_tree: Optional[object] = None
    # 质量数据来源. 取值:
    # 'checkm' -> 真实运行 CheckM 得到的质量值;
    # 'default_no_nucleotide' -> 蛋白-only 输入回退(无.fna);
    # 'default_no_checkm' -> CheckM 未安装/运行失败/无质量文件,使用默认估计;
    # 'precomputed' -> 用户通过 --checkm-results 提供的结果文件.
    quality_source: str = "checkm"


@dataclass
class AdaptiveParams:
    min_hmm_score: float = 20.0
    min_occupancy: float = 0.75
    max_markers: int = 60
    missing_data_strategy: str = "strict_gap"
    quality_weighted: bool = True
    allow_partial_hits: bool = False


@dataclass
class PreprocessingResult:
    quality_results: Optional[QualityData] = None
    adaptive_params: Optional[AdaptiveParams] = None
    layers: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class OptimizedMatrix:
    matrix: Optional[object] = None
    selected_markers: List[str] = field(default_factory=list)
    selected_species: List[str] = field(default_factory=list)
    effective_information: float = 0.0
    strategy: str = "sparse_optimized"


@dataclass
class WeightedAlignment:
    alignment: Optional[object] = None
    weight_file: str = ""
    weights: Dict[str, float] = field(default_factory=dict)


@dataclass
class Allele:
    id: str = ""
    representative: Optional[object] = None
    frequency: float = 0.0
    members: List[str] = field(default_factory=list)


@dataclass
class HeterogeneityResult:
    selected_sequence: Optional[object] = None
    heterogeneity_score: float = 0.0
    n_alleles: int = 0
    allele_frequencies: Dict[str, float] = field(default_factory=dict)
    confidence: str = "high"


# MissingPatternReport removed by: declared but never assigned.
