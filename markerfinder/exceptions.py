"""Custom exception hierarchy for MarkerFinder.

All public exceptions inherit from PhyloToolError.
"""


class PhyloToolError(Exception):
    """Base exception for all MarkerFinder errors."""


class PhyloToolUnavailable(PhyloToolError):
    """A phylogenetic Python dependency (ete3,...) is not importable in the
    current interpreter.

    Raised by:func:`markerfinder.utils.etree.require_ete3` (the single import
    entry point). Callers must either fall back to an explicitly
    declared pure-Python measurement path or surface this as "not measured"
    — never silently substitute a numeric placeholder.
    """


class UnsupportedCriterion(PhyloToolError):
    """A decision criterion was requested whose prerequisites are not met.

    Used by the marker-level consistency screen gate (
    ): ``--hgt-mode consistency|hybrid`` refuses to start when the
    required prerequisite work packages' observable artifacts are absent.
    """


class AssertionFailureError(PhyloToolError):
    """A severity=fail output assertion fired at end of run.

    Mapped to exit code 4 (EXIT_ASSERTION_FAILED) by the CLI so "output
    unreasonable" is distinguishable from a runtime error. Suppressed only
    by ``--allow-assertion-failure``, and the bypass is written to the report.
    """


class ExternalToolError(PhyloToolError):
    """An external bioinformatics tool (MAFFT, trimAl, FastTree, IQ-TREE3,
    ASTRAL-III, HMMER/hmmsearch, CheckM,...) failed or was not found on PATH.

    Raising this (instead of silently returning an empty/placeholder result)
    prevents a tool failure from being masqueraded as a valid biological
    conclusion downstream.
    """


class PhyloFormatError(PhyloToolError):
    """Invalid phylogenetic file format (Newick, Nexus, NHX, FASTA, etc.)."""


class TaxonomyConflictError(PhyloToolError):
    """Conflicting taxonomy information from multiple sources."""


class MonophylyError(PhyloToolError):
    """Error during monophyly check (e.g., taxon not found on tree)."""


class TreeValidationError(PhyloToolError):
    """Tree structure validation failure (negative branch, duplicate tips, etc.)."""


class SequenceValidationError(PhyloToolError):
    """Sequence file validation failure (invalid alphabet, duplicate IDs, etc.)."""


class CrossValidationError(PhyloToolError):
    """Mismatch between tree tip labels and sequence IDs."""


class InputError(PhyloToolError):
    """Invalid input (missing file, wrong path type, empty file)."""


class ConfigError(PhyloToolError):
    """Configuration file error (parse failure, invalid key)."""


class MultiTreeError(PhyloToolError):
    """Multiple trees detected in input with no handling mode specified."""
