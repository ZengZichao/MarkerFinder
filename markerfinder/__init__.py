"""MarkerFinder: Adaptive HGT-Aware Phylogenomic Pipeline.

Author: Zengzichao (zengzichao@sjtu.edu.cn, ORCID 0000-0001-6553-970X)

Public API:
    parse_taxonomy(label, mode="reverse") -> dict
    is_monophyletic(tree_newick, taxon_label, taxonomy_map, level=None) -> bool
    merge_taxonomy(embedded, from_table, priority="table") -> dict
    load_taxonomy_table(file_path, sep=None, taxonomy_format="table", delimiter_mode="reverse") -> dict
    find_special_identifier(tree_newick, taxonomy_map, identifier) -> node
    load_tree(path, validate=True) -> Tree
    cross_validate(tree_path, seq_path, strict=True) -> ValidationReport
"""

from markerfinder._version import (
    __author__,
    __email__,
    __orcid__,
    __version__,
    get_version_string,
)

# Must precede any ete3 import: Python 3.13 removed the stdlib ``cgi`` module
# that ete3 imports at package-import time. Importing the package installs a
# minimal stand-in so that `import ete3` works on 3.13+ for every caller
# (including third-party and test code that imports ete3 directly).
from markerfinder import _cgi_compat as _cgi_compat  # noqa: F401

from markerfinder.taxonomy import (
    parse_taxonomy,
    is_monophyletic,
    merge_taxonomy,
    load_taxonomy_table,
    find_special_identifier,
    parse_custom_levels,
    get_merged_level_maps,
)
from markerfinder.validation import load_tree, cross_validate, ValidationReport, check_control_characters, validate_string_safety
from markerfinder.exceptions import (
    PhyloToolError,
    PhyloFormatError,
    TaxonomyConflictError,
    MonophylyError,
    TreeValidationError,
    SequenceValidationError,
    CrossValidationError,
    InputError,
    ConfigError,
    MultiTreeError,
)

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__orcid__",
    "get_version_string",
    "parse_taxonomy",
    "is_monophyletic",
    "merge_taxonomy",
    "load_taxonomy_table",
    "find_special_identifier",
    "parse_custom_levels",
    "get_merged_level_maps",
    "load_tree",
    "cross_validate",
    "ValidationReport",
    "check_control_characters",
    "validate_string_safety",
    "PhyloToolError",
    "PhyloFormatError",
    "TaxonomyConflictError",
    "MonophylyError",
    "TreeValidationError",
    "SequenceValidationError",
    "CrossValidationError",
    "InputError",
    "ConfigError",
    "MultiTreeError",
]
