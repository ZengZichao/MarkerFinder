"""Execute the ``infer`` subcommand (Phase 3, phylogenetic inference)."""

from __future__ import annotations

from markerfinder.cli.commands import _execute_pipeline


def execute(args, pipeline, tree, table_taxa, taxonomy_origin, logger) -> int:
    """Run the infer step (supermatrix + coalescent species-tree inference)."""
    return _execute_pipeline("infer", args, pipeline, tree, table_taxa, taxonomy_origin, logger)
