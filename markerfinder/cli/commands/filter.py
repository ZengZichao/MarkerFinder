"""Execute the ``filter`` subcommand (Phase 2, HGT-aware marker filtering)."""

from __future__ import annotations

from markerfinder.cli.commands import _execute_pipeline


def execute(args, pipeline, tree, table_taxa, taxonomy_origin, logger) -> int:
    """Run the filter step (phylogenetic HGT screening + level partitioning)."""
    return _execute_pipeline("filter", args, pipeline, tree, table_taxa, taxonomy_origin, logger)
