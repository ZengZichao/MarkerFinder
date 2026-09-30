"""Execute the ``scan`` subcommand (Phase 0+1+1.5)."""

from __future__ import annotations

from markerfinder.cli.commands import _execute_pipeline


def execute(args, pipeline, tree, table_taxa, taxonomy_origin, logger) -> int:
    """Run the scan step (quality preprocessing + marker scanning + extraction)."""
    return _execute_pipeline("scan", args, pipeline, tree, table_taxa, taxonomy_origin, logger)
