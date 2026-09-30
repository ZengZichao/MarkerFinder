"""Execute the ``report`` subcommand (Phase 4, report generation)."""

from __future__ import annotations

from markerfinder.cli.commands import _execute_pipeline


def execute(args, pipeline, tree, table_taxa, taxonomy_origin, logger) -> int:
    """Run the report step (HTML/text reports and metadata emission)."""
    return _execute_pipeline("report", args, pipeline, tree, table_taxa, taxonomy_origin, logger)
