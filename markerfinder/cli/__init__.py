"""MarkerFinder command-line interface package.

Re-exports the argument parser builder and the CLI ``main`` entry point so that
``markerfinder.__main__`` (and the ``markerfinder`` console script declared in
``pyproject.toml``) can import them directly.
"""

from __future__ import annotations

from markerfinder.cli.parser import _build_parser as build_parser
from markerfinder.cli.main import main

__all__ = ["build_parser", "main"]
