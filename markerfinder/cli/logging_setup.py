"""Signal handling and log-level resolution for the MarkerFinder CLI."""

from __future__ import annotations

import logging
import signal
import sys

from markerfinder.cli.constants import EXIT_INTERRUPT

logger = logging.getLogger(__name__)

_interrupted = False


def _signal_handler(signum: int, frame: object) -> None:
    global _interrupted
    if _interrupted:
        sys.stderr.write("\nForced exit.\n")
        sys.exit(EXIT_INTERRUPT)
    _interrupted = True
    logger.warning("Interrupt received, shutting down gracefully...")
    sys.exit(EXIT_INTERRUPT)


def _register_signal_handlers() -> None:
    try:
        signal.signal(signal.SIGINT, _signal_handler)
        signal.signal(signal.SIGTERM, _signal_handler)
    except (OSError, ValueError):
        pass


def resolve_log_level(args) -> int:
    """Map parsed ``-v`` verbosity onto a logging level.

    Mirrors the historical behaviour in ``__main__.main``: 0 → WARNING,
    1 → INFO, 2+ → DEBUG.
    """
    verbose = getattr(args, "verbose", 0)
    if verbose >= 2:
        return logging.DEBUG
    if verbose >= 1:
        return logging.INFO
    return logging.WARNING
