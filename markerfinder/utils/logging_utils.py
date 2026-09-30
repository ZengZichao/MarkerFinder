"""Logging utilities: ISO 8601 timestamps, flush, file persistence.

Log format: 2025-03-21T10:15:30.123 | INFO | message
No color output. Console and file have identical format.
File handler includes extra debug info (stack traces on CRITICAL).
"""

from __future__ import annotations

import logging
import os
import sys
from typing import Optional


class ISO8601Formatter(logging.Formatter):
    """Custom formatter with ISO 8601 timestamps to millisecond precision."""

    def formatTime(self, record: logging.LogRecord, datefmt: Optional[str] = None) -> str:
        import datetime
        dt = datetime.datetime.fromtimestamp(record.created, tz=datetime.timezone.utc)
        return dt.strftime("%Y-%m-%dT%H:%M:%S.") + f"{record.msecs:03.0f}"

    def format(self, record: logging.LogRecord) -> str:
        timestamp = self.formatTime(record)
        level = f"{record.levelname:<8s}"
        module_func = f"{record.module}/{record.funcName}"
        msg = record.getMessage()

        formatted = f"{timestamp} | {level} | [{module_func}] {msg}"

        if record.exc_info and record.exc_info[0] is not None:
            formatted += "\n" + self.formatException(record.exc_info)

        return formatted


class FlushStreamHandler(logging.StreamHandler):
    """Stream handler that flushes after every emit."""

    def emit(self, record: logging.LogRecord) -> None:
        super().emit(record)
        self.flush()


class UTF8FileHandler(logging.FileHandler):
    """File handler with explicit UT encoding and flush-on-emit."""

    def __init__(self, filename: str, mode: str = "a") -> None:
        super().__init__(filename, mode=mode, encoding="utf-8")

    def emit(self, record: logging.LogRecord) -> None:
        super().emit(record)
        self.flush()


def setup_logging(
    level: int = logging.INFO,
    log_file: Optional[str] = None,
) -> logging.Logger:
    """Configure the root logger with ISO 8601 timestamps and flush.

    Args:
        level: logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        log_file: optional path to also write logs to a UT file.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    formatter = ISO8601Formatter()

    console_handler = FlushStreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if log_file:
        # Ensure the parent directory exists; the output directory is created
        # By the pipeline, which runs after logging is configured.
        log_dir = os.path.dirname(log_file)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        file_handler = UTF8FileHandler(log_file, mode="a")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    return root_logger


def get_logger(name: str, level: Optional[int] = None) -> logging.Logger:
    """Get a named logger (does not reconfigure handlers)."""
    logger = logging.getLogger(name)
    if level is not None:
        logger.setLevel(level)
    return logger
