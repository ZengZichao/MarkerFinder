"""MarkerFinder CLI entry point.

Exit codes:
  0: Success
  1: Runtime error (missing dependency, pipeline failure)
  2: Argument or configuration error
  3: Input data format or content error
  130: User interrupt (SIGINT/SIGTERM)

This module is now a thin dispatch shim. The argument parser and the ``main``
orchestrator live in:mod:`markerfinder.cli`; the step-state helpers that the
unit tests import are re-exported here so the public import contract
(``from markerfinder.__main__ import main`` and the four ``_handle_step_state``
/ ``_load_step_state`` / ``_step_status`` / ``_truncate_state_for_redo``
helpers) is preserved exactly.
"""

from markerfinder.cli import build_parser, main
from markerfinder.cli.step_state import (
    _handle_step_state,
    _load_step_state,
    _step_status,
    _truncate_state_for_redo,
)
from markerfinder.cli.constants import (
    EXIT_SUCCESS,
    EXIT_RUNTIME_ERROR,
    EXIT_ARG_ERROR,
    EXIT_DATA_ERROR,
    EXIT_INTERRUPT,
)

if __name__ == "__main__":
    import sys
    sys.exit(main())
