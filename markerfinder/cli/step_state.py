"""Step-state helpers for the pipeline subcommands.

Verbatim relocation of the state-inspection/redo logic that used to live in
``markerfinder.__main__``. The four public helpers
(``_load_step_state``, ``_truncate_state_for_redo``, ``_step_status``,
``_handle_step_state``) are re-exported from ``markerfinder.__main__`` so the
existing unit tests continue to import them unchanged.
"""

from __future__ import annotations

import logging
import sys
from typing import Dict, Optional, Tuple

from markerfinder.utils.state_codec import (
    StateSchemaError,
    load_pipeline_state,
    truncate_pipeline_state,
)
from markerfinder.cli.constants import EXIT_SUCCESS, EXIT_ARG_ERROR

logger = logging.getLogger(__name__)


def _load_step_state(output_dir: str) -> Optional[dict]:
    """Load persisted step state, returning None if missing.

    Uses the explicit versioned schema (``.markerfinder/.pipeline_state.json``).
    Raises ``StateSchemaError`` when a state file is present but has an
    incompatible ``schema_version`` (so corrupt/foreign state fails loudly
    instead of silently being treated as "no state").
    """
    return load_pipeline_state(output_dir)


def _truncate_state_for_redo(output_dir: str, step: str) -> None:
    """Reset state so that ``step`` and any later steps are no longer marked done.

    This makes re-running a completed step safe: downstream steps will be
    recomputed from the new results of the redone step. Object files are left
    intact; only the ``completed_steps`` index is rewritten.
    """
    from markerfinder.pipeline import MarkerFinderPipeline

    step_order = MarkerFinderPipeline.STEP_ORDER
    try:
        truncate_pipeline_state(output_dir, step, step_order)
    except StateSchemaError as e:
        logger.warning(f"  Failed to update state for redo: {e}")
    except OSError as e:
        logger.warning(f"  Failed to update state for redo: {e}")


def _step_status(command: str, output_dir: str) -> dict:
    """Inspect subcommand state and return a status dictionary.

    Returns keys:
      - state_exists: bool
      - completed_steps: list
      - requested_step: str
      - requested_index: int
      - already_completed: bool (requested step already done)
      - prerequisite_ok: bool (previous step done, or command == scan)
      - missing_prerequisites: list of step names that need to run first
      - next_needed_step: the first missing prerequisite, or None
    """
    from markerfinder.pipeline import MarkerFinderPipeline

    state = _load_step_state(output_dir)
    completed = state.get("completed_steps", []) if state else []
    step_order = MarkerFinderPipeline.STEP_ORDER
    requested_idx = step_order.index(command)
    already_completed = command in completed

    missing = []
    if command != "scan":
        for step in step_order[:requested_idx]:
            if step not in completed:
                missing.append(step)

    return {
        "state_exists": state is not None,
        "completed_steps": completed,
        "requested_step": command,
        "requested_index": requested_idx,
        "already_completed": already_completed,
        "prerequisite_ok": len(missing) == 0,
        "missing_prerequisites": missing,
        "next_needed_step": missing[0] if missing else None,
    }


def _prompt_choice(message: str, choices: Dict[str, str]) -> Optional[str]:
    """Interactive prompt for TTY; returns None otherwise.

    choices maps single-letter key to description.
    """
    if not sys.stdin.isatty():
        return None
    choice_keys = "/".join(choices.keys())
    print(f"\n{message}")
    for key, desc in choices.items():
        print(f"  [{key}] {desc}")
    try:
        answer = input(f"Choice ({choice_keys}): ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        return None
    return answer if answer in choices else None


def _handle_step_state(command: str, args) -> Tuple[Optional[str], int]:
    """Inspect state for a subcommand and decide how to proceed.

    Returns a tuple of (effective_command, exit_code). ``effective_command``
    is the step to run, or None when the user chooses to abort. ``exit_code``
    is used when effective_command is None.
    """
    output_dir = args.output or "./output"
    status = _step_status(command, output_dir)

    if not status["state_exists"]:
        # No prior state: just run as requested.
        return command, EXIT_SUCCESS

    # Case 1: requested step already completed.
    if status["already_completed"]:
        if args.redo:
            return command, EXIT_SUCCESS
        completed_str = ", ".join(status["completed_steps"]) or "none"
        msg = (
            f"Step '{command}' was already completed. "
            f"Completed so far: {completed_str}."
        )
        logger.warning(msg)
        if args.resume:
            # --resume with an already-completed step means "keep going"
            # To the next step after it, if any.
            step_order = ["scan", "filter", "infer", "report"]
            idx = step_order.index(command)
            if idx + 1 < len(step_order):
                next_step = step_order[idx + 1]
                if next_step == command:
                    logger.info("  --resume: no further steps available.")
                    return None, EXIT_SUCCESS
                logger.info(f"  --resume: advancing to next step '{next_step}'")
                return next_step, EXIT_SUCCESS
            else:
                logger.info("  --resume: requested step is the final step; nothing to resume.")
                return None, EXIT_SUCCESS
        choice = _prompt_choice(
            "What would you like to do?",
            {
                "r": "redo this step (overwrite previous results)",
                "n": "start a new run in a different output directory (re-run with -o <new_dir>)",
                "s": "skip and exit",
            },
        )
        if choice == "r":
            args.redo = True
            return command, EXIT_SUCCESS
        elif choice == "n":
            logger.warning("  To start a new run, specify a different -o/--output directory.")
            return None, EXIT_SUCCESS
        else:
            logger.warning("  Exiting. Use --redo to re-run this step non-interactively.")
            return None, EXIT_SUCCESS

    # Case 2: prerequisite(s) missing.
    if not status["prerequisite_ok"]:
        if args.resume:
            return command, EXIT_SUCCESS
        missing = status["missing_prerequisites"]
        completed_str = ", ".join(status["completed_steps"]) or "none"
        next_step = status["next_needed_step"]
        logger.warning(
            f"Step '{command}' requires prior step(s) to be completed first. "
            f"Completed so far: {completed_str}. "
            f"Missing: {', '.join(missing)}."
        )
        choice = _prompt_choice(
            "What would you like to do?",
            {
                "c": f"continue by running '{next_step}' first, then proceed to '{command}'",
                "a": "abort",
            },
        )
        if choice == "c":
            args.resume = True
            return command, EXIT_SUCCESS
        else:
            logger.warning(
                f"  Exiting. Run 'markerfinder {next_step} -o {output_dir}' first, "
                "or use --resume to auto-run missing steps non-interactively."
            )
            return None, EXIT_ARG_ERROR

    # Case 3: prerequisites OK and step not yet done.
    return command, EXIT_SUCCESS
