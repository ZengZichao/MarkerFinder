"""Shared runner for the scan/filter/infer/report pipeline subcommands.

``_execute_pipeline`` is a verbatim relocation of the single-step branch that
used to live in ``markerfinder.__main__.main``. Each subcommand module
(``scan.py`` / ``filter.py`` / ``infer.py`` / ``report.py``) exposes a thin
``execute`` wrapper that targets its step; ``main`` dispatches to the right
one. Behaviour — step-state handling, ``--resume`` / ``--redo`` logic, the
ordered ``run_step`` loop, and the exit codes — is identical to the original.
"""

from __future__ import annotations

import logging

from markerfinder.exceptions import (
    AssertionFailureError,
    PhyloToolError,
    UnsupportedCriterion,
)
from markerfinder.utils.io import load_genomes_from_directory
from markerfinder.cli.constants import (
    EXIT_SUCCESS,
    EXIT_ARG_ERROR,
    EXIT_ASSERTION_FAILED,
    EXIT_INCONCLUSIVE,
    EXIT_DATA_ERROR,
    EXIT_RUNTIME_ERROR,
)
from markerfinder.cli.step_state import (
    _handle_step_state,
    _step_status,
    _truncate_state_for_redo,
)

logger = logging.getLogger(__name__)


def _execute_pipeline(command, args, pipeline, tree, table_taxa, taxonomy_origin, logger) -> int:
    """Run one or more ordered pipeline steps for a subcommand.

    Behaviour is identical to the historical single-step branch of ``main``:
    resolve step state, honour ``--resume`` / ``--redo``, then execute each
    step via ``pipeline.run_step``.
    """
    effective_command, abort_code = _handle_step_state(command, args)
    if effective_command is None:
        return abort_code

    # Determine the ordered list of steps to execute. When --resume is set,
    # This includes any missing prerequisites plus the requested step (or the
    # Next step after an already-completed one).
    step_order = ["scan", "filter", "infer", "report"]
    if args.resume:
        # Recalculate dynamically; _handle_step_state may have redirected us
        # To a different step.
        target_idx = step_order.index(effective_command)
        status = _step_status(
            effective_command, args.output or "./output",
        )
        # Run from the first missing prerequisite up to the effective command.
        start_idx = step_order.index(status["next_needed_step"]) if status["next_needed_step"] else target_idx
        steps_to_run = step_order[start_idx:target_idx + 1]
    else:
        steps_to_run = [effective_command]

    # If the user explicitly asked to redo the originally requested step,
    # Remove that step and any later steps from the persisted state so that
    # Downstream steps are recomputed from the new results.
    if args.redo and effective_command == command:
        _truncate_state_for_redo(args.output or "./output", command)

    # --resume 会把缺失的前置步骤(通常是 scan)展开进 steps_to_run;
    # 非 scan 子命令不强制 -i, 此时必须在启动前给出可操作的报错, 而不是让
    # Load_genomes_from_directory(None) 以晦涩的 TypeError 崩掉.
    if "scan" in steps_to_run and not getattr(args, "input", None):
        logger.error(
            "--resume needs to run the 'scan' step first, which requires "
            "-i/--input. Re-run with: markerfinder "
            f"{command} -i <genome_dir> -o {args.output or './output'} --resume"
        )
        return EXIT_ARG_ERROR

    for step in steps_to_run:
        genomes = None
        reference_tree = None
        table_taxa_step = None
        if step == "scan":
            genomes = load_genomes_from_directory(args.input)
            if not genomes:
                logger.error(f"No genomes found in {args.input}")
                return EXIT_DATA_ERROR
            reference_tree = tree
            table_taxa_step = table_taxa

        try:
            step_result = pipeline.run_step(
                step,
                genomes=genomes,
                reference_tree=reference_tree,
                table_taxa=table_taxa_step,
                taxonomy_origin=taxonomy_origin,
            )
        except RuntimeError as e:
            # Missing prerequisite step (e.g. 'filter' before 'scan')
            logger.error(str(e))
            return EXIT_ARG_ERROR
        except UnsupportedCriterion as e:
            # Gated criterion refused to start.
            logger.error(str(e))
            return EXIT_ARG_ERROR
        except AssertionFailureError as e:
            # Output assertions failed — distinct exit code 4.
            logger.error(f"Assertion failure: {e}")
            return EXIT_ASSERTION_FAILED
        except PhyloToolError as e:
            logger.error(f"Pipeline error: {e}")
            return EXIT_DATA_ERROR
        except Exception as e:
            logger.error(f"Step '{step}' failed: {e}")
            return EXIT_RUNTIME_ERROR

        # A finished 'infer' step whose recommendation layer could
        # Not name a resolved tree answers with exit 5 here too. Without this the
        # One-shot run and the stepwise route report the same scientific state with
        # Different codes, and a scheduler driving the steps would treat an
        # Inconclusive result as a success.
        if step == "infer":
            from markerfinder.cli.constants import recommendation_is_inconclusive

            phylo = (step_result or {}).get("phylo_result")
            if recommendation_is_inconclusive(
                getattr(phylo, "tree_recommendation", None)
            ):
                logger.error(
                    "Step 'infer' completed, but the recommendation layer is "
                    "INCONCLUSIVE (no resolved tree): exit code 5, not a success."
                )
                return EXIT_INCONCLUSIVE

    return EXIT_SUCCESS
