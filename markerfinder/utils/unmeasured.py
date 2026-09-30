"""A run has to be able to say which metrics it did *not* measure, and why.

The requirement is explicit: "外部工具/依赖不可用时，日志与报告必须显式声明哪些指标
因此未被测量，不得静默降级为数值" (P0, evidence site ``hgt_filter.py:52-54``).

Without this module a run says two things and neither of them is the truth:
``evidence_coverage: 0.30`` in the summary (a fraction, which names no metric and
no cause) and per-marker ``None`` values that render as ``NA``. A reader could see
that 70% of markers were unmeasured without learning that the reason was
"ete3 unusable in this interpreter" or "no CheckM results supplied" — and the two
have completely different remedies. That gap is how a dependency quietly becoming
optional turns into a wrong conclusion.

The ledger is process-global on purpose: the call sites that discover an
unmeasurable metric are deep inside per-marker loops, and threading a collector
through six signatures would either be forgotten in the next refactor (an easy
pattern to slip back into) or push the recording back into the callers that
cannot see the cause. Determinism is preserved by sorting on render and by
clearing at the start of every run; tests must clear in a fixture.
"""

from __future__ import annotations

import logging
import threading
from typing import Dict, List, Optional, Set, Tuple

_LOCK = threading.Lock()
# (metric, cause) -> marker ids; an entry with an empty id set is a run-level gap.
_ENTRIES: Dict[Tuple[str, str], Set[str]] = {}

# How many marker ids to name per line before switching to a count.
MAX_NAMED_MARKERS = 8


def record(metric: str, cause: str, marker_id: Optional[str] = None) -> None:
    """Note that ``metric`` was not produced, because of ``cause``."""
    metric = (metric or "").strip() or "<unnamed metric>"
    cause = (cause or "").strip() or "no cause reported"
    with _LOCK:
        ids = _ENTRIES.setdefault((metric, cause), set())
        if marker_id:
            ids.add(str(marker_id))


def rollup() -> List[Tuple[str, str, int, List[str]]]:
    """Sorted ``(metric, cause, n_markers, named_marker_ids)``."""
    with _LOCK:
        snapshot = list(_ENTRIES.items())
    rows = []
    for (metric, cause), ids in snapshot:
        ordered = sorted(ids)
        rows.append((metric, cause, len(ordered), ordered[:MAX_NAMED_MARKERS]))
    return sorted(rows, key=lambda row: (row[0], row[1]))


def clear() -> None:
    with _LOCK:
        _ENTRIES.clear()


def render_lines() -> List[str]:
    """The declaration block, identical for the log and the summary product.

    Empty when nothing was recorded: silence here means "every declared metric was
    measured", which is only trustworthy because the call sites that could not
    measure are the ones writing into this ledger.
    """
    rows = rollup()
    if not rows:
        return []
    lines = [
        "METRICS NOT MEASURED THIS RUN — reported as NA/None, never as a",
        "number, and every conclusion below is computed on what WAS measured:",
    ]
    for metric, cause, n_markers, named in rows:
        if n_markers == 0:
            scope = "whole run"
        elif n_markers > len(named):
            scope = (
                f"{n_markers} marker(s): {', '.join(named)} "
                f"(+{n_markers - len(named)} more)"
            )
        else:
            scope = f"{n_markers} marker(s): {', '.join(named)}"
        lines.append(f"  - {metric}: {cause} [{scope}]")
    lines.append("")
    return lines


def declare_in_log(logger: logging.Logger, level: int = logging.WARNING) -> int:
    """Emit:func:`render_lines` through ``logger``; return how many lines went.

     names the log and the report as two separate carriers, so the same
    block has to reach both. Keeping the emission here (instead of inline in
    ``run``) makes it testable without a five-phase run with external tools.
    """
    lines = render_lines()
    for line in lines:
        logger.log(level, line)
    return len(lines)
