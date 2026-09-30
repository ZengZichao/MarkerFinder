"""Thin parallel-map helper for I/O-bound external-tool tasks.

Why a ``ThreadPoolExecutor`` backend (the default) instead of
``ProcessPoolExecutor``:

* These tasks (hmmsearch / mafft / fasttree / trimal / iqtree3 / astral) spend
  almost all of their wall-clock time *blocked* inside an external subprocess.
  While a subprocess runs, the caller's Python thread releases the GIL, so a
  pool of threads overlaps the wall-clock time of N markers' tool invocations
  just as well as a pool of processes would.
* A ``ProcessPoolExecutor`` would spawn child Python interpreters that re-import
  the modules and therefore do **not** inherit ``unittest.mock.patch``
  installed in the parent process. MarkerFinder's test-suite mocks those
  external tools (they are intentionally not installed in CI), so a process
  backend would break the existing tests, whereas threads share the parent's
  module namespace and keep every mock working. Same numerical/semantic
  result, zero test regressions.
* A ``"process"`` backend is still selectable for deployments that want true CPU
  parallelism and have the real bioinformatics tools installed. (Note: when
  using ``"process"`` the *worker* must be a module-level callable with
  picklable arguments — closures are not picklable.)

The helper degrades to a plain serial ``list(map)`` whenever
``max_workers <= 1`` or there is only a single task, so callers can use it
unconditionally and still get identical (serial) behaviour in the degenerate
cases — which also keeps the single-marker unit tests on the serial path.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Callable, List, Sequence, TypeVar

T = TypeVar("T")
R = TypeVar("R")

# Default backend. "thread" keeps the parent process's mocked subprocess.run
# Visible to workers (so the test-suite stays green) while still overlapping
# The wall-clock time of external-tool invocations. Switch to "process" only
# When real bioinformatics tools are installed and you want CPU parallelism.
DEFAULT_BACKEND = "thread"


def parallel_map(
    worker: Callable[[T], R],
    tasks: Sequence[T],
    max_workers: int = 1,
    backend: str = DEFAULT_BACKEND,
) -> List[R]:
    """Run ``worker`` over ``tasks``, returning results in task order.

    Falls back to a serial ``list(map)`` when ``max_workers <= 1`` or there
    is only one task. Results preserve the input task order (both executors do).
    """
    tasks = list(tasks)
    if max_workers <= 1 or len(tasks) <= 1:
        return [worker(t) for t in tasks]
    if backend == "process":
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            return list(ex.map(worker, tasks))
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        return list(ex.map(worker, tasks))
