"""Explicit, versioned schema serialization for MarkerFinder subcommand state.

This module replaces the legacy whole-object ``pickle`` persistence used by
``MarkerFinderPipeline.run_step`` and the ``__main__`` step-state helpers.

Design goals
------------------------------------
* **Single source of truth**: one on-disk layout, shared by ``pipeline.py`` and
  ``__main__.py``.
* **Explicit schema**: a ``PipelineState`` index (dataclass, ``schema_version``)
  records only JSON-safe scalars/paths. Heavy objects are *never* pickled as a
  blob.
* **Trees are never JSON-ified**: ``ete3``/model ``Tree`` objects are written to
  NEWICK files and referenced by relative path; the codec refuses to embed them
  in JSON.
* **schema_version validation**: a present-but-incompatible state file raises a
  clear ``StateSchemaError`` instead of silently loading garbage / dirty state.

On-disk layout (under ``<output_dir>/.markerfinder``)
----------------------------------------------------
* ``.pipeline_state.json`` — the explicit ``PipelineState`` index.
* ``objects/<key>.json`` — one typed-JSON encoded object per logical key.
* ``objects/tree_<n>.nwk`` — NEWICK files for persisted trees (referenced by
  relative path from the index).

The typed-JSON codec records the concrete runtime type of every value (so the
many ``Optional[object]`` fields on the data models reconstruct unambiguously)
and special-cases ``Tree`` to a NEWICK file reference.
"""

from __future__ import annotations

import logging

import base64
import importlib
import json
import threading
from dataclasses import dataclass, field, fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

from markerfinder._version import get_version_string

__all__ = [
    "SCHEMA_VERSION",
    "STATE_DIR_NAME",
    "STATE_FILE_NAME",
    "PipelineState",
    "StateSchemaError",
    "state_dir_for",
    "save_pipeline_state",
    "load_pipeline_state",
    "read_completed_steps",
    "truncate_pipeline_state",
    "prev_step_product_rel",
]

logger = logging.getLogger(__name__)

# V2 adds per-marker decision-card provenance fields.
# Old state MUST NOT be silently loaded; _assert_compatible
# Raises with an actionable re-run hint.
SCHEMA_VERSION = 2
STATE_DIR_NAME = ".markerfinder"
STATE_FILE_NAME = ".pipeline_state.json"
OBJECTS_DIR_NAME = "objects"

# Scalar fields stored directly on the PipelineState index (not as object files).
_SCALAR_KEYS = {"completed_steps", "start_time", "taxonomy_map"}

# Maps a step to the primary product object key produced by the *previous* step.
# Used by the atomicity guard in ``run_step`` so a corrupt/missing prior artifact
# Fails loudly instead of producing dirty downstream state.
_PREV_STEP_PRODUCT: Dict[str, str] = {
    "filter": "marker_sequences",
    "infer": "hgt_report",
    "report": "phylo_result",
}

_TAG = "__mf__"


class StateSchemaError(Exception):
    """Raised when persisted state cannot be loaded (bad schema / version)."""


@dataclass
class PipelineState:
    schema_version: int = SCHEMA_VERSION
    completed_steps: List[str] = field(default_factory=list)
    start_time: float = 0.0
    config_version: str = ""
    taxonomy_map: Optional[dict] = None
    objects: Dict[str, str] = field(default_factory=dict)

    def to_payload(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "completed_steps": self.completed_steps,
            "start_time": self.start_time,
            "config_version": self.config_version,
            "taxonomy_map": self.taxonomy_map,
            "objects": self.objects,
        }


class _CodecContext:
    """Per-save state for the recursive codec (tree file naming, imports)."""

    def __init__(self, state_dir: Path):
        self.state_dir = state_dir
        self._lock = threading.Lock()
        self._tree_counter = 0

    def next_tree_path(self) -> Path:
        with self._lock:
            self._tree_counter += 1
            return self.state_dir / OBJECTS_DIR_NAME / f"tree_{self._tree_counter}.nwk"


_CLASS_CACHE: Dict[str, Any] = {}


def _import_fqn(fqn: str) -> Any:
    cls = _CLASS_CACHE.get(fqn)
    if cls is not None:
        return cls
    module_name, attr = fqn.rsplit(".", 1)
    module = importlib.import_module(module_name)
    cls = getattr(module, attr)
    _CLASS_CACHE[fqn] = cls
    return cls


def _fqn(obj: Any) -> str:
    return f"{type(obj).__module__}.{type(obj).__qualname__}"


def _encode(value: Any, ctx: _CodecContext) -> Any:
    """Recursively encode ``value`` into a JSON-safe, type-tagged node."""
    if value is None:
        return {_TAG: True, "t": "none", "d": None}
    if isinstance(value, bool):
        return {_TAG: True, "t": "bool", "d": value}
    if isinstance(value, int):
        return {_TAG: True, "t": "int", "d": value}
    if isinstance(value, float):
        return {_TAG: True, "t": "float", "d": value}
    if isinstance(value, str):
        return {_TAG: True, "t": "str", "d": value}
    if isinstance(value, bytes):
        return {_TAG: True, "t": "bytes", "d": base64.b64encode(value).decode("ascii")}

    # Tree (ete3-backed model) -> NEWICK file reference (never embedded in JSON).
    from markerfinder.models.tree import Tree

    if isinstance(value, Tree):
        tree_path = ctx.next_tree_path()
        tree_path.parent.mkdir(parents=True, exist_ok=True)
        tree_path.write_text(value.newick, encoding="utf-8", newline="\n")
        rel = str(tree_path.relative_to(ctx.state_dir))
        return {_TAG: True, "t": "tree", "d": rel}

    if isinstance(value, Enum):
        return {
            _TAG: True,
            "t": "enum",
            "d": {"cls": _fqn(value), "name": value.name, "value": value.value},
        }
    if isinstance(value, tuple):
        return {_TAG: True, "t": "tuple", "d": [_encode(v, ctx) for v in value]}
    if isinstance(value, set):
        return {_TAG: True, "t": "set", "d": [_encode(v, ctx) for v in value]}
    if isinstance(value, list):
        return {_TAG: True, "t": "list", "d": [_encode(v, ctx) for v in value]}
    if isinstance(value, dict):
        return {
            _TAG: True,
            "t": "dict",
            "d": [[_encode(k, ctx), _encode(v, ctx)] for k, v in value.items()],
        }
    if is_dataclass(value) and not isinstance(value, type):
        encoded_fields = {
            f.name: _encode(getattr(value, f.name), ctx) for f in fields(value)
        }
        return {
            _TAG: True,
            "t": "dataclass",
            "d": {"cls": _fqn(value), "fields": encoded_fields},
        }

    # Last-resort: anything natively JSON-serializable. Should not be reached for
    # The known MarkerFinder object graph.
    try:
        json.dumps(value)
    except (TypeError, ValueError):
        raise StateSchemaError(
            f"Cannot encode value of type {type(value).__name__!r} into the "
            f"explicit pipeline-state schema."
        )
    return {_TAG: True, "t": "json", "d": value}


def _decode(node: Any, ctx: _CodecContext) -> Any:
    """Recursively decode a type-tagged node back into a Python object."""
    if not (isinstance(node, dict) and node.get(_TAG) is True):
        # Be lenient with untagged values (defensive; our own files are tagged).
        return node
    t = node["t"]
    d = node["d"]
    if t == "none":
        return None
    if t == "bool":
        return bool(d)
    if t == "int":
        return int(d)
    if t == "float":
        return float(d)
    if t == "str":
        return str(d)
    if t == "bytes":
        return base64.b64decode(d)
    if t == "tree":
        from markerfinder.models.tree import Tree

        return Tree.read(str(ctx.state_dir / d))
    if t == "enum":
        cls = _import_fqn(d["cls"])
        try:
            return cls[d["name"]]
        except (KeyError, ValueError):
            return cls(d["value"])
    if t == "tuple":
        return tuple(_decode(x, ctx) for x in d)
    if t == "set":
        return set(_decode(x, ctx) for x in d)
    if t == "list":
        return [_decode(x, ctx) for x in d]
    if t == "dict":
        return {_decode(k, ctx): _decode(v, ctx) for k, v in d}
    if t == "dataclass":
        cls = _import_fqn(d["cls"])
        kwargs = {fname: _decode(v, ctx) for fname, v in d["fields"].items()}
        return cls(**kwargs)
    if t == "json":
        return d
    raise StateSchemaError(f"Unknown pipeline-state schema type tag: {t!r}")


def state_dir_for(output_dir: str) -> Path:
    """Return the ``.markerfinder`` state directory for ``output_dir``."""
    return Path(output_dir) / STATE_DIR_NAME


def save_pipeline_state(output_dir: str, state: Dict[str, Any]) -> None:
    """Persist a ``run_step`` ``state`` dict using the explicit schema.

    Scalar fields (``completed_steps``, ``start_time``, ``taxonomy_map``) go into
    the ``PipelineState`` index; every other key is encoded to its own object file.
    Trees inside those objects are written as NEWICK files and referenced by path.
    """
    sdir = state_dir_for(output_dir)
    obj_dir = sdir / OBJECTS_DIR_NAME
    obj_dir.mkdir(parents=True, exist_ok=True)
    ctx = _CodecContext(sdir)

    completed_steps = list(state.get("completed_steps", []))
    start_time = float(state.get("start_time", 0.0) or 0.0)
    taxonomy_map = state.get("taxonomy_map")
    if taxonomy_map is not None and not isinstance(taxonomy_map, dict):
        taxonomy_map = dict(taxonomy_map)

    objects: Dict[str, str] = {}
    for key, value in state.items():
        if key in _SCALAR_KEYS:
            continue
        rel = f"{OBJECTS_DIR_NAME}/{key}.json"
        payload = _encode(value, ctx)
        (sdir / rel).write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
            newline="\n",
        )
        objects[key] = rel

    ps = PipelineState(
        schema_version=SCHEMA_VERSION,
        completed_steps=completed_steps,
        start_time=start_time,
        config_version=get_version_string(),
        taxonomy_map=taxonomy_map,
        objects=objects,
    )
    (sdir / STATE_FILE_NAME).write_text(
        json.dumps(ps.to_payload(), ensure_ascii=False, indent=2),
        encoding="utf-8",
        newline="\n",
    )


def load_pipeline_state(output_dir: str) -> Optional[Dict[str, Any]]:
    """Load a persisted ``run_step`` state dict, or ``None`` if absent.

    Raises ``StateSchemaError`` when a state file is present but unreadable or has
    an incompatible ``schema_version`` (so callers fail loudly, never on dirty
    state).
    """
    sdir = state_dir_for(output_dir)
    state_file = sdir / STATE_FILE_NAME
    if not state_file.exists():
        return None
    try:
        raw = json.loads(state_file.read_text(encoding="utf-8"))
    except Exception as e:  # Noqa: BLE001 - surface any read/parse failure clearly
        raise StateSchemaError(f"Failed to read pipeline state: {e}")
    _assert_compatible(raw)
    ctx = _CodecContext(sdir)
    state: Dict[str, Any] = {
        "completed_steps": list(raw.get("completed_steps", [])),
        "start_time": float(raw.get("start_time", 0.0) or 0.0),
        "taxonomy_map": raw.get("taxonomy_map"),
    }
    for key, rel in raw.get("objects", {}).items():
        node = json.loads((sdir / rel).read_text(encoding="utf-8"))
        state[key] = _decode(node, ctx)
    return state


def read_completed_steps(output_dir: str) -> Optional[List[str]]:
    """Return completed steps, or ``None`` if no (compatible) state exists."""
    sdir = state_dir_for(output_dir)
    state_file = sdir / STATE_FILE_NAME
    if not state_file.exists():
        return None
    try:
        raw = json.loads(state_file.read_text(encoding="utf-8"))
    except Exception as e:  # Noqa: BLE001 - unreadable state is logged, not hidden
        logger.warning(
            f"Pipeline state at {state_file} is unreadable ({e}); "
            "treating as no state. Re-run with --redo if this persists."
        )
        return None
    if raw.get("schema_version") != SCHEMA_VERSION:
        # A version mismatch is NOT "no state" — silently
        # Dropping it would discard --resume context without a trace.
        raise StateSchemaError(
            f"Pipeline state at {state_file} was written by schema_version="
            f"{raw.get('schema_version')!r} but this MarkerFinder expects "
            f"{SCHEMA_VERSION} (v2 adds decision-card provenance). The state is "
            f"NOT loadable; re-run the affected steps with --redo (steps are "
            f"recomputed from scratch, no data is lost)."
        )
    return list(raw.get("completed_steps", []))


def truncate_pipeline_state(output_dir: str, step: str, step_order: List[str]) -> None:
    """Drop ``step`` and later steps from ``completed_steps`` (for --redo).

    Only the index is rewritten; object files are left intact so a redone step
    recomputes and overwrites them. Raises ``StateSchemaError`` on incompatible
    state.
    """
    sdir = state_dir_for(output_dir)
    state_file = sdir / STATE_FILE_NAME
    if not state_file.exists():
        return
    try:
        raw = json.loads(state_file.read_text(encoding="utf-8"))
    except Exception as e:  # Noqa: BLE001
        raise StateSchemaError(f"Failed to read pipeline state for redo: {e}")
    _assert_compatible(raw)
    if step not in step_order:
        return
    idx = step_order.index(step)
    raw["completed_steps"] = [
        s for s in raw.get("completed_steps", []) if step_order.index(s) < idx
    ]
    state_file.write_text(
        json.dumps(raw, ensure_ascii=False, indent=2),
        encoding="utf-8",
        newline="\n",
    )


def prev_step_product_rel(step: str) -> Optional[str]:
    """Relative path (from the state dir) of the previous step's primary product.

    Returns ``None`` for ``scan`` (no predecessor).
    """
    key = _PREV_STEP_PRODUCT.get(step)
    if key is None:
        return None
    return f"{OBJECTS_DIR_NAME}/{key}.json"


def _assert_compatible(raw: dict) -> None:
    version = raw.get("schema_version")
    if version != SCHEMA_VERSION:
        raise StateSchemaError(
            f"This state was written by schema_version={version!r} but this "
            f"MarkerFinder expects {SCHEMA_VERSION} (v2 adds decision-card "
            f"provenance fields). The state cannot be resumed; re-run with "
            f"--redo to recompute the affected steps (Phase 2 onwards if you "
            f"have scan output), or delete {STATE_DIR_NAME}/ to start fresh."
        )
