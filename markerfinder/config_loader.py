"""Configuration file support: YAML, TOML, and JSON loading.

Provides --config <path> with CLI > Config > Defaults priority.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

from markerfinder.exceptions import ConfigError

logger = logging.getLogger(__name__)


def _load_yaml(path: str) -> Dict[str, Any]:
    """Load a YAML configuration file."""
    try:
        import yaml
    except ImportError:
        raise ConfigError(
            "PyYAML is required for YAML config files. "
            "Install with: pip install pyyaml"
        )

    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"YAML config must be a mapping, got {type(data).__name__}")
    return data


def _load_toml(path: str) -> Dict[str, Any]:
    """Load a TOML configuration file."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            raise ConfigError(
                "tomllib (Python 3.11+) or tomli is required for TOML config files. "
                "Install with: pip install tomli"
            )

    with open(path, "rb") as f:
        data = tomllib.load(f)

    if not isinstance(data, dict):
        raise ConfigError(f"TOML config must be a mapping, got {type(data).__name__}")
    return data


def _load_json(path: str) -> Dict[str, Any]:
    """Load a JSON configuration file (e.g. run_config.json)."""
    import json
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ConfigError(f"JSON config must be a mapping, got {type(data).__name__}")
    return data


def load_config_file(path: str) -> Dict[str, Any]:
    """Load a configuration file (YAML, TOML, or JSON).

    A top-level ``parameters`` mapping (e.g. the ``run_config.json`` written by
    a previous run) is flattened automatically so that
    ``markerfinder --config previous_run/Phase5_metadata/run_config.json``
    actually replays the recorded parameters instead of silently ignoring them
    (README「可复现性」承诺).

    Args:
        path: path to the configuration file.

    Returns:
        Dictionary of configuration values.

    Raises:
        ConfigError: if file does not exist, cannot be parsed, or has wrong format.
    """
    config_path = Path(path)

    if not config_path.exists():
        raise ConfigError(f"Configuration file not found: {path}")

    suffix = config_path.suffix.lower()

    if suffix in (".yaml", ".yml"):
        data = _load_yaml(path)
    elif suffix == ".toml":
        data = _load_toml(path)
    elif suffix == ".json":
        data = _load_json(path)
    else:
        try:
            data = _load_yaml(path)
        except ConfigError:
            data = None
        if data is None:
            try:
                data = _load_toml(path)
            except ConfigError:
                data = None
        if data is None:
            try:
                data = _load_json(path)
            except ConfigError:
                data = None
        if data is None:
            raise ConfigError(
                f"Cannot determine config format for '{path}'. "
                f"Use .yaml, .yml, .toml, or .json extension."
            )

    data = _flatten_parameters(data)
    return _flatten_snapshot(data) if _is_snapshot(data) else data


# Path-valued fields a run_config snapshot records. The writer resolves them to
# Absolute paths (see ``Pipeline._record_paths_resolved``), because a relative
# Path in a snapshot cannot be interpreted later without knowing which directory
# It was made relative to: measured on a replay of a snapshot whose ``db_dir``
# Was recorded as ``./db``, the value silently pointed at the previous run's
# Output directory.
SNAPSHOT_PATH_KEYS = frozenset({
    "input_dir", "output_dir", "tmp_dir", "db_dir", "template_dir",
    "checkm_results", "taxonomy_table",
    # The same fields under the names the command line uses.
    "input", "output", "marker_hmm_dir", "gtdb_markers_dir", "species_tree",
    "sequences", "cog_category_map",
})
SNAPSHOT_DIR_NAME = "Phase5_metadata"
SNAPSHOT_MARKERS = frozenset({
    "parameters", "database_versions", "run_duration_seconds",
    "markerfinder_version",
})

# The recording fields of a snapshot that are facts about the run rather than
# Parameters of it; they have no command-line equivalent and must not be
# Reported as "unrecognized keys" every time a snapshot is replayed. The names
# Of the envelope itself (``parameters``, ``database_versions``, …) are already
# In ``SNAPSHOT_MARKERS``; structural leftovers such as ``databases`` are left
# Alone so the CLI keeps reporting what it cannot express.
SNAPSHOT_ENVELOPE_KEYS = frozenset({
    "timestamp", "run_duration_seconds", "markerfinder_version", "tmp_dir_auto",
}) | SNAPSHOT_MARKERS

# Scratch, not configuration: the recorded temporary directory belongs to the
# Run that wrote the snapshot. Replaying it would reuse another run's working
# Files (and report ``tmp_dir_auto=false``, which changes whether they are
# Cleaned up), so a replay allocates its own.
SNAPSHOT_SCRATCH_KEYS = frozenset({"tmp_dir"})

# Fields the snapshot records under a name the command line does not use.
# ``PipelineConfig.to_dict`` is written from the dataclasses, while ``--config``
# Feeds the argparse namespace, so the two vocabularies differ in exactly these
# Places; without the mapping a replay silently loses the marker directory, the
# Taxonomy table and the risk bands, and reports them as unrecognised.
SNAPSHOT_KEY_ALIASES = {
    "input_dir": "input",
    "output_dir": "output",
    "cpus": "threads",
    "sequences_path": "sequences",
    "level1_max": "hgt_threshold_l1l2",
    "level2_max": "hgt_threshold_l2l3",
    "adaptive_far_thresholds": "hgt_adaptive_thresholds",
    "ufboot_replicates": "ufboot",
}

# One more level of nesting that carries real parameters.
SNAPSHOT_NESTED_BLOCKS = frozenset({"level_thresholds"})

# Where the snapshot keeps the *requested* value next to the value Phase 0
# Adapted it into. A replay must reproduce what the user asked for; feeding back
# The adapted value would let a budget of 4 replay as a budget of 150.
SNAPSHOT_INTENT_FIELDS = {
    "max_markers": "user_max_markers",
    "min_occupancy": "user_min_occupancy",
}

# Other fields that mean the same thing as an intent field. The MAG module keeps
# Its own copy of the occupancy floor under the deprecated name, and its recorded
# Default would otherwise arrive as a *request* on replay — measured as a replay
# Whose occupancy floor was 0.3 where the recorded run had applied 0.75.
SNAPSHOT_INTENT_ALIASES = {
    "min_occupancy": ("min_marker_coverage",),
    "max_markers": (),
}


def _flatten_snapshot(data: Dict[str, Any]) -> Dict[str, Any]:
    """Turn a recorded ``run_config.json`` back into flat configuration.

    The snapshot nests the per-module dataclasses (``selection_config``,
    ``hgt_config``, ``taxonomy_config``, …) and renames a handful of fields, so
    reading it verbatim leaves the values that actually drive a run —
    ``gtdb_markers_dir``, ``taxonomy_table``, the risk bands — buried inside
    blocks the command line has no name for. Those blocks are flattened here and
    the renamed keys restored, so ``--config <previous run>`` reproduces the
    previous run instead of falling back to defaults while claiming otherwise.
    """
    flat: Dict[str, Any] = {}

    def put(key: str, value: Any) -> None:
        key = SNAPSHOT_KEY_ALIASES.get(key, key)
        flat.setdefault(key, value)

    for key, value in data.items():
        if key == "parameters" and isinstance(value, dict):
            # A snapshot whose 'parameters' block was already lifted out by
            # _flatten_parameters still carries the sub-blocks at top level.
            for inner, inner_value in value.items():
                put(inner, inner_value)
            continue
        if key.endswith("_config") and isinstance(value, dict):
            for inner, inner_value in value.items():
                if inner in SNAPSHOT_NESTED_BLOCKS and isinstance(inner_value, dict):
                    for deep_key, deep_value in inner_value.items():
                        put(deep_key, deep_value)
                    continue
                put(inner, inner_value)
            continue
        put(key, value)

    for adapted, requested in SNAPSHOT_INTENT_FIELDS.items():
        if requested not in flat:
            continue
        intent = flat.pop(requested)
        if intent is not None:
            # What the user asked for is the reproducible fact.
            flat[adapted] = intent
        else:
            # Nothing was asked for: the value recorded under ``adapted`` is what
            # Phase 0 derived for that dataset. Feeding a derived value back as
            # If it had been requested turns an adaptive default into a pinned
            # One, so the field (and any module's copy of it under the
            # Deprecated alias name) is dropped and re-derived by the replay.
            flat.pop(adapted, None)
            for alias in SNAPSHOT_INTENT_ALIASES.get(adapted, ()):
                flat.pop(alias, None)

    dropped = SNAPSHOT_ENVELOPE_KEYS | SNAPSHOT_SCRATCH_KEYS
    if SNAPSHOT_SCRATCH_KEYS & set(flat):
        logger.debug(
            "Replaying a recorded run_config.json: %d field(s) recovered from "
            "the nested parameter blocks; %s dropped as per-run scratch",
            len(flat), ",".join(sorted(SNAPSHOT_SCRATCH_KEYS)),
        )
    return {k: v for k, v in flat.items() if k not in dropped}


def _is_snapshot(data: Any) -> bool:
    return isinstance(data, dict) and bool(SNAPSHOT_MARKERS & set(data))


def _flatten_parameters(data: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten a nested ``parameters`` mapping into the top level.

    ``Phase5_metadata/run_config.json`` stores the parameter snapshot under a
    ``parameters`` key. Without flattening, ``merge_config_with_args`` only
    sees unrelated top-level keys (version/timestamp/...) and the replay is a
    silent no-op. Non-mapping ``parameters`` values are left untouched.
    """
    if not isinstance(data, dict):
        return data
    params = data.get("parameters")
    if not isinstance(params, dict):
        return data
    flattened = {k: v for k, v in data.items() if k != "parameters"}
    flattened.update(params)
    logger.debug(
        "Config file: flattened %d entr(ies) from nested 'parameters' mapping",
        len(params),
    )
    return flattened


def merge_config_with_args(
    config_data: Dict[str, Any],
    args_dict: Dict[str, Any],
    defaults: Dict[str, Any],
) -> Dict[str, Any]:
    """Merge configuration sources with priority: CLI > Config > Defaults.

    A CLI argument is considered "explicitly set" only if its value differs
    from the argparse default. This prevents non-None defaults (e.g.
    ``--threads 1``) from silently overriding config-file values.

    Args:
        config_data: values from config file.
        args_dict: values from CLI arguments (should contain only explicitly
            set values; caller should filter out defaults).
        defaults: default values.

    Returns:
        Merged dictionary.
    """
    result = dict(defaults)

    for key, value in config_data.items():
        cli_key = key.replace("-", "_")
        if cli_key not in result or value is not None:
            result[cli_key] = value

    for key, value in args_dict.items():
        if value is not None and value != defaults.get(key):
            result[key] = value

    return result
