"""Taxonomy parsing, external table loading, and monophyly checking.

Supports two formats:
  Format A (embedded): GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_...
    Delimiters: _d_ _p_ _c_ _o_ _f_ _g_ (single letter between underscores)
  Format B (table-style): d__Bacteria;p__Cyanobacteriota;c__...
    Pairs: {level}__{value} separated by; or custom separator

Provides public API:
  parse_taxonomy(label, mode="reverse") -> dict
  is_monophyletic(tree, taxon_label, rooted=True) -> bool
  load_tree(path, validate=True) -> Tree
"""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from markerfinder.exceptions import (
    MonophylyError,
    TaxonomyConflictError,
)
from markerfinder.validation import check_control_characters

logger = logging.getLogger(__name__)

# Track which missing-level DEBUG messages have already been emitted so the
# Same message is not repeated for every label parsed in a run.
_MISSING_LEVELS_LOGGED: Set[str] = set()

LEVEL_ORDER_A = ["g", "f", "o", "c", "p", "d"]
LEVEL_MAP_A = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
LEVEL_MAP_B = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
STANDARD_LEVELS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# ``--monophyly-rank`` vocabulary.
#
# ``AUTO_RANK`` is the default: the rank the monophyly proportion is measured at
# Is derived from the tree (one rank below the taxonomic scope of its tips).
# Naming a rank instead makes that rank authoritative. ``AUTO_NO_SCOPE_RANK`` is
# What ``auto`` measures when the tips are not cohesive at any rank (e.g. they
# Span several domains), so the option still resolves to a single number.
AUTO_RANK = "auto"
AUTO_NO_SCOPE_RANK = "genus"


def _empty_taxonomy() -> Dict[str, Optional[str]]:
    return {level: None for level in STANDARD_LEVELS}


def parse_custom_levels(levels_str: str) -> Tuple[Dict[str, str], List[str]]:
    """解析自定义分类级别字符串。

    Args:
        levels_str: 格式为 'level:prefix' 的逗号分隔列表，如 'kingdom:_k_'
                     prefix 可以是 Format A 格式（_k_）或 Format B 格式（k）

    Returns:
        (custom_level_map, custom_level_names)
        custom_level_map: {prefix_code: level_name}（统一使用单字母代码）
        custom_level_names: [level_name,...]
    """
    custom_map: Dict[str, str] = {}
    custom_names: List[str] = []

    for pair in levels_str.split(","):
        pair = pair.strip()
        if not pair or ":" not in pair:
            continue
        name, _, prefix = pair.partition(":")
        name = name.strip()
        prefix = prefix.strip()
        if not name or not prefix:
            continue
        # 统一提取单字母代码：_k_ → k, k → k
        code = prefix.strip("_")
        if not code:
            continue
        custom_map[code] = name
        custom_names.append(name)
        logger.debug(f"Registered custom taxonomy level: '{name}' with code '{code}'")

    return custom_map, custom_names


def get_merged_level_maps(
    custom_levels: Optional[str] = None,
) -> Tuple[Dict[str, str], Dict[str, str], List[str]]:
    """获取合并了自定义级别的级别映射表。

    Args:
        custom_levels: 自定义级别字符串（如 'kingdom:_k_'）

    Returns:
        (merged_map_a, merged_map_b, merged_level_names)
        map_a 和 map_b 均使用单字母代码作为键（如 'k'）
    """
    map_a = dict(LEVEL_MAP_A)
    map_b = dict(LEVEL_MAP_B)
    level_names = list(STANDARD_LEVELS)

    if custom_levels:
        custom_map, custom_names = parse_custom_levels(custom_levels)
        # Custom_map 的键已经是单字母代码，可直接用于两个 map
        map_a.update(custom_map)
        map_b.update(custom_map)
        level_names.extend(custom_names)

    return map_a, map_b, level_names


def parse_taxonomy_embedded(
    label: str,
    mode: str = "reverse",
    custom_levels: Optional[str] = None,
) -> Dict[str, Optional[str]]:
    """Parse Format A (embedded taxonomy in label).

    Args:
        label: e.g. 'GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_...'
        mode: 'reverse' (right-to-left), 'greedy' (first _d_ to end), 'segment'
        custom_levels: custom level prefixes (e.g. 'kingdom:_k_')

    Returns:
        Standard taxonomy dict with keys domain..species (plus custom levels).
    """
    level_map_a, _, level_names = get_merged_level_maps(custom_levels)
    result: Dict[str, Optional[str]] = {name: None for name in level_names}

    pattern = re.compile(r"_([a-z])_")
    all_matches = list(pattern.finditer(label))

    if not all_matches:
        return result

    if mode == "reverse":
        seen_levels: List[str] = []
        positions: List[Tuple[str, int]] = []

        for m in reversed(all_matches):
            level_code = m.group(1)
            if level_code not in level_map_a:
                continue
            if level_code in seen_levels:
                logger.warning(
                    f"Duplicate level _{level_code}_ in label '{label}', "
                    f"falling back to greedy mode"
                )
                return parse_taxonomy_embedded(label, mode="greedy", custom_levels=custom_levels)
            seen_levels.append(level_code)
            positions.append((level_code, m.start()))

        positions.reverse()

        for idx, (level_code, start) in enumerate(positions):
            prefix_end = start + 3
            if idx + 1 < len(positions):
                next_start = positions[idx + 1][1]
                value = label[prefix_end:next_start]
            else:
                value = label[prefix_end:]

            value = value.strip("_")
            if not value:
                logger.debug(f"Empty value for level _{level_code}_ in '{label}'")
                continue

            level_name = level_map_a[level_code]
            result[level_name] = value

    elif mode == "greedy":
        d_match = re.search(r"_d_", label)
        if d_match is None:
            return result

        remainder = label[d_match.start():]
        parts = re.split(r"_([a-z])_", remainder)

        current_level = None
        for i, part in enumerate(parts):
            if i == 0:
                continue
            if i % 2 == 1:
                if part in level_map_a:
                    current_level = part
                else:
                    current_level = None
            else:
                if current_level and part:
                    result[level_map_a[current_level]] = part
                    current_level = None

    elif mode == "segment":
        segments = re.split(r"_(?:[a-z])_", label)
        level_codes = re.findall(r"_([a-z])_", label)

        prefix = segments[0] if segments else ""
        for idx, code in enumerate(level_codes):
            if code in level_map_a and idx + 1 < len(segments):
                value = segments[idx + 1]
                if value:
                    result[level_map_a[code]] = value
                    if value != segments[idx + 1].strip("_"):
                        logger.warning(
                            f"Ambiguous segment for _{code}_ in '{label}'"
                        )
    else:
        raise ValueError(f"Unknown taxonomy delimiter mode: {mode}")

    for level_name in level_names:
        if result[level_name] is None and level_name not in _MISSING_LEVELS_LOGGED:
            _MISSING_LEVELS_LOGGED.add(level_name)
            logger.debug(f"Missing level '{level_name}' in embedded label '{label}'")

    return result


def parse_taxonomy_table(
    taxonomy_string: str,
    sep: str = ";",
    custom_levels: Optional[str] = None,
) -> Dict[str, Optional[str]]:
    """Parse Format B (table-style taxonomy string).

    Args:
        taxonomy_string: e.g. 'd__Bacteria;p__Cyanobacteriota;c__...;s__'
        sep: separator between level-value pairs (default ';')
        custom_levels: custom level prefixes (e.g. 'kingdom:k')

    Returns:
        Standard taxonomy dict.
    """
    _, level_map_b, level_names = get_merged_level_maps(custom_levels)
    result: Dict[str, Optional[str]] = {name: None for name in level_names}

    if not taxonomy_string or not taxonomy_string.strip():
        return result

    pairs = taxonomy_string.strip().split(sep)

    for pair in pairs:
        pair = pair.strip()
        if not pair:
            continue

        if "__" not in pair:
            logger.warning(f"Invalid taxonomy pair (missing __): '{pair}'")
            continue

        level_code, _, value = pair.partition("__")

        if level_code not in level_map_b:
            logger.warning(f"Unknown level prefix '{level_code}' in pair '{pair}'")
            continue

        if ";" in value or "__" in value:
            logger.error(
                f"Invalid taxonomy value contains separator or __: '{pair}'"
            )
            continue

        level_name = level_map_b[level_code]
        if value:
            result[level_name] = value
        else:
            logger.debug(f"Empty value for level '{level_name}' (code '{level_code}')")
            result[level_name] = None

    for level_name in level_names:
        if result[level_name] is None and level_name not in _MISSING_LEVELS_LOGGED:
            _MISSING_LEVELS_LOGGED.add(level_name)
            logger.debug(f"Missing level '{level_name}' in table string")

    return result


def parse_taxonomy(
    label: str,
    mode: str = "reverse",
    custom_levels: Optional[str] = None,
) -> Dict[str, Optional[str]]:
    """Public API: parse taxonomy from a label string.

    Auto-detects Format A (embedded) vs Format B (table-style):
      - If label contains '{level}__' patterns, treat as Format B
      - Otherwise, treat as Format A

    Args:
        label: taxonomy string in either format
        mode: parsing mode for Format A ('reverse', 'greedy', 'segment')
        custom_levels: custom level prefixes (e.g. 'kingdom:_k_')

    Returns:
        Standard taxonomy dict with keys domain..species (None for missing).
    """
    if re.search(r"[a-z]__", label):
        return parse_taxonomy_table(label, custom_levels=custom_levels)
    return parse_taxonomy_embedded(label, mode=mode, custom_levels=custom_levels)


def merge_taxonomy(
    embedded: Dict[str, Optional[str]],
    from_table: Dict[str, Optional[str]],
    priority: str = "table",
) -> Dict[str, Optional[str]]:
    """Merge taxonomy from two sources with conflict detection.

    Args:
        embedded: taxonomy parsed from embedded label (Format A)
        from_table: taxonomy from external table (Format B)
        priority: 'embedded' or 'table' — which source takes precedence on conflict

    Returns:
        Merged taxonomy dict.
    """
    result = _empty_taxonomy()

    for level in STANDARD_LEVELS:
        emb_val = embedded.get(level)
        tab_val = from_table.get(level)

        if emb_val is None and tab_val is None:
            continue
        elif emb_val is None:
            result[level] = tab_val
        elif tab_val is None:
            result[level] = emb_val
        elif emb_val == tab_val:
            result[level] = emb_val
        else:
            preferred = emb_val if priority == "embedded" else tab_val
            logger.warning(
                f"Taxonomy conflict at '{level}': embedded='{emb_val}' vs "
                f"table='{tab_val}'. Using '{priority}' priority: '{preferred}'"
            )
            result[level] = preferred

    return result


def _detect_table_separator(file_path: str) -> str:
    """Auto-detect whether a table file uses tab or comma separators."""
    with open(file_path, encoding="utf-8", newline="") as f:
        first_line = f.readline()
    if "\t" in first_line:
        return "\t"
    if "," in first_line:
        return ","
    return "\t"


def _has_level_markers(tax_string: str, taxonomy_format: str = "table") -> bool:
    """Does this cell carry taxonomic level markers at all?

    Used only to recognise a header line. Format B needs a ``{code}__`` pair;
    Format A (embedded segment) needs a ``_x_`` delimiter. Anything else is
    prose, e.g. the ``taxonomy`` column label of a spreadsheet export.
    """
    if not tax_string:
        return False
    if taxonomy_format == "embedded":
        return re.search(r"_[a-z]_", tax_string) is not None
    return "__" in tax_string


def load_taxonomy_table(
    file_path: str,
    sep: Optional[str] = None,
    taxonomy_format: str = "table",
    delimiter_mode: str = "reverse",
    ignore_malformed: bool = False,
    custom_levels: Optional[str] = None,
) -> Dict[str, Dict[str, Optional[str]]]:
    """Load an external taxonomy table file.

    Table format:
      Column 1: tip label (must match tree tip names)
      Column 2: taxonomy string (Format B or Format A segment)
      Separator: auto-detected (tab or comma), or forced via --table-sep

    Args:
        file_path: path to the taxonomy table
        sep: separator override (None for auto-detect)
        taxonomy_format: 'table' (Format B) or 'embedded' (Format A segment)
        delimiter_mode: parsing mode for Format A ('reverse', 'greedy', 'segment')
        ignore_malformed: skip malformed rows instead of terminating
        custom_levels: ``--taxonomy-levels`` value ('name:prefix' pairs, e.g.
            ``'kingdom:k__'``) extending the built-in ranks; without it a table
            carrying a rank the built-in set does not know is read with that
            rank dropped.

    Returns:
        {tip_label: taxonomy_dict} — the dict carries the built-in ranks plus
        any rank named by *custom_levels*.

    Raises:
        TaxonomyConflictError: if malformed row found and ignore_malformed=False
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Taxonomy table not found: {file_path}")

    if path.stat().st_size == 0:
        raise ValueError(f"Taxonomy table is empty: {file_path}")

    if sep is None:
        sep = _detect_table_separator(file_path)

    if custom_levels:
        # Announced once per table, at the level a reader of the run log sees:
        # Which ranks were parsed is a property of the result, not a debug
        # Detail, and silently dropping one is how a kingdom-level table used
        # To read as a domain..species table.
        logger.info(
            f"Extending the built-in ranks while reading {file_path}: "
            f"{custom_levels}"
        )

    result: Dict[str, Dict[str, Optional[str]]] = {}
    malformed_count = 0

    def _rows_with_parse_errors(rdr):
        """Yield ``(row, line_num, parse_error)``.

        The ``csv`` parser itself can refuse a line (NUL byte, unterminated
        quote). CPython raises ``csv.Error`` for an embedded NUL on 3.10–3.12
        but tolerates it on 3.13+, which would let the same taxonomy table
        either crash with a raw parser error or silently reach the
        control-character check below depending on the interpreter. Surfacing
        the failure as a value keeps behaviour identical across the supported
        range and routes it through the module's documented malformed-row
        policy (failures must be loud and distinguishable).
        """
        n = 0
        while True:
            try:
                row = next(rdr)
            except StopIteration:
                return
            except csv.Error as exc:
                n += 1
                yield n, None, str(exc)
                continue
            n += 1
            yield n, row, None

    with open(file_path, encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter=sep)
        csv_error_streak = 0
        for line_num, row, parse_error in _rows_with_parse_errors(reader):
            if parse_error is not None:
                csv_error_streak += 1
                msg = (
                    f"Unparseable line {line_num} in taxonomy table "
                    f"{file_path}: {parse_error}"
                )
                if not ignore_malformed:
                    raise TaxonomyConflictError(msg)
                if csv_error_streak >= 3:
                    # Refuse to truncate the table silently.
                    raise TaxonomyConflictError(
                        msg + " (3rd consecutive parser error — not continuing "
                              "with a partially read taxonomy table)"
                    )
                logger.error(msg + " — skipping (ignore_malformed=True)")
                malformed_count += 1
                continue
            csv_error_streak = 0
            if not row or (len(row) >= 1 and row[0].startswith("#")):
                continue
            if len(row) < 2:
                logger.warning(
                    f"Skipping line {line_num} in taxonomy table: insufficient columns"
                )
                continue

            tip_label = row[0].strip()
            tax_string = row[1].strip()

            # Header sniffing. Every spreadsheet export puts a column label on
            # Line 1 (``genome_id\ttaxonomy``), and the reader used to parse it
            # Like any other row: the only trace was a WARNING about the pair
            # 'taxonomy' missing '__', and a phantom tip named ``genome_id``
            # With an all-None taxonomy entered the mapping — which then joins
            # Every downstream tip-set comparison, including the taxonomic
            # Scope detection that picks the monophyly rank. A first line whose
            # Taxonomy cell is non-empty but carries no level marker at all is
            # A header. (An EMPTY cell is not: that is a real tip with missing
            # Taxonomy, and dropping a tip silently would be worse than the
            # Phantom-tip bug this fixes.)
            if (
                line_num == 1
                and tax_string
                and not _has_level_markers(tax_string, taxonomy_format)
            ):
                logger.info(
                    f"Treating line 1 of taxonomy table {file_path} as a header "
                    f"row (no level markers in {tax_string!r}) — skipped"
                )
                continue

            if not tip_label:
                logger.warning(
                    f"Skipping line {line_num}: empty tip label"
                )
                continue

            # 检查控制字符和 Unicode 双向覆盖字符
            safety_errors = check_control_characters(tip_label, f"taxonomy table line {line_num} tip label")
            safety_errors += check_control_characters(tax_string, f"taxonomy table line {line_num} taxonomy string")
            if safety_errors:
                if ignore_malformed:
                    logger.error(
                        f"Skipping malformed line {line_num}: {safety_errors[:3]}"
                    )
                    malformed_count += 1
                    continue
                else:
                    raise TaxonomyConflictError(
                        f"Unsafe characters in taxonomy table line {line_num}: {safety_errors[:5]}"
                    )

            # 检查分类学循环依赖
            if taxonomy_format == "table":
                parsed = parse_taxonomy_table(tax_string, custom_levels=custom_levels)
                _check_taxonomy_cycle(tip_label, parsed, result, ignore_malformed, line_num)
            else:
                parsed = parse_taxonomy_embedded(
                    tax_string, mode=delimiter_mode, custom_levels=custom_levels)

            result[tip_label] = parsed

    if malformed_count > 0:
        logger.warning(f"Skipped {malformed_count} malformed lines in taxonomy table")

    logger.info(f"Loaded taxonomy for {len(result)} tips from {file_path}")
    return result


def _check_taxonomy_cycle(
    tip_label: str,
    parsed: Dict[str, Optional[str]],
    existing: Dict[str, Dict[str, Optional[str]]],
    ignore_malformed: bool,
    line_num: int,
) -> None:
    """检查分类学循环依赖（如 A→B, B→A）。

    Raises:
        TaxonomyConflictError: if cycle detected and ignore_malformed=False
    """
    # 检查当前 tip 的分类值是否引用了另一个 tip
    for level, value in parsed.items():
        if value and value in existing:
            other_tax = existing[value]
            for other_level, other_value in other_tax.items():
                if other_value and other_value == tip_label:
                    msg = (
                        f"Taxonomy cycle detected at line {line_num}: "
                        f"'{tip_label}' references '{value}', "
                        f"but '{value}' references '{tip_label}'"
                    )
                    if ignore_malformed:
                        logger.error(msg)
                    else:
                        raise TaxonomyConflictError(msg)


def validate_table_tree_consistency(
    table_taxa: Dict[str, Dict],
    tree_tips: Set[str],
) -> Tuple[Set[str], Set[str]]:
    """Check consistency between taxonomy table and tree tip labels.

    Returns:
        (tips_only_in_table, tips_only_in_tree)
    """
    table_labels = set(table_taxa.keys())
    only_in_table = table_labels - tree_tips
    only_in_tree = tree_tips - table_labels

    if only_in_table:
        logger.warning(
            f"Taxonomy table contains {len(only_in_table)} labels not in tree: "
            f"{sorted(list(only_in_table))[:10]}"
        )
    if only_in_tree:
        logger.debug(
            f"Tree contains {len(only_in_tree)} tips without table taxonomy: "
            f"{sorted(list(only_in_tree))[:10]}"
        )

    return only_in_table, only_in_tree


def _get_tips_for_taxon(
    taxonomy_map: Dict[str, Dict[str, Optional[str]]],
    taxon_label: str,
    level: Optional[str] = None,
) -> Set[str]:
    """Collect all tip labels belonging to a given taxon.

    If level is None, searches across all levels.
    """
    tips: Set[str] = set()
    for tip, tax_dict in taxonomy_map.items():
        if level:
            if tax_dict.get(level) == taxon_label:
                tips.add(tip)
        else:
            for lev, val in tax_dict.items():
                if val == taxon_label:
                    tips.add(tip)
    return tips


def is_monophyletic(
    tree_newick: str,
    taxon_label: str,
    taxonomy_map: Dict[str, Dict[str, Optional[str]]],
    level: Optional[str] = None,
) -> bool:
    """Check if a taxon forms a monophyletic group on the tree.

    Algorithm:
      1. Collect all tips belonging to the taxon.
      2. If empty -> ERROR (raise MonophylyError).
      3. Compute MRCA of these tips.
      4. Get all descendant tips of MRCA.
      5. If descendant set == taxon tip set -> monophyletic.
      6. If taxon tips appear in multiple disjoint subtrees -> non-monophyletic.

    Args:
        tree_newick: Newick string of the tree.
        taxon_label: the taxonomic name to check (e.g. 'Cyanobacteriota').
        taxonomy_map: {tip_label: {domain:..., phylum:...,...}}
        level: optional specific level to search (e.g. 'phylum').

    Returns:
        True if monophyletic, False otherwise.

    Raises:
        MonophylyError: if the taxon is not found on the tree at all.
    """
    try:
        from markerfinder.utils.etree import require_ete3

        EteTree = require_ete3().Tree
        have_ete3 = True
    except Exception:
        # Ete3 缺失时改走 utils.etree 纯 Python 单系回退——
        # 只覆盖"检验"，不覆盖"定根"（mad_root 仍需 ete3）。
        have_ete3 = False

    taxon_tips = _get_tips_for_taxon(taxonomy_map, taxon_label, level)

    if not taxon_tips:
        raise MonophylyError(
            f"Taxon '{taxon_label}' not found in taxonomy mapping"
        )

    if not have_ete3:
        from markerfinder.utils import etree as _etree

        all_tree_tips = _etree.tip_set(tree_newick)
        taxon_tips_on_tree = taxon_tips & all_tree_tips
        if not taxon_tips_on_tree:
            raise MonophylyError(
                f"Taxon '{taxon_label}' tips not present on tree"
            )
        return _etree.is_clade_monophyletic(tree_newick, taxon_tips_on_tree)

    t = EteTree(tree_newick, format=1)

    all_tree_tips = {n.name for n in t.get_leaves() if n.name}
    taxon_tips_on_tree = taxon_tips & all_tree_tips

    if not taxon_tips_on_tree:
        raise MonophylyError(
            f"Taxon '{taxon_label}' tips not present on tree"
        )

    mrca = t.get_common_ancestor(list(taxon_tips_on_tree))
    mrca_descendants = {n.name for n in mrca.get_leaves() if n.name}

    if mrca_descendants == taxon_tips_on_tree:
        return True

    return False


def monophyly_proportion(
    tree_newick: str,
    taxonomy_map: Dict[str, Dict[str, Optional[str]]],
    level: str,
) -> Tuple[Optional[float], int, int]:
    """Proportion of *informative* taxa at *level* that the tree groups together.

    For every distinct taxon label at the given rank, only taxa represented by
    at least two tips on the tree are considered (a single tip can trivially
    form a clade and would otherwise inflate the proportion). A taxon enters the
    denominator only when the tips *outside* it also number at least two: with
    one tip on the other side the claim "this taxon is monophyletic" cannot be
    falsified by any topology, so measuring it would turn "not measurable" into a
    number.

    A taxon counts as monophyletic when the tree carries the split
    ``{taxon | rest}`` — equivalently, when the taxon *or* its complement is a
    clade. That is the root-invariant reading the manual documents: the verdict
    must not depend on where ``mad_root`` happened to place the root.

    Args:
        tree_newick: Newick string of the (rooted or unrooted) tree.
        taxonomy_map: ``{tip_label: {domain:..., genus:...,...}}`` — the same
            mapping:func:`load_taxonomy_table` returns.
        level: taxonomic rank to summarize (e.g. ``"genus"``).

    Returns:
        (proportion, n_total_taxa, n_monophyletic) where ``proportion`` is
        ``n_monophyletic / n_total_taxa`` (``None`` when no taxon at *level* has
        an informative >=2-vs->=2 split to be tested on).
    """
    try:
        from markerfinder.utils.etree import require_ete3

        EteTree = require_ete3().Tree
        have_ete3 = True
    except Exception:
        # Ete3 缺失时走 utils.etree 纯 Python 回退，
        # 复用完全相同的 clade-成员判定语义。
        have_ete3 = False

    from collections import defaultdict

    if not have_ete3:
        from markerfinder.utils import etree as _etree

        all_tree_tips = _etree.tip_set(tree_newick)
        # Clade-成员集合：与 ete3 路径的 node_for_clade 同义（内部节点叶集）。
        node_for_clade: Dict[frozenset, object] = {
            clade: True for clade in _etree.clades_of(tree_newick)
        }
    else:
        t = EteTree(tree_newick, format=1)
        all_tree_tips = {n.name for n in t.get_leaves() if n.name}

        # Pre-compute the descendant-leaf set of every internal node exactly once
        # (O(tree)). A taxon's monophyly then reduces to a single dict lookup: the
        # Taxon is monophyletic iff some internal node's leaf clade equals exactly
        # Its tips on the tree — equivalently, the MRCA of those tips has no extra
        # Leaves. This replaces the per-taxon ``t.get_common_ancestor`` call that
        # Previously re-walked the whole tree for every taxon (O(taxa) traversals on
        # Large datasets). The returned proportion is unchanged.
        node_for_clade: Dict[frozenset, object] = {}
        for node in t.traverse():
            if node.is_leaf():
                continue
            clade = frozenset(n.name for n in node.get_leaves() if n.name)
            # A leaf-clade set maps to a unique node; keep the first (deepest) seen.
            node_for_clade.setdefault(clade, node)

    # Group tree-relevant tips by their taxon label at the chosen rank.
    taxon_to_tips: Dict[str, set] = defaultdict(set)
    for tip, tax_dict in taxonomy_map.items():
        val = tax_dict.get(level)
        if val:
            taxon_to_tips[val].add(tip)

    total = 0
    monophyletic = 0
    for label, tips in taxon_to_tips.items():
        tips_on_tree = tips & all_tree_tips
        # Need ≥2 representatives to form a non-trivial clade.
        if len(tips_on_tree) < 2:
            continue
        group = frozenset(tips_on_tree)
        other = frozenset(all_tree_tips) - group
        # A split needs both sides to carry ≥2 tips, otherwise the claim
        # "taxon X is (not) monophyletic" has no resolving power: with a single
        # Tip on the other side, {X | one tip} is a split of *every* unrooted
        # Topology, so nothing about the gene tree can falsify it. Counting such
        # A taxon would encode "not measurable" as a number — in the
        # Root-dependent reading it even produced risk 1.0 for a perfectly
        # Concordant tree. Skipped here; the caller reports the empty denominator
        # As NOT_MEASURABLE.
        if len(other) < 2:
            continue
        total += 1
        # Root-invariant monophyly: the taxon is a clade for *some* rooting, i.e.
        # The gene tree carries the split {taxon | rest}. Testing only "is a
        # Clade of the rooted tree" would make the verdict depend on where the
        # Rooting rule happened to put the root — the very artefact the
        # Companion note in the manual says monophyly must not have.
        if group in node_for_clade or other in node_for_clade:
            monophyletic += 1

    if total == 0:
        return None, 0, 0
    return monophyletic / total, total, monophyletic


def detect_scope_rank(
    taxonomy_map: Dict[str, Dict[str, Optional[str]]],
    tip_labels: Iterable[str],
) -> Optional[str]:
    """Deepest taxonomic rank at which *all* given tips share one label.

    This identifies the taxonomic *scope* of a gene tree: e.g. if every tip
    belongs to the same phylum, the tree's scope is ``"phylum"``. The HGT
    monophyly screen should then be applied one rank *below* the scope (see
:func:`monophyly_proportion`).

    A rank is only considered cohesive when *every* tip has a non-empty label at
    that rank and all labels are identical; the deepest (finest) such rank wins.
    Cohesion is therefore transitive: if tips differ at phylum they also differ
    at every finer rank, so the cohesive ranks always form a prefix of
:data:`STANDARD_LEVELS`.

    Args:
        taxonomy_map: ``{tip: {rank: value,...}}`` mapping.
        tip_labels: iterable of tip names (e.g. genome ids) present on the tree.

    Returns:
        The deepest cohesive rank (one of:data:`STANDARD_LEVELS`), or ``None``
        when the tips are not cohesive at any rank (e.g. they span multiple
        domains, so no single rank from domain..species unifies them).
    """
    tips = [t for t in tip_labels if t in taxonomy_map]
    if not tips:
        return None
    n = len(tips)
    present = {lvl: 0 for lvl in STANDARD_LEVELS}
    labels = {lvl: set() for lvl in STANDARD_LEVELS}
    for tip in tips:
        tax = taxonomy_map[tip]
        for lvl in STANDARD_LEVELS:
            v = tax.get(lvl)
            if v:
                present[lvl] += 1
                labels[lvl].add(v)
    # Deepest (finest) rank first; the first rank with a single label shared by
    # Every tip is the tree's scope.
    for lvl in reversed(STANDARD_LEVELS):
        if present[lvl] == n and len(labels[lvl]) == 1:
            return lvl
    return None


def find_special_identifier(
    tree_newick: str,
    taxonomy_map: Dict[str, Dict[str, Optional[str]]],
    identifier: str,
) -> Optional[object]:
    """Find the node for a special identifier.

    Special identifiers (case-sensitive, uppercase only):
      - root: tree root
      - LUCA: MRCA of all Bacteria + Archaea tips
      - LACA: MRCA of all Archaea tips
      - LBCA: MRCA of all Bacteria tips

    Args:
        tree_newick: Newick string.
        taxonomy_map: taxonomy mapping.
        identifier: one of 'root', 'LUCA', 'LACA', 'LBCA'.

    Returns:
        ete3 Tree node, or None if identifier cannot be resolved.

    Raises:
        MonophylyError: if required domain is missing from the tree.
    """
    try:
        from markerfinder.utils.etree import require_ete3

        EteTree = require_ete3().Tree
        have_ete3 = True
    except Exception:
        # Ete3 缺失时走纯 Python MRCA 回退，返回一个暴露
        # ``get_leaves`` / ``name`` / ``is_leaf`` 的最小节点替身。
        have_ete3 = False

    if identifier != identifier.upper() and identifier != "root":
        raise MonophylyError(
            f"Special identifier '{identifier}' must be UPPERCASE (except 'root')"
        )

    if not have_ete3:
        return _find_special_identifier_fallback(
            tree_newick, taxonomy_map, identifier
        )

    t = EteTree(tree_newick, format=1)

    if identifier == "root":
        return t

    bacteria_tips = _get_tips_for_taxon(taxonomy_map, "Bacteria", "domain")
    archaea_tips = _get_tips_for_taxon(taxonomy_map, "Archaea", "domain")

    all_tree_tips = {n.name for n in t.get_leaves() if n.name}
    bacteria_on_tree = bacteria_tips & all_tree_tips
    archaea_on_tree = archaea_tips & all_tree_tips

    if identifier == "LUCA":
        if not bacteria_on_tree or not archaea_on_tree:
            raise MonophylyError(
                "LUCA requires both Bacteria and Archaea on the tree"
            )
        combined = list(bacteria_on_tree | archaea_on_tree)
        return t.get_common_ancestor(combined)

    if identifier == "LBCA":
        if not bacteria_on_tree:
            raise MonophylyError("LBCA requires Bacteria on the tree")
        return t.get_common_ancestor(list(bacteria_on_tree))

    if identifier == "LACA":
        if not archaea_on_tree:
            raise MonophylyError("LACA requires Archaea on the tree")
        return t.get_common_ancestor(list(archaea_on_tree))

    raise MonophylyError(f"Unknown special identifier: '{identifier}'")


class _FallbackLeaf:
    """Minimal leaf stand-in for the no-ete3 path of find_special_identifier."""

    __slots__ = ("name",)

    def __init__(self, name: str) -> None:
        self.name = name


class _FallbackNode:
    """Minimal MRCA node stand-in (no-ete3 path of find_special_identifier).

    Exposes the subset of the ete3 node API the function's consumers use:
    ``get_leaves``, ``name``, ``is_leaf``.
    """

    def __init__(self, names: Iterable[str]) -> None:
        self._leaves = [_FallbackLeaf(n) for n in sorted(names)]

    def get_leaves(self) -> list:
        return list(self._leaves)

    def is_leaf(self) -> bool:
        return len(self._leaves) <= 1

    @property
    def name(self) -> str:
        return "|".join(leaf.name for leaf in self._leaves)


def _find_special_identifier_fallback(
    tree_newick: str,
    taxonomy_map: Dict[str, Dict[str, Optional[str]]],
    identifier: str,
):
    """Pure-Python MRCA lookup when ete3 is unavailable."""
    from markerfinder.utils import etree as _etree

    all_tree_tips = _etree.tip_set(tree_newick)

    if identifier == "root":
        return _FallbackNode(all_tree_tips)

    def _domain_tips_on_tree(domain: str) -> set:
        tips = _get_tips_for_taxon(taxonomy_map, domain, "domain")
        return tips & all_tree_tips

    if identifier == "LUCA":
        combined = _domain_tips_on_tree("Bacteria") | _domain_tips_on_tree("Archaea")
        if not _domain_tips_on_tree("Bacteria") or not _domain_tips_on_tree("Archaea"):
            raise MonophylyError(
                "LUCA requires both Bacteria and Archaea on the tree"
            )
        return _FallbackNode(combined)

    if identifier == "LBCA":
        bacteria = _domain_tips_on_tree("Bacteria")
        if not bacteria:
            raise MonophylyError("LBCA requires Bacteria on the tree")
        return _FallbackNode(bacteria)

    if identifier == "LACA":
        archaea = _domain_tips_on_tree("Archaea")
        if not archaea:
            raise MonophylyError("LACA requires Archaea on the tree")
        return _FallbackNode(archaea)

    raise MonophylyError(f"Unknown special identifier: '{identifier}'")
