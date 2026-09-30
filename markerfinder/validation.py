"""Input file validation: tree files, sequence files, cross-validation.

Public API:
  load_tree(path, validate=True) -> Tree
  cross_validate(tree_path, seq_path, strict=True) -> ValidationReport
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from markerfinder.exceptions import (
    CrossValidationError,
    PhyloFormatError,
    SequenceValidationError,
    TreeValidationError,
)

logger = logging.getLogger(__name__)

# Unicode bidirectional text override characters (security risk)
BIDI_CONTROLS = {
    "\u202E",  # RIGHT-TO-LEFT OVERRIDE
    "\u202D",  # LEFT-TO-RIGHT OVERRIDE
    "\u202B",  # RIGHT-TO-LEFT EMBEDDING
    "\u202A",  # LEFT-TO-RIGHT EMBEDDING
    "\u202C",  # POP DIRECTIONAL FORMATTING
    "\u2066",  # LEFT-TO-RIGHT ISOLATE
    "\u2067",  # RIGHT-TO-LEFT ISOLATE
    "\u2068",  # FIRST STRONG ISOLATE
    "\u2069",  # POP DIRECTIONAL ISOLATE
    "\u200E",  # LEFT-TO-RIGHT MARK
    "\u200F",  # RIGHT-TO-LEFT MARK
}

# Control characters that are not allowed (except common whitespace)
_ALLOWED_WHITESPACE = {"\n", "\r", "\t", " "}


def check_control_characters(text: str, field_name: str = "input") -> List[str]:
    """检查文本中是否包含控制字符或 Unicode 双向覆盖字符。

    Args:
        text: 要检查的字符串。
        field_name: 字段名称（用于错误消息）。

    Returns:
        错误消息列表（空列表表示无问题）。
    """
    errors: List[str] = []

    for i, ch in enumerate(text):
        # 检查 Unicode 双向覆盖字符
        if ch in BIDI_CONTROLS:
            errors.append(
                f"Unicode bidi override character U+{ord(ch):04X} found in {field_name} at position {i}"
            )
            continue

        # 检查控制字符（排除允许的空白符）
        if ord(ch) < 0x20 and ch not in _ALLOWED_WHITESPACE:
            errors.append(
                f"Control character U+{ord(ch):04X} found in {field_name} at position {i}"
            )
        elif ord(ch) == 0x7F:  # DEL
            errors.append(
                f"Control character U+{ord(ch):04X} (DEL) found in {field_name} at position {i}"
            )

    return errors


def validate_string_safety(text: str, field_name: str = "input") -> bool:
    """验证字符串安全性，拒绝包含控制字符或双向覆盖字符的输入。

    Args:
        text: 要验证的字符串。
        field_name: 字段名称。

    Returns:
        True 如果字符串安全。

    Raises:
        SequenceValidationError: 如果发现危险字符。
    """
    errors = check_control_characters(text, field_name)
    if errors:
        raise SequenceValidationError(
            f"Unsafe characters detected in {field_name}: {errors[:5]}"
        )
    return True


class ValidationReport:
    """Result of a cross-validation between tree and sequence file."""

    def __init__(self) -> None:
        self.tree_tips: Set[str] = set()
        self.seq_ids: Set[str] = set()
        self.only_in_tree: Set[str] = set()
        self.only_in_seq: Set[str] = set()
        self.is_valid: bool = True
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def __repr__(self) -> str:
        return (
            f"ValidationReport(valid={self.is_valid}, "
            f"tree_tips={len(self.tree_tips)}, seq_ids={len(self.seq_ids)}, "
            f"only_in_tree={len(self.only_in_tree)}, only_in_seq={len(self.only_in_seq)})"
        )


def validate_newick_string(newick: str) -> Tuple[List[str], List[str]]:
    """Validate a Newick string for structural correctness.

    Checks:
      - Parenthesis balance
      - Non-negative branch lengths
      - No empty node names for tips
      - Duplicate tip names
      - No self-referencing

    Returns:
        A 2-tuple ``(errors, warnings)``.
        ``errors`` are hard structural problems (unbalanced parens, negative
        branch lengths, empty/duplicate tip names) that should block tree use.
        ``warnings`` are non-fatal advisories (e.g. duplicate internal-node
        names). Callers that only need a yes/no verdict should check ``errors``;
        do NOT treat ``warnings`` as errors.
    """
    errors: List[str] = []
    warnings: List[str] = []

    if not newick or not newick.strip():
        errors.append("Empty Newick string")
        return errors, warnings

    depth = 0
    for i, ch in enumerate(newick):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                errors.append(f"Unbalanced parenthesis at position {i}")
                break
    if depth != 0:
        errors.append(f"Unbalanced parenthesis: final depth={depth}")

    # Negative branch lengths. Support scientific notation (e.g.:-1e-3) as well
    # As plain decimals (-0.5).
    neg_matches = re.findall(r":-\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", newick)
    if neg_matches:
        errors.append(f"Negative branch lengths found: {neg_matches[:5]}")

    # Newick/NHX quoted labels: a label may be wrapped in single or double quotes
    # So it can contain '(' ')' ',' ':' etc. (e.g. 'Taxon (A)', "Node 1"). The
    # Whole quoted token is captured as one label so duplicate/empty detection
    # And quoted-label parsing work correctly.
    _QUOTED_LABEL = r"'(?:[^']*)'|\"(?:[^\"]*)\""
    tip_pattern = re.findall(
        r"(?:\(|,)\s*(" + _QUOTED_LABEL + r"|[^(),:;\s]*)\s*(?=[,:)\]])",
        newick,
    )
    # 检测冒号前为空的叶名（如 "(0.1" 或 ",:0.2"），正则无法捕获，单独处理
    empty_colon_tips = re.findall(r"(?:\(|,)\s*:", newick)
    if empty_colon_tips:
        errors.append(
            f"Empty tip name found in tree ({len(empty_colon_tips)} occurrence(s))"
        )
    tips = [t.strip() for t in tip_pattern if t.strip()]
    seen: Dict[str, int] = {}
    for tip in tips:
        if not tip:
            errors.append("Empty tip name found in tree")
            continue
        seen[tip] = seen.get(tip, 0) + 1

    duplicates = {name: count for name, count in seen.items() if count > 1}
    if duplicates:
        errors.append(f"Duplicate tip names: {duplicates}")

    # Internal-node labels follow ')' and precede ',' '(' or ';' (i.e. node
    # Names that are NOT followed by a branch length). Quoted labels supported.
    internal_pattern = re.findall(
        r"\)\s*(" + _QUOTED_LABEL + r"|[^(),;]+?)\s*(?=[(),;])",
        newick,
    )
    internal_seen: Dict[str, int] = {}
    for name in internal_pattern:
        name = name.strip()
        if name and re.match(r"^-?\d", name):
            continue
        if name:
            internal_seen[name] = internal_seen.get(name, 0) + 1
    internal_dup = {n: c for n, c in internal_seen.items() if c > 1}
    if internal_dup:
        warnings.append(f"Duplicate internal node names: {internal_dup}")

    return errors, warnings


def detect_tree_format(file_path: str) -> str:
    """Detect tree file format: 'newick', 'nexus', or 'nhx'.

    Based on file extension and content.
    """
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in (".nex", ".nexus"):
        return "nexus"
    if suffix == ".nhx":
        return "nhx"

    try:
        with open(file_path, encoding="utf-8", newline="") as f:
            first_line = f.readline(1024).strip()
    except (FileNotFoundError, OSError):
        return "newick"

    if first_line.lower().startswith("#nexus"):
        return "nexus"
    if "&&NHX" in first_line:
        return "nhx"

    return "newick"


def count_trees_in_file(file_path: str, fmt: str = "newick") -> int:
    """Count the number of trees in a file.

    For Nexus: count TREE/TREE_* statements.
    For Newick: count non-empty lines ending with semicolons.
    """
    try:
        with open(file_path, encoding="utf-8", newline="") as f:
            content = f.read()
    except (FileNotFoundError, OSError):
        # — legitimate business value, not a disguised failure: this
        # Function counts trees, and an unreadable file contributes zero trees
        # To the count. Its only caller (cli/validation.py) compares the number
        # Against the expected tree count and reports the file as unusable, so
        # 0 leads to a visible error rather than a silent pass.
        return 0

    if fmt == "nexus":
        tree_blocks = re.findall(r"^\s*tree\s+\S+", content, re.IGNORECASE | re.MULTILINE)
        return len(tree_blocks)

    lines = content.strip().split("\n")
    count = 0
    for line in lines:
        line = line.strip()
        if line and line.endswith(";"):
            count += 1
    # Return the true count (0 is a legitimate "no trees" result). The caller is
    # Responsible for deciding how to handle an empty tree file; previously this
    # Returned ``max(count, 1)`` which silently reported an empty file as
    # Containing one tree and let downstream parsing crash.
    return count


def split_newick_trees(text: str) -> List[str]:
    """Split multi-Newick content into individual tree strings.

    Quote-aware: semicolons inside single/double-quoted labels do not split.
    Returns only non-empty segments, each terminated with ``;``.
    """
    trees: List[str] = []
    depth = 0
    start = 0
    in_quote = False
    quote_char = ""
    for i, ch in enumerate(text):
        if in_quote:
            if ch == quote_char:
                in_quote = False
            continue
        if ch in ("'", '"'):
            in_quote = True
            quote_char = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        elif ch == ";" and depth <= 0:
            seg = text[start:i + 1].strip()
            if seg:
                trees.append(seg)
            start = i + 1
    tail = text[start:].strip()
    if tail:
        trees.append(tail if tail.endswith(";") else tail + ";")
    return trees


def _extract_trees(content: str, fmt: str) -> List[str]:
    """Extract individual tree strings from file content (Nexus or Newick)."""
    if fmt == "nexus":
        matches = re.findall(
            r"tree\s+\S+\s*=\s*([;()\[\]A-Za-z0-9_.:,{}\/\-]+);",
            content,
            re.IGNORECASE,
        )
        return [m.strip() + ";" for m in matches if m.strip()]
    return split_newick_trees(content)


def parse_nhx_annotations(newick: str) -> Dict[str, Dict[str, str]]:
    """Parse NHX annotations from a Newick string.

    Returns:
        {node_name_or_clade: {nhx_key: nhx_value}}
    """
    annotations: Dict[str, Dict[str, str]] = {}
    nhx_pattern = re.findall(r"\[&&NHX:([^\]]+)\]", newick)
    for i, nhx_str in enumerate(nhx_pattern):
        props: Dict[str, str] = {}
        for pair in nhx_str.split(":"):
            if "=" in pair:
                k, _, v = pair.partition("=")
                props[k] = v
        annotations[f"nhx_block_{i}"] = props
    return annotations


def strip_nhx_annotations(newick: str) -> str:
    """Remove NHX annotation blocks from a Newick string."""
    return re.sub(r"\[&&NHX:[^\]]*\]", "", newick)


def _validate_newick_content(
    newick: str,
    strip_annotations: bool = False,
) -> Tuple[str, List[str]]:
    """Validate a single Newick string; raise:class:`TreeValidationError` on critical errors.

    Shared by:func:`validate_tree_file` (whole-file content) and
:func:`load_tree` (per-tree selection via ``tree_index``) so both paths
    apply identical annotation-stripping, safety and criticality rules.

    Critical (raising) errors: unbalanced parentheses, negative branch lengths,
    duplicate tip names, and empty tip names — anything else is returned as
    soft ``errors + warnings`` for the caller to log.
    """
    if not newick.endswith(";"):
        newick += ";"

    if strip_annotations:
        newick = strip_nhx_annotations(newick)

    safety_errors = check_control_characters(newick, "tree file content")
    if safety_errors:
        raise TreeValidationError(
            f"Unsafe characters in tree file: {safety_errors[:5]}"
        )

    errors, warnings = validate_newick_string(newick)

    # 括号不平衡与空 tip 名同样是结构性硬伤: 放行只会把失败推迟到 ete3
    # 内部并产生难解的报错, 因此与负分支/重复 tip 一并列为 critical.
    critical_markers = ("Negative", "Duplicate tip", "Unbalanced", "Empty tip")
    critical_errors = [e for e in errors if any(m in e for m in critical_markers)]
    if critical_errors:
        raise TreeValidationError(
            f"Critical tree validation errors: {critical_errors}"
        )

    return newick, errors + warnings


def validate_tree_file(
    file_path: str,
    strip_annotations: bool = False,
) -> Tuple[str, List[str]]:
    """Validate a tree file and return its content.

    Args:
        file_path: path to the tree file.
        strip_annotations: remove NHX/bootstrap annotations.

    Returns:
        (newick_content, list_of_error_strings)

    Raises:
        TreeValidationError: if validation fails with CRITICAL errors.
        FileNotFoundError: if file does not exist.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Tree file not found: {file_path}")

    if path.is_dir():
        raise TreeValidationError(f"Expected a file, got directory: {file_path}")

    if path.stat().st_size == 0:
        raise TreeValidationError(f"Tree file is empty: {file_path}")

    fmt = detect_tree_format(file_path)

    with open(file_path, encoding="utf-8", newline="") as f:
        content = f.read().strip()

    if fmt == "nexus":
        tree_match = re.search(
            r"tree\s+\S+\s*=\s*([;()\[\]A-Za-z0-9_.:,{}\/\-]+);",
            content,
            re.IGNORECASE,
        )
        if tree_match:
            newick = tree_match.group(1) + ";"
        else:
            raise TreeValidationError(
                f"No valid tree found in Nexus file: {file_path}"
            )
    else:
        newick = content

    return _validate_newick_content(newick, strip_annotations)


def load_tree(
    path: str,
    validate: bool = True,
    strip_annotations: bool = False,
    tree_index: Optional[int] = None,
) -> "object":
    """Public API: load and optionally validate a tree file.

    Args:
        path: path to tree file.
        validate: run deep validation (default True).
        strip_annotations: remove NHX/bootstrap annotations.
        tree_index: 0-based index of the tree to load when the file holds
            multiple trees (``--multi-tree-mode first/last/random`` resolves to
            an index before calling). ``None`` keeps the historical behaviour of
            reading the whole file content as one tree.

    Returns:
        markerfinder.models.tree.Tree object.

    Raises:
        TreeValidationError: on critical validation failures or an
            out-of-range ``tree_index``.
        FileNotFoundError: if file does not exist.
    """
    from markerfinder.models.tree import Tree

    if tree_index is not None:
        tree_path = Path(path)
        if not tree_path.exists():
            raise FileNotFoundError(f"Tree file not found: {path}")
        if tree_path.stat().st_size == 0:
            raise TreeValidationError(f"Tree file is empty: {path}")
        fmt = detect_tree_format(path)
        with open(path, encoding="utf-8", newline="") as f:
            content = f.read()
        trees = _extract_trees(content.strip(), fmt)
        if not trees:
            raise TreeValidationError(f"No trees found in {path}")
        if tree_index < 0 or tree_index >= len(trees):
            raise TreeValidationError(
                f"tree_index {tree_index} out of range: file holds {len(trees)} tree(s)"
            )
        if validate:
            newick, errors = _validate_newick_content(trees[tree_index], strip_annotations)
            for err in errors:
                if "Duplicate" in err:
                    logger.warning(f"Tree validation warning: {err}")
                else:
                    logger.info(f"Tree validation: {err}")
            return Tree(newick=newick)
        return Tree(newick=trees[tree_index])

    if validate:
        newick, errors = validate_tree_file(path, strip_annotations)
        for err in errors:
            if "Duplicate" in err:
                logger.warning(f"Tree validation warning: {err}")
            else:
                logger.info(f"Tree validation: {err}")
        return Tree(newick=newick)

    with open(path, encoding="utf-8", newline="") as f:
        newick = f.read().strip()

    return Tree(newick=newick)


def detect_sequence_format(file_path: str) -> str:
    """Detect sequence file format: 'fasta', 'fastq', 'phylip', or 'stockholm'.

    Based on file extension and content.
    """
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in (".fq", ".fastq"):
        return "fastq"

    if suffix in (".phylip", ".phy", ".ph"):
        return "phylip"

    if suffix in (".sto", ".stockholm"):
        return "stockholm"

    try:
        with open(file_path, encoding="utf-8", newline="") as f:
            first_line = f.readline(1024).strip()
    except (FileNotFoundError, OSError):
        return "fasta"

    # Stockholm 格式检测
    if first_line.startswith("# STOCKHOLM"):
        return "stockholm"

    # PHYLIP 格式检测：首行为两个整数（序列数 长度）
    if re.match(r"^\s*\d+\s+\d+\s*$", first_line):
        return "phylip"

    if first_line.startswith("@"):
        return "fastq"
    if first_line.startswith(">"):
        return "fasta"

    return "fasta"


def validate_sequence_file(
    file_path: str,
    mol_type: Optional[str] = None,
    skip_length_check: bool = False,
) -> Tuple[Dict[str, str], List[str]]:
    """Validate a sequence file (FASTA/FASTQ).

    Args:
        file_path: path to sequence file.
        mol_type: 'DNA', 'RNA', or 'protein' for alphabet checking.
                   None for auto-detect.
        skip_length_check: skip alignment length consistency check.

    Returns:
        ({seq_id: sequence_string}, list_of_error_strings)

    Raises:
        SequenceValidationError: on critical failures (duplicate IDs, etc.)
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Sequence file not found: {file_path}")

    if path.is_dir():
        raise SequenceValidationError(
            f"Expected a file, got directory: {file_path}"
        )

    if path.stat().st_size == 0:
        raise SequenceValidationError(f"Sequence file is empty: {file_path}")

    fmt = detect_sequence_format(file_path)

    if fmt in ("phylip", "stockholm"):
        raise SequenceValidationError(
            f"Unsupported sequence format '{fmt}' detected in {file_path}. "
            f"Please convert to FASTA format before running MarkerFinder."
        )
    sequences: Dict[str, str] = {}
    errors: List[str] = []
    lengths: Dict[str, int] = {}

    seen_ids: Dict[str, int] = {}
    with open(file_path, encoding="utf-8", newline="") as f:
        current_id = ""
        current_seq: List[str] = []
        line_num = 0

        for line in f:
            line_num += 1
            line = line.rstrip("\n\r")

            if fmt == "fastq":
                if line_num % 4 == 1:
                    if not line.startswith("@"):
                        errors.append(f"Line {line_num}: expected '@' header in FASTQ")
                    current_id = line[1:].split()[0] if line else ""
                elif line_num % 4 == 2:
                    current_seq = [line]
                else:
                    continue
            elif line.startswith(">"):
                if current_id:
                    seq_str = "".join(current_seq)
                    sequences[current_id] = seq_str
                    lengths[current_id] = len(seq_str)
                current_id = line[1:].split()[0] if len(line) > 1 else ""
                current_seq = []
                if current_id:
                    safety_errors = check_control_characters(current_id, "sequence ID")
                    if safety_errors:
                        raise SequenceValidationError(
                            f"Unsafe characters in sequence ID at line {line_num}: {safety_errors[:3]}"
                        )
                    seen_ids[current_id] = seen_ids.get(current_id, 0) + 1
            else:
                current_seq.append(line.strip())

        if current_id:
            seq_str = "".join(current_seq)
            sequences[current_id] = seq_str
            lengths[current_id] = len(seq_str)

    duplicate_ids = [sid for sid, count in seen_ids.items() if count > 1]
    if duplicate_ids:
        raise SequenceValidationError(
            f"Duplicate sequence IDs: {duplicate_ids}"
        )

    if mol_type is None and sequences:
        mol_type = _auto_detect_mol_type(sequences)
        logger.info(f"Auto-detected molecule type: {mol_type}")

    if mol_type and sequences:
        alphabet_errors = _validate_alphabet(sequences, mol_type)
        errors.extend(alphabet_errors)

    if not skip_length_check and len(set(lengths.values())) > 1:
        length_vals = list(lengths.values())
        errors.append(
            f"Sequence lengths not consistent: min={min(length_vals)}, max={max(length_vals)}"
        )
        logger.warning("Sequence length inconsistency detected")

    return sequences, errors


DNA_ALPHABET = set("ATGCNRYKMSWBDHV.-")
RNA_ALPHABET = set("AUGCNRYKMSWBDHV.-")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBJZX.*-")


def _auto_detect_mol_type(sequences: Dict[str, str]) -> str:
    """Auto-detect molecule type based on character distribution."""
    all_chars = set()
    total_len = 0
    for seq in list(sequences.values())[:10]:
        all_chars.update(seq.upper())
        total_len += len(seq)

    if not all_chars or total_len == 0:
        return "protein"

    dna_only = all_chars - DNA_ALPHABET
    rna_only = all_chars - RNA_ALPHABET
    prot_only = all_chars - PROTEIN_ALPHABET

    if not dna_only:
        return "DNA"
    if not rna_only:
        return "RNA"
    return "protein"


def _validate_alphabet(
    sequences: Dict[str, str], mol_type: str
) -> List[str]:
    """Validate sequences against the expected alphabet."""
    errors: List[str] = []

    if mol_type == "DNA":
        valid = DNA_ALPHABET
    elif mol_type == "RNA":
        valid = RNA_ALPHABET
    elif mol_type == "protein":
        valid = PROTEIN_ALPHABET
    else:
        return errors

    for seq_id, seq in sequences.items():
        invalid = set(seq.upper()) - valid
        if invalid:
            errors.append(
                f"Invalid characters in '{seq_id}' for {mol_type}: {sorted(invalid)}"
            )

    return errors


def cross_validate(
    tree_path: str,
    seq_path: str,
    strict: bool = True,
    mol_type: Optional[str] = None,
    skip_length_check: bool = False,
) -> ValidationReport:
    """Public API: cross-validate tree tip labels against sequence IDs.

    Args:
        tree_path: path to tree file.
        seq_path: path to sequence file.
        strict: if True (default), mismatches cause ERROR.
        mol_type: the molecule type the user forced with ``--mol-type``.
            ``None`` (not given) keeps auto-detection. When it IS given, an
            alphabet violation is an error, not a warning: the option exists to
            assert an alphabet, and asserting one that the data breaks has to
            stop the run.
        skip_length_check: forward ``--skip-length-check`` to the sequence
            validator instead of letting it sit unused at this call site.

    Returns:
        ValidationReport with mismatch details.
    """
    report = ValidationReport()

    try:
        tree_obj = load_tree(tree_path, validate=True)
        report.tree_tips = set(tree_obj.get_tips())
    except Exception as e:
        report.errors.append(f"Failed to load tree: {e}")
        report.is_valid = False
        return report

    try:
        seqs, seq_errors = validate_sequence_file(
            seq_path, mol_type=mol_type, skip_length_check=skip_length_check)
        report.seq_ids = set(seqs.keys())
        if mol_type:
            logger.info(
                f"Sequence alphabet asserted from --mol-type: {mol_type} "
                f"({len(seqs)} sequence(s) checked)."
            )
        for err in seq_errors:
            if mol_type:
                report.errors.append(f"Sequence: {err}")
            else:
                report.warnings.append(f"Sequence: {err}")
    except Exception as e:
        report.errors.append(f"Failed to load sequences: {e}")
        report.is_valid = False
        return report

    if mol_type:
        # The alphabet assertion has to bite: ``--mol-type`` exists to state what
        # The input must be, and a run that reports the violation while carrying
        # On would treat a stated precondition as advice.
        violated = [e for e in report.errors if e.startswith("Sequence:")]
        if violated:
            report.is_valid = False
            if strict:
                raise CrossValidationError(
                    f"Molecule-type assertion failed (--mol-type {mol_type}): "
                    + "; ".join(violated[:3])
                )

    report.only_in_tree = report.tree_tips - report.seq_ids
    report.only_in_seq = report.seq_ids - report.tree_tips

    if report.only_in_tree or report.only_in_seq:
        msg_parts = []
        if report.only_in_tree:
            msg_parts.append(
                f"Only in tree ({len(report.only_in_tree)}): "
                f"{sorted(list(report.only_in_tree))[:10]}"
            )
        if report.only_in_seq:
            msg_parts.append(
                f"Only in sequences ({len(report.only_in_seq)}): "
                f"{sorted(list(report.only_in_seq))[:10]}"
            )
        msg = "Tree-sequence label mismatch. " + "; ".join(msg_parts)
        report.errors.append(msg)

        if strict:
            report.is_valid = False
            raise CrossValidationError(msg)
        else:
            report.warnings.append(msg)

    return report
