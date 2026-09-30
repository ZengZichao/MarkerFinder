"""Self-test runner: verify dependencies and built-in examples.

The historical ``_run_self_test`` function (≈266 lines) is split into one
``_test_*`` helper per scenario group. Each helper returns a list of
``(name, status, detail)`` tuples; ``_run_self_test`` only aggregates them,
prints the table, and returns the overall exit code. The set of scenarios,
their order, and the printed output are byte-for-byte identical to the
original monolithic implementation — only the internal structure changed.
"""

from __future__ import annotations

import logging
import sys
from typing import List

from markerfinder.exceptions import PhyloToolUnavailable

from markerfinder.banner import print_banner
from markerfinder.cli.constants import EXIT_SUCCESS, EXIT_RUNTIME_ERROR

logger = logging.getLogger(__name__)

# Mirrors ``requires-python = ">=3.10,<3.13"`` in pyproject.toml.
# Running outside the declared range is an environment error and must be
# Visible in ``--check`` ( ②).
SUPPORTED_PY_MIN = (3, 10)
SUPPORTED_PY_MAX_EXCLUSIVE = (3, 13)


def interpreter_in_range(
    version_tuple: tuple,
    *,
    min_version: tuple = SUPPORTED_PY_MIN,
    max_exclusive: tuple = SUPPORTED_PY_MAX_EXCLUSIVE,
) -> bool:
    """True iff the interpreter version lies in the declared support range."""
    v = tuple(int(x) for x in version_tuple[:2])
    return min_version <= v < max_exclusive


def _version_tuple(version: str) -> tuple:
    """把版本号字符串解析成可比较的整数元组.

    字符串序比较会得出 "1.100" < "1.99" 这类错误结论; 取数字段构造
    元组后按数值逐段比较. 解析失败的段记 0.
    """
    import re as _re
    parts = _re.findall(r"\d+", str(version))
    return tuple(int(p) for p in parts[:3]) if parts else (0,)


def _test_dependencies() -> List[tuple]:
    results: List[tuple] = []
    # Interpreter version must lie inside the declared support range.
    v = sys.version_info
    # Index (not attribute) access so the check also works with a plain
    # (major, minor, micro) tuple — keeps the must-fail control patchable.
    triplet = (v[0], v[1], v[2])
    in_range = interpreter_in_range(triplet)
    results.append((
        f"Interpreter {triplet[0]}.{triplet[1]}.{triplet[2]} in >=3.10,<3.13",
        "PASS" if in_range else "FAIL",
        "supported" if in_range else (
            "out of range: pyproject.toml requires-python = \">=3.10,<3.13\"; "
            "create a 3.10–3.12 environment (conda env create -f environment.yml)"
        ),
    ))
    # The declared runtime imports, version-checked. A library the package does
    # not use is not probed here: --check reports what a run needs, nothing else.
    deps = [
        ("biopython", "Bio", "1.81"),
    ]
    for name, module, min_ver in deps:
        try:
            mod = __import__(module)
            ver = getattr(mod, "__version__", "0.0")
            passed = _version_tuple(ver) >= _version_tuple(min_ver)
            results.append((
                f"Import {name} (>= {min_ver})",
                "PASS" if passed else "FAIL", ver,
            ))
        except ImportError:
            results.append((f"Import {name} (>= {min_ver})", "FAIL", "not installed"))

    try:
        # Single sanctioned entry point.
        from markerfinder.utils.etree import require_ete3

        EteTree = require_ete3().Tree
        t = EteTree("((A,B),(C,D));")
        results.append(("ete3 Tree parse", "PASS", "ok"))
    except (ImportError, PhyloToolUnavailable):
        # The wording of this row is a documented contract, so the migration to
        # require_ete3() kept it verbatim: the dependency is missing (reported
        # honestly), and the pure-Python split-set fallback keeps
        # RF/quartet/monophyly measurable — the reader must see that.
        results.append((
            "ete3 import", "FAIL",
            "not installed — pure-Python split-set fallback active for "
            "RF/quartet/monophyly (method='splits-python')",
        ))
    except Exception as e:
        results.append(("ete3 Tree parse", "FAIL", str(e)))
    return results


def _test_external_tools() -> List[tuple]:
    """Disclose every external binary the package can invoke.

    CheckM is the reason this exists: ``mag_optimization._run_checkm`` calls
    ``checkm lineage_wf`` and falls back to default quality estimates when it is
    absent, but the tool was not in ``EXTERNAL_TOOLS``, so nothing before the run
    said the completeness/contamination numbers might be estimated.

    Every row is INFO, never PASS/FAIL. The scored ``--check`` item count is a
    quoted number (README/MANUAL state it and
    ``tests/unit/test_docs_numbers_are_current.py`` pins it), and ``PATH`` is
    per-machine -- scoring it would make the documented count drift with the
    host. Criticality is still enforced where it belongs: a missing ``hmmsearch``
    hard-fails at pipeline construction.
    """
    from markerfinder.utils.dependency_check import EXTERNAL_TOOLS, check_external_tools

    results: List[tuple] = []
    try:
        availability = check_external_tools()
    except Exception as exc:  # A PATH probe must never break --check itself
        return [("External tool survey", "INFO", f"not executed ({exc})")]

    for name, purpose, is_critical in EXTERNAL_TOOLS:
        found = availability.get(name, False)
        kind = "critical" if is_critical else "optional"
        if found:
            detail = "on PATH"
        elif is_critical:
            detail = f"NOT FOUND — the run fails at startup. Purpose: {purpose}"
        elif name == "checkm":
            detail = (
                f"not found — Phase 0 falls back to default quality estimates "
                f"(the report labels them quality_source='default_no_checkm'). "
                f"Purpose: {purpose}"
            )
        else:
            detail = f"not found — that step degrades or is skipped. Purpose: {purpose}"
        results.append((f"External tool: {name} ({kind})", "INFO", detail))
    return results


def _test_format_a() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.taxonomy import parse_taxonomy

        r = parse_taxonomy("GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_Cyanobacteriia_o_Cyanobacteriales_f_Prochloraceae_g_Prochloron", mode="reverse")
        results.append(("Format A parse (reverse)", "PASS" if r.get("domain") == "Bacteria" else "FAIL", str(r.get("domain"))))

        # 缺失级别
        r2 = parse_taxonomy("GB_GCA_001_d_Bacteria_p_Firmicutes", mode="reverse")
        results.append(("Format A parse (missing levels)", "PASS" if r2.get("domain") == "Bacteria" and r2.get("family") is None else "FAIL", f"d={r2.get('domain')}, f={r2.get('family')}"))

        # 含下划线分类名
        r3 = parse_taxonomy("ID_d_Bacteria_p_Proteo_bacteria_c_Gamma", mode="reverse")
        results.append(("Format A parse (underscore in name)", "PASS" if r3.get("domain") == "Bacteria" else "FAIL", str(r3.get("domain"))))

        # Greedy 模式
        r4 = parse_taxonomy("ID_d_Bacteria_p_Firmicutes_c_Bacilli", mode="greedy")
        results.append(("Format A parse (greedy)", "PASS" if r4.get("domain") == "Bacteria" else "FAIL", str(r4.get("domain"))))
    except Exception as e:
        results.append(("Format A parse", "FAIL", str(e)))
    return results


def _test_format_b() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.taxonomy import parse_taxonomy_table

        r = parse_taxonomy_table("d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Cyanobacteriales;f__Prochloraceae;g__Prochloron;s__")
        results.append(("Format B parse (normal)", "PASS" if r.get("domain") == "Bacteria" else "FAIL", str(r.get("domain"))))

        # 缺失值
        r2 = parse_taxonomy_table("d__Bacteria;p__Cyanobacteriota;s__")
        results.append(("Format B parse (missing values)", "PASS" if r2.get("domain") == "Bacteria" and r2.get("species") is None else "FAIL", f"d={r2.get('domain')}, s={r2.get('species')}"))

        # 分隔符冲突（值中包含分号）
        r3 = parse_taxonomy_table("d__Bacteria;p__C Firm;icutes")
        results.append(("Format B parse (separator in value)", "PASS" if r3.get("domain") == "Bacteria" else "FAIL", str(r3.get("domain"))))
    except Exception as e:
        results.append(("Format B parse", "FAIL", str(e)))
    return results


def _test_merge_taxonomy() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.taxonomy import merge_taxonomy, parse_taxonomy_embedded, parse_taxonomy_table

        emb = parse_taxonomy_embedded("ID_d_Bacteria_p_Firmicutes")
        tab = parse_taxonomy_table("d__Archaea;p__Proteobacteria")
        merged = merge_taxonomy(emb, tab, priority="table")
        results.append(("Merge taxonomy (table priority)", "PASS" if merged.get("domain") == "Archaea" else "FAIL", str(merged.get("domain"))))

        merged2 = merge_taxonomy(emb, tab, priority="embedded")
        results.append(("Merge taxonomy (embedded priority)", "PASS" if merged2.get("domain") == "Bacteria" else "FAIL", str(merged2.get("domain"))))
    except Exception as e:
        results.append(("Merge taxonomy", "FAIL", str(e)))
    return results


def _test_monophyly() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.taxonomy import is_monophyletic

        tree_nwk = "((B1:0.1,B2:0.2),(A1:0.3,A2:0.4));"
        tax_map = {
            "B1": {"domain": "Bacteria"}, "B2": {"domain": "Bacteria"},
            "A1": {"domain": "Archaea"}, "A2": {"domain": "Archaea"},
        }
        result = is_monophyletic(tree_nwk, "Bacteria", tax_map, "domain")
        results.append(("Monophyly (Bacteria, monophyletic)", "PASS" if result is True else "FAIL", str(result)))

        tree_non_mono = "((B1:0.1,A1:0.2),(B2:0.3,A2:0.4));"
        result2 = is_monophyletic(tree_non_mono, "Bacteria", tax_map, "domain")
        results.append(("Monophyly (Bacteria, non-monophyletic)", "PASS" if result2 is False else "FAIL", str(result2)))
    except ImportError:
        # Is_monophyletic internally falls back to the pure-Python
        # Etree path when ete3 is missing; reaching ImportError here means
        # The package itself is broken, which is a FAIL, not a SKIP.
        results.append(("Monophyly check", "FAIL", "markerfinder.taxonomy not importable"))
    except Exception as e:
        # Both measurement paths (ete3 + pure-Python fallback) failed.
        results.append(("Monophyly check", "FAIL", str(e)))
    return results


def _test_newick_validation() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.validation import validate_newick_string

        errs, _ = validate_newick_string("((A:0.1,B:0.2),(C:0.3,D:0.4));")
        results.append(("Newick validation (valid)", "PASS" if not any("Negative" in e or "Duplicate" in e for e in errs) else "FAIL", str(errs)))

        errs2, _ = validate_newick_string("((A:-1,B),(C,D));")
        results.append(("Newick validation (negative branch)", "PASS" if any("Negative" in e for e in errs2) else "FAIL", "detected" if any("Negative" in e for e in errs2) else "missed"))

        errs3, _ = validate_newick_string("((A:0.1,B:0.2),(A:0.3,D:0.4));")
        results.append(("Newick validation (duplicate tip)", "PASS" if any("Duplicate" in e for e in errs3) else "FAIL", "detected" if any("Duplicate" in e for e in errs3) else "missed"))

        errs4, _ = validate_newick_string("((:0.1,B:0.2),(C:0.3,D:0.4));")
        results.append(("Newick validation (empty tip name)", "PASS" if any("Empty" in e for e in errs4) else "FAIL", "detected" if any("Empty" in e for e in errs4) else "missed"))
    except Exception as e:
        results.append(("Tree validation", "FAIL", str(e)))
    return results


def _test_sequence_validation() -> List[tuple]:
    results: List[tuple] = []
    try:
        import tempfile as _tempfile
        import os as _os
        from markerfinder.validation import validate_sequence_file

        # 正常 FASTA
        fd, fasta_path = _tempfile.mkstemp(suffix=".fasta")
        _os.close(fd)
        _os.chmod(fasta_path, 0o600)
        from pathlib import Path as _Path
        _Path(fasta_path).resolve().write_text(">seq1\nATGCGTACGT\n>seq2\nATGCGTACGT\n", encoding="utf-8")
        seqs, errs = validate_sequence_file(fasta_path, mol_type="DNA")
        results.append(("Sequence validation (valid DNA)", "PASS" if not errs else "FAIL", str(errs) if errs else "ok"))
        _os.unlink(fasta_path)

        # 无效字母
        fd, fasta_path = _tempfile.mkstemp(suffix=".fasta")
        _os.close(fd)
        _os.chmod(fasta_path, 0o600)
        _Path(fasta_path).resolve().write_text(">seq1\nATGCGTXYZ!\n", encoding="utf-8")
        seqs, errs = validate_sequence_file(fasta_path, mol_type="DNA")
        results.append(("Sequence validation (invalid alphabet)", "PASS" if errs else "FAIL", "detected" if errs else "missed"))
        _os.unlink(fasta_path)

        # 重复 ID
        fd, fasta_path = _tempfile.mkstemp(suffix=".fasta")
        _os.close(fd)
        _os.chmod(fasta_path, 0o600)
        _Path(fasta_path).resolve().write_text(">seq1\nATGC\n>seq1\nATGC\n", encoding="utf-8")
        try:
            validate_sequence_file(fasta_path, mol_type="DNA")
            results.append(("Sequence validation (duplicate ID)", "FAIL", "not detected"))
        except Exception:
            results.append(("Sequence validation (duplicate ID)", "PASS", "detected"))
        _os.unlink(fasta_path)

        # 长度不一致
        fd, fasta_path = _tempfile.mkstemp(suffix=".fasta")
        _os.close(fd)
        _os.chmod(fasta_path, 0o600)
        _Path(fasta_path).resolve().write_text(">seq1\nATGCGT\n>seq2\nATGCGTACGT\n", encoding="utf-8")
        seqs, errs = validate_sequence_file(fasta_path, mol_type="DNA")
        results.append(("Sequence validation (length mismatch)", "PASS" if any("length" in e.lower() for e in errs) else "FAIL", str(errs) if errs else "missed"))
        _os.unlink(fasta_path)
    except Exception as e:
        results.append(("Sequence validation", "FAIL", str(e)))
    return results


def _test_control_chars() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.validation import check_control_characters, validate_string_safety
        from markerfinder.exceptions import SequenceValidationError

        # 检测双向覆盖字符
        bidi_str = "normal_text\u202Ereversed"
        errs = check_control_characters(bidi_str, "test_field")
        results.append(("Control char detection (bidi)", "PASS" if errs else "FAIL", "detected" if errs else "missed"))

        # 检测控制字符
        ctrl_str = "text\x00null\x01char"
        errs2 = check_control_characters(ctrl_str, "test_field")
        results.append(("Control char detection (NUL)", "PASS" if errs2 else "FAIL", "detected" if errs2 else "missed"))

        # 正常字符串不误报
        clean_str = "normal_label_no_issues"
        errs3 = check_control_characters(clean_str, "test_field")
        results.append(("Control char detection (clean)", "PASS" if not errs3 else "FAIL", "clean" if not errs3 else "false positive"))

        # Validate_string_safety 应抛异常
        try:
            validate_string_safety("\u202Eevil", "test")
            results.append(("String safety rejection", "FAIL", "not rejected"))
        except SequenceValidationError:
            results.append(("String safety rejection", "PASS", "rejected"))
    except Exception as e:
        results.append(("Control char detection", "FAIL", str(e)))
    return results


def _test_taxonomy_cycle() -> List[tuple]:
    results: List[tuple] = []
    try:
        import tempfile as _tempfile
        import os as _os
        from markerfinder.taxonomy import load_taxonomy_table
        from markerfinder.exceptions import TaxonomyConflictError

        fd, taxa_path = _tempfile.mkstemp(suffix=".tsv")
        _os.close(fd)
        _os.chmod(taxa_path, 0o600)
        from pathlib import Path as _Path
        _Path(taxa_path).resolve().write_text("A\td__A;p__B\nB\td__B;p__A\n", encoding="utf-8")
        try:
            load_taxonomy_table(taxa_path)
            results.append(("Taxonomy cycle detection", "FAIL", "not detected"))
        except TaxonomyConflictError:
            results.append(("Taxonomy cycle detection", "PASS", "detected"))
        _os.unlink(taxa_path)
    except Exception as e:
        results.append(("Taxonomy cycle detection", "FAIL", str(e)))
    return results


def _test_custom_levels() -> List[tuple]:
    results: List[tuple] = []
    try:
        from markerfinder.taxonomy import parse_taxonomy, get_merged_level_maps

        map_a, map_b, names = get_merged_level_maps("kingdom:_k_,subspecies:_ss_")
        # Map_a 的键是去掉下划线的单字母/多字母代码（'k'、'ss'），与 Format A/B 解析一致
        has_kingdom = "k" in map_a and "kingdom" in names
        has_subspecies = "ss" in map_a and "subspecies" in names
        results.append(("Custom taxonomy levels", "PASS" if has_kingdom and has_subspecies else "FAIL", f"kingdom={has_kingdom}, subspecies={has_subspecies}"))

        r = parse_taxonomy("ID_k_Eukarya_d_Animalia_p_Chordata", custom_levels="kingdom:_k_")
        results.append(("Custom level parsing (Format A)", "PASS" if r.get("kingdom") == "Eukarya" else "FAIL", str(r.get("kingdom"))))

        r2 = parse_taxonomy("k__Eukarya;d__Animalia;p__Chordata", custom_levels="kingdom:k")
        results.append(("Custom level parsing (Format B)", "PASS" if r2.get("kingdom") == "Eukarya" else "FAIL", str(r2.get("kingdom"))))
    except Exception as e:
        results.append(("Custom taxonomy levels", "FAIL", str(e)))
    return results


def _test_composition() -> List[tuple]:
    """The composition diagnostics get a self-test of their own.

    Not a smoke test. Each item can FAIL on its own: RCV has to actually
    discriminate a biased set from a uniform one (a metric that returns a
    constant is what the zero-discrimination probe measures on the reference tree),
    and the protein path has
    to answer "not applicable" rather than a numeric zero -- the GC-bias of an
    amino-acid alignment is undefined, and 0.0 there reads as "no bias detected".
    """
    results: List[tuple] = []
    try:
        from markerfinder.modules.composition import (
            gc_or_codon_bias,
            rcv,
        )

        # RCV with labels is the WORST between-group deviation from the pooled
        # Composition, so the fixture needs two groups: one label would make every
        # Set score exactly 0 by construction (my first version did that and the
        # Check correctly failed on its own bad control).
        labels = ["g1"] * 6 + ["g2"] * 6
        uniform = ["AC" * 12] * 6 + ["CA" * 12] * 6      # Same composition, 2 groups
        skewed = ["A" * 24] * 6 + ["C" * 24] * 6          # Opposite compositions
        flat = rcv(uniform, labels)
        biased = rcv(skewed, labels)
        if biased > flat + 1e-9:
            results.append((
                "Composition RCV discriminates", "PASS",
                f"skewed {biased:.4f} > uniform {flat:.4f}",
            ))
        else:
            results.append((
                "Composition RCV discriminates", "FAIL",
                f"skewed {biased:.4f} vs uniform {flat:.4f}: the metric is flat",
            ))
    except Exception as e:
        results.append(("Composition RCV discriminates", "FAIL", str(e)))

    try:
        bias = gc_or_codon_bias(["ACDEFGHIK", "KIHGFEDCA"], protein=True)
        if bias is None:
            results.append((
                "Composition: protein GC bias is NA", "PASS",
                "None, rendered NA (not 0.0)",
            ))
        else:
            results.append((
                "Composition: protein GC bias is NA", "FAIL",
                f"returned {bias!r}; a number here is a placeholder (G1)",
            ))
    except Exception as e:
        results.append(("Composition: protein GC bias is NA", "FAIL", str(e)))

    try:
        empty = rcv([], [])
        if empty == 0.0:
            results.append((
                "Composition: empty input is a defined zero", "PASS",
                "0.0 (no sites, not an unmeasurable marker)",
            ))
        else:
            results.append((
                "Composition: empty input is a defined zero", "FAIL",
                f"returned {empty!r}",
            ))
    except Exception as e:
        results.append(("Composition: empty input is a defined zero", "FAIL", str(e)))
    return results


def _test_assertions() -> List[tuple]:
    """Every registered assertion must fire on its negative
    fixture and pass on the healthy state. An assertion that cannot go red is
    treated as absent."""
    results: List[tuple] = []
    try:
        from markerfinder.assertions import (
            REGISTRY,
            _good_state,
            build_negative_state,
            run_assertions,
        )

        good = run_assertions(_good_state(), mode="check")
        bad_good = [
            r.assertion_id for r in good.results
            if not r.passed and r.severity == "fail"
        ]
        results.append((
            "Assertions (healthy state)",
            "PASS" if not bad_good else "FAIL",
            f"{good.pass_count}/{len(good.results)} pass" if not bad_good
            else f"falsely fired: {bad_good}",
        ))

        for spec in REGISTRY:
            report = run_assertions(build_negative_state(spec.id), mode="check")
            fired = [
                r for r in report.results
                if not r.passed and r.assertion_id == spec.id
            ]
            collateral = [
                r.assertion_id for r in report.results
                if not r.passed and r.assertion_id != spec.id
            ]
            if fired and not collateral:
                results.append((f"Assertion {spec.id} must-fail", "PASS", "fired"))
            elif fired:
                results.append((
                    f"Assertion {spec.id} must-fail", "FAIL",
                    f"fired but collateral: {collateral}",
                ))
            else:
                results.append((
                    f"Assertion {spec.id} must-fail", "FAIL",
                    "did NOT fire on its negative fixture",
                ))
    except Exception as e:
        results.append(("Assertions (registry)", "FAIL", str(e)))
    return results


def _test_database_hashes(db_dir: Optional[str] = None) -> List[tuple]:
    """A declared database hash must actually be compared.

    ``db/expected_hashes.json`` ships with empty values, so on a fresh install
    nothing is compared. That state must be reported as INFO with an explicit
    "not a pass" note: a comparison that never ran reading as green is the
    exact failure mode this refactor exists to remove.
    """
    from pathlib import Path as _P

    from markerfinder.utils.db_versioning import DatabaseVersionManager

    results: List[tuple] = []
    try:
        base = _P(db_dir) if db_dir else _P(__file__).resolve().parents[2] / "db"
        manager = DatabaseVersionManager(str(base))
        declared = {
            name: value.strip()
            for name, value in manager._load_expected_hashes().items()  # Noqa: SLF001
            if isinstance(value, str) and value.strip()
            and not name.startswith("_")
        }
        if not declared:
            results.append((
                "Database SHA-256 comparison", "INFO",
                "db/expected_hashes.json declares no populated entries — "
                "comparison NOT RUN (this is not a pass; fill the file to "
                "enable database-hash verification)",
            ))
            return results
        for name, expected in sorted(declared.items()):
            target = base / name
            if not target.exists():
                results.append((
                    f"Database hash [{name}]", "INFO",
                    f"expected hash declared but {target} is not present in "
                    f"this install — NOT CHECKED",
                ))
                continue
            info = manager.verify_database(name, str(target))
            if info.get("hash_match") is True:
                results.append((
                    f"Database hash [{name}]", "PASS",
                    f"SHA-256 matches the expected value ({expected[:12]}…)",
                ))
            elif info.get("hash_match") is False:
                results.append((
                    f"Database hash [{name}]", "FAIL",
                    f"SHA-256 {str(info.get('hash'))[:12]}… != expected "
                    f"{expected[:12]}… — the database changed since it was "
                    f"pinned since it was recorded",
                ))
            else:
                results.append((
                    f"Database hash [{name}]", "FAIL",
                    "comparison produced no verdict (hash_match is None despite"
                    " a declared expected value)",
                ))
    except Exception as e:  # Noqa: BLE001 - self-test never crashes the CLI
        results.append(("Database SHA-256 comparison", "FAIL", str(e)))
    return results


def _run_self_test(db_dir: Optional[str] = None) -> int:
    """Run self-test: verify dependencies and built-in examples.

    ``db_dir`` is threaded through to the hash check: ``--check`` used to
    call ``_test_database_hashes`` with no argument, so ``--check --db-dir X``
    verified the repository's default ``db/`` and reported it as if it had
    examined the directory the user named. A check that looks at the wrong
    artifact is worse than no check, because it prints a verdict.

    测试场景矩阵覆盖：
      - 依赖检查
      - 格式 A 解析（正常/缺失级别/含下划线分类名）
      - 格式 B 解析（正常/分隔符冲突/缺失值）
      - 外部表格补充与混合优先级
      - 单系群判定
      - Newick 验证（正常/负分支/重复tip/空节点）
      - 序列验证（字母表异常/ID重复/长度不一致）
      - 恶意字符注入拒绝
      - 分类学循环依赖
      - 多棵树处理
    """
    print_banner()
    results: List[tuple] = []

    # --- 1. 依赖检查 ---
    results += _test_dependencies()
    # --- 1b. 外部工具面披露：INFO 项，不参与计分 ---
    results += _test_external_tools()

    # --- 2. 格式 A 解析 ---
    results += _test_format_a()

    # --- 3. 格式 B 解析 ---
    results += _test_format_b()

    # --- 4. 混合优先级 ---
    results += _test_merge_taxonomy()

    # --- 5. 单系群判定 ---
    results += _test_monophyly()
    results += _test_database_hashes(db_dir)

    # --- 6. 树文件验证 ---
    results += _test_newick_validation()

    # --- 7. 序列验证 ---
    results += _test_sequence_validation()

    # --- 8. 恶意字符注入拒绝 ---
    results += _test_control_chars()

    # --- 9. 分类学循环依赖 ---
    results += _test_taxonomy_cycle()

    # --- 10. 自定义级别扩展 ---
    results += _test_custom_levels()

    # --- 11. 输出断言 must-fail 控制 ---
    results += _test_composition()
    results += _test_assertions()

    # --- 12. 接线状态披露 ---
    # Phase 1.5 ortholog resolution is not wired into the pipeline; disclose
    # It here rather than letting the "adaptive" narrative overstate it.
    results.append((
        "Phase 1.5 ortholog resolution", "INFO",
        "not wired in — the module exists but the pipeline skips it",
    ))

    # --- 输出结果 ---
    max_name_len = max(len(r[0]) for r in results) if results else 20
    print(f"\n{'Test':<{max_name_len + 2}} {'Status':<8} {'Detail'}")
    print("-" * (max_name_len + 30))
    all_pass = True
    for name, status, detail in results:
        print(f"{name:<{max_name_len + 2}} [{status}]   {detail}")
        if status == "FAIL":
            all_pass = False

    print()
    return EXIT_SUCCESS if all_pass else EXIT_RUNTIME_ERROR
