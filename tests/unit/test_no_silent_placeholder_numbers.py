"""沉默即失败 — an ``except`` must not return a fabricated number.

"任何 except 分支返回数值的路径必须要么改为 None，要么显式声明为合法的业务值并注释说明."

Without this test the rule is never checked at all; scanning the
package found four ``except`` blocks returning numeric literals, none of them
annotated, and one of them (``_gene_tree_mean_support``) returned exactly the
0.5 "neutral placeholder" that goal G1 and forbid.

The scan is deliberately narrow: an ``except`` whose first real statement
returns a numeric literal. A site may be exempted only by carrying an explicit
justification marker (合法的业务值) near the return, so exemptions are
reviewable rather than invisible.
"""

from __future__ import annotations

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parents[2]
PACKAGE = REPO / "markerfinder"

RETURN_LITERAL = re.compile(
    r"return\s+(-?\d+(?:\.\d+)?|float\(['\"]?nan['\"]?\))\b"
)
EXCEPT_LINE = re.compile(r"^\s*except\b")
JUSTIFICATION = ("合法的业务值", "legitimate business value",
                 "not a placeholder")


def _scan_directory(root: pathlib.Path):
    """Run the same rule over an arbitrary tree of.py files."""
    offenders, total = [], 0
    for path in sorted(root.rglob("*.py")):
        lines = path.read_text(encoding="utf-8", errors="replace").split("\n")
        for index, line in enumerate(lines):
            if not EXCEPT_LINE.match(line):
                continue
            indent = len(line) - len(line.lstrip())
            for offset in range(index + 1, min(index + 10, len(lines))):
                body = lines[offset]
                if not body.strip():
                    continue
                if len(body) - len(body.lstrip()) <= indent:
                    break
                match = RETURN_LITERAL.search(body)
                if match:
                    total += 1
                    window = "\n".join(lines[max(0, index - 6):offset + 1])
                    if not any(m in window for m in JUSTIFICATION):
                        offenders.append(
                            f"{path.relative_to(root)}:{index + 1} "
                            f"{match.group(0)}"
                        )
                    break
    return total, offenders


def _scan():
    return _scan_directory(PACKAGE)


def test_scanner_detects_a_known_bad_site(tmp_path):
    """Positive control, self-contained.

    Counting how many sites the real package has would drift every time one is
    fixed (it did), so the control instead plants a known violation and a known
    justified one and requires the scanner to tell them apart. Without this, a
    regex that matched nothing would make "zero offenders" pass vacuously.
    """
    (tmp_path / "bad.py").write_text(
        "def f(x):\n"
        "    try:\n"
        "        return x.calc()\n"
        "    except Exception:\n"
        "        return 0.5\n",
        encoding="utf-8",
    )
    (tmp_path / "good.py").write_text(
        "def f(x):\n"
        "    try:\n"
        "        return x.count()\n"
        "    except OSError:\n"
        "        # 合法的业务值 — this is a count of zero items.\n"
        "        return 0\n",
        encoding="utf-8",
    )
    total, offenders = _scan_directory(tmp_path)
    assert total == 2, f"scanner matched {total} sites; it is not seeing code"
    assert any("bad.py" in o for o in offenders), offenders
    assert not any("good.py" in o for o in offenders), offenders


def test_no_unjustified_except_returns_a_number():
    total, offenders = _scan()
    assert not offenders, (
        "these except blocks turn a failure into a number "
        "without declaring it a legitimate business value:\n  "
        + "\n  ".join(offenders)
    )


def test_the_fixed_sites_now_declare_unmeasurable():
    """The four previously-silent sites behave as the doc requires."""
    from markerfinder.models.tree import Tree
    from markerfinder.pipeline import _gene_tree_mean_support

    # G1 /: no neutral 0.5 placeholder on the informativeness path.
    assert _gene_tree_mean_support(Tree(newick="((A,B),(C,D));")) is None
    # A reserved tree API asked about taxa that are not on the tree: the MRCA
    # Cannot be computed, which must not read as depth 0 (or any depth at all).
    assert Tree(newick="((A,B),(C,D));").get_mrca_depth({"A", "NOT_IN_TREE"}) is None
