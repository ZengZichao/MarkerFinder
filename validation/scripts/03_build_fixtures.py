#!/usr/bin/env python3
"""Step 3: build every remaining input fixture the validation cases read.

Everything here is derived from ``data/PROVENANCE.tsv`` (which is itself read
from NCBI), so no lineage, tree topology or quality number in the bundle is
typed by hand. The generated artefacts are:

``taxonomy/`` Format B tables (with and without a header line), Format A
                   (embedded-segment) table, a malformed-row table.
``trees/`` reference species tree per set (topology implied by the
                   NCBI lineages), an NHX-annotated variant, a multi-tree file,
                   and the illegal references the CLI must refuse.
``checkm/`` pre-computed CheckM table for ``--checkm-results``, using the
                   completeness/contamination NCBI itself records for each
                   assembly (``checkm_info`` in the datasets report).
``cog/`` marker -> functional-category map, where the category text is
                   the DESCRIPTION line of the bundled HMM profile.
``mustpass/`` taxonomy must-pass baseline for ``--taxonomy-mustpass``.
``configs/`` YAML / TOML / JSON config files for ``--config``.
``sequences/`` protein FASTA for the tree-vs-sequence cross-validation.
``MANIFEST.sha256`` checksum of every file under ``data/`` (provenance audit).

Usage
-----
    python 03_build_fixtures.py [--force]
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent
DATA = ROOT / "data"
PROVENANCE = DATA / "PROVENANCE.tsv"

REPO_DB = ROOT.parent / "db" / "gtdb_markers" / "bac120"

TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_CODE = {
    "domain": "d", "phylum": "p", "class": "c", "order": "o",
    "family": "f", "genus": "g", "species": "s",
}
SETS_FILE = DATA / "sets"


def read_provenance() -> List[Dict[str, str]]:
    if not PROVENANCE.exists():
        raise SystemExit(f"{PROVENANCE} missing — run 01_fetch_genomes.py first.")
    lines = PROVENANCE.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def read_sets() -> Dict[str, List[str]]:
    return {
        p.stem: [x.strip() for x in p.read_text(encoding="utf-8").splitlines() if x.strip()]
        for p in sorted(SETS_FILE.glob("*.txt"))
    }


def read_fasta(path: Path) -> Dict[str, str]:
    """Minimal FASTA reader: header id (first token) -> sequence."""
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif line.strip() and name is not None:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def lineage_b(row: Dict[str, str], upto: str = "genus") -> str:
    """Format B string ``d__X;p__Y;...`` up to and including ``upto``."""
    parts = []
    for rank in TAX_RANKS:
        code = RANK_CODE[rank]
        value = row.get(rank, "")
        parts.append(f"{code}__{value}")
        if rank == upto:
            break
    return ";".join(parts)


def lineage_a(row: Dict[str, str]) -> str:
    """Format A embedded segment ``_d_X_p_Y_..._g_Z`` (spaces -> underscores)."""
    parts = []
    for rank in TAX_RANKS[:-1]:
        value = (row.get(rank) or "NA").replace(" ", "_")
        parts.append(f"_{RANK_CODE[rank]}_{value}")
    return "".join(parts).lstrip("_")


# --------------------------------------------------------------------------
# Taxonomy tables
# --------------------------------------------------------------------------
def build_taxonomy_tables(rows: List[Dict[str, str]], sets: Dict[str, List[str]]) -> None:
    out = DATA / "taxonomy"
    out.mkdir(exist_ok=True)
    by_acc = {r["accession"]: r for r in rows}

    for name, members in sets.items():
        if name == "single1":
            continue
        lines = ["genome_id\ttaxonomy"]
        for acc in members:
            lines.append(f"{acc}\t{lineage_b(by_acc[acc])}")
        path = out / f"taxonomy_{name}.tsv"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"[write] {path} ({len(members)} rows)")

    members = sets["small8"]
    # Headerless twin: the reader must treat both shapes identically.
    body = (out / f"taxonomy_{members and 'small8'}.tsv").read_text(
        encoding="utf-8").splitlines()[1:]
    (out / "taxonomy_small8_noheader.tsv").write_text(
        "\n".join(body) + "\n", encoding="utf-8")
    print(f"[write] {out / 'taxonomy_small8_noheader.tsv'}")

    lines = ["genome_id\ttaxonomy_segment"]
    for acc in members:
        lines.append(f"{acc}\t{lineage_a(by_acc[acc])}")
    (out / "taxonomy_small8_formatA.tsv").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    print(f"[write] {out / 'taxonomy_small8_formatA.tsv'}")

    # Malformed rows: one row with no level prefix at all, one with a stray
    # Unknown prefix, and a truncated row. All must be reported, never ignored.
    bad = [
        "genome_id\ttaxonomy",
        members[0] + "\t" + lineage_b(by_acc[members[0]]),
        "NOT_A_TAXONOMY_ROW\tthis is prose, not a lineage",
        members[1] + "\tX__magic;d__Bacteria;p__Bacillota",
        members[2],
    ]
    (out / "taxonomy_malformed.tsv").write_text("\n".join(bad) + "\n", encoding="utf-8")
    print(f"[write] {out / 'taxonomy_malformed.tsv'}")

    # Comma-separated twin, to exercise --table-sep.
    csv_lines = ["genome_id,taxonomy"]
    for acc in members:
        csv_lines.append(f"{acc},{lineage_b(by_acc[acc])}")
    (out / "taxonomy_small8.csv").write_text("\n".join(csv_lines) + "\n",
                                             encoding="utf-8")
    print(f"[write] {out / 'taxonomy_small8.csv'}")


# --------------------------------------------------------------------------
# Reference trees
# --------------------------------------------------------------------------
def _build_newick(tips: Dict[str, str], ranks: List[str]) -> str:
    """Nested Newick whose topology is implied by the lineage strings.

    ``tips`` maps tip label -> ``;``-separated lineage (Format B). Internal
    nodes are unnamed; every edge gets one constant length so that no tool
    downstream sees a zero-length branch and calls it a polytomy.

    A rank that does not split its tips adds no node and therefore no edge:
    earlier code appended a length at *every* rank a leaf passed through, which
    produced ``GCF_000009045.1:0.01:0.01:0.01`` — parseable by nobody, and read
    by ete3 as a node literally named ``:0.01``.

    Where one taxonomic level holds more than two children, the remaining order
    is not implied by the lineage strings, and it is resolved into a balanced
    binary subtree in sorted-label order rather than left as a polytomy: the
    reference validator refuses internal nodes of degree > 3 (POLYTOMY), because
    an unresolved reference cannot measure incongruence. The resolution is
    therefore arbitrary but deterministic, and it never crosses a taxonomic
    boundary — every clade in the output is a taxon.
    """
    bl = 0.01
    groups: Dict[Tuple[str, ...], List[str]] = {}
    for tip, lineage in tips.items():
        levels = {}
        for pair in lineage.split(";"):
            code, _, value = pair.partition("__")
            levels[code] = value
        key = tuple(levels.get(RANK_CODE[r], "") for r in ranks)
        groups.setdefault(key, []).append(tip)

    def resolve(parts: List[str]) -> str:
        """Join child expressions into one binary subtree."""
        if len(parts) == 1:
            return parts[0]
        if len(parts) == 2:
            return f"({parts[0]},{parts[1]}):{bl}"
        mid = len(parts) // 2
        return f"({resolve(parts[:mid])},{resolve(parts[mid:])}):{bl}"

    def nest(items: List[Tuple[Tuple[str, ...], List[str]]],
             depth: int) -> str:
        """The subtree spanned by ``items`` at rank ``depth``."""
        if depth >= len(ranks):
            return resolve([f"{tip}:{bl}" for tip in
                            sorted(t for _k, ts in items for t in ts)])
        by_value: Dict[str, List[Tuple[Tuple[str, ...], List[str]]]] = {}
        for key, ts in items:
            by_value.setdefault(key[depth], []).append((key, ts))
        children = [nest(sub, depth + 1) for _v, sub in sorted(by_value.items())]
        return resolve(children)

    root = nest(sorted(groups.items()), 0)
    if not root.startswith("("):
        root = f"({root})"
    return root + ";"


_DOUBLE_LENGTH = re.compile(r":\d+(?:\.\d+)?\s*:")


def _node_degrees(newick: str) -> List[int]:
    """Child counts of every internal node, by scanning the Newick tokens."""
    degrees: List[int] = []
    stack: List[int] = []
    for ch in newick:
        if ch == "(":
            stack.append(1)
        elif ch == ",":
            if stack:
                stack[-1] += 1
        elif ch == ")":
            if stack:
                degrees.append(stack.pop())
    return degrees


def _validate_reference(newick: str, members: List[str], path: Path) -> None:
    """Refuse to ship a reference tree that is not a resolved tree over these genomes.

    Every consumer of these files reads tips as genome ids, so a malformed
    Newick does not merely look wrong — it silently changes which labels exist.
    And an unresolved node makes the reference unusable for incongruence
    measurement (the validator reports POLYTOMY), so that is caught here rather
    than discovered by a case that then has nothing to compare against.
    """
    if _DOUBLE_LENGTH.search(newick):
        raise SystemExit(
            f"{path}: generated tree repeats a branch length "
            f"({_DOUBLE_LENGTH.search(newick).group(0)!r}); the writer is broken"
        )
    found = set(re.findall(r"([A-Za-z0-9_.\-]+):(?:\d+(?:\.\d+)?)", newick))
    missing = set(members) - found
    extra = found - set(members)
    if missing or extra:
        raise SystemExit(
            f"{path}: generated tree does not cover the set "
            f"(missing={sorted(missing)} unexpected={sorted(extra)})"
        )
    if newick.count(";") != 1:
        raise SystemExit(f"{path}: not exactly one Newick statement")
    polytomies = [d for d in _node_degrees(newick) if d > 2]
    if len(members) > 1 and polytomies:
        raise SystemExit(
            f"{path}: generated tree leaves {len(polytomies)} unresolved node(s) "
            f"(max degree {max(polytomies)}); the reference validator rejects "
            "polytomies, so the tree could not measure anything"
        )


def build_trees(rows: List[Dict[str, str]], sets: Dict[str, List[str]]) -> None:
    out = DATA / "trees"
    out.mkdir(exist_ok=True)
    by_acc = {r["accession"]: r for r in rows}

    for name, members in sets.items():
        tips = {acc: lineage_b(by_acc[acc]) for acc in members}
        newick = _build_newick(tips, ["domain", "phylum", "class", "order",
                                      "family", "genus"])
        path = out / f"reference_{name}.nwk"
        _validate_reference(newick, members, path)
        path.write_text(newick + "\n", encoding="utf-8")
        print(f"[write] {path}")
        if name != "small8":
            continue
        # NHX-annotated twin: same topology, per-tip annotations to strip.
        nhx = newick
        for acc in members:
            genus = by_acc[acc]["genus"].replace(" ", "_")
            comp = by_acc[acc]["checkm_completeness"]
            nhx = nhx.replace(
                f"{acc}:0.01",
                f"{acc}[&&NHX:taxon={genus}&completeness={comp}]:0.01",
            )
        (out / "reference_small8.nhx.nwk").write_text(nhx + "\n", encoding="utf-8")
        print(f"[write] {out / 'reference_small8.nhx.nwk'}")

        multi = "\n".join([newick, newick.replace(":0.01", ":0.02"), newick])
        (out / "reference_small8_multi.nwk").write_text(multi + "\n",
                                                        encoding="utf-8")
        print(f"[write] {out / 'reference_small8_multi.nwk'}")

    # Illegal references, generated (not copied from tests/) so the bundle is
    # Self-contained: too few tips, duplicated tips, negative branch length.
    (out / "illegal_three_tips.nwk").write_text(
        "((GCF_000009045.1:0.1,GCF_000007825.1:0.1):0.1,GCF_000009645.1:0.1);\n",
        encoding="utf-8")
    (out / "illegal_duplicate_tips.nwk").write_text(
        "((A:0.1,B:0.1)N1:0.1,(A:0.1,C:0.1)N2:0.1)ROOT;\n", encoding="utf-8")
    (out / "illegal_negative_branch.nwk").write_text(
        "(((GCF_000009045.1:0.1,GCF_000007825.1:0.1):0.1,"
        "GCF_000009645.1:0.1):-0.5,GCF_000196035.1:0.1);\n", encoding="utf-8")
    for f in ("illegal_three_tips.nwk", "illegal_duplicate_tips.nwk",
              "illegal_negative_branch.nwk"):
        print(f"[write] {out / f}")


# --------------------------------------------------------------------------
# CheckM / COG / must-pass / config / sequence fixtures
# --------------------------------------------------------------------------
def build_checkm(rows: List[Dict[str, str]], sets: Dict[str, List[str]]) -> None:
    out = DATA / "checkm"
    out.mkdir(exist_ok=True)
    by_acc = {r["accession"]: r for r in rows}
    header = "Bin Id\tMarker lineages\tCompleteness\tContamination\tQuality"
    for name, members in sets.items():
        lines = [header]
        skipped = []
        for acc in members:
            r = by_acc[acc]
            raw_comp, raw_cont = r["checkm_completeness"], r["checkm_contamination"]
            if raw_comp in ("", "None") or raw_cont in ("", "None"):
                # NCBI records no checkm_info for this assembly. Nothing is
                # Invented: the row is omitted and the reason is written into
                # The file, so the case can assert the documented fallback.
                skipped.append(acc)
                continue
            comp = float(raw_comp)
            cont = float(raw_cont)
            quality = round(comp - 5 * cont, 2)
            lines.append(
                f"{acc}\t[unassigned]\t{comp:.2f}\t{cont:.2f}\t{quality:.2f}"
            )
        for acc in skipped:
            lines.insert(1, f"# {acc}: NCBI reports no checkm_info; row omitted")
        path = out / f"checkm_results_{name}.tsv"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"[write] {path} "
              f"({len(members) - len(skipped)} rows, {len(skipped)} omitted)")


def hmm_description(hmm_path: Path) -> str:
    """The functional-category text carried by an HMM profile.

    HMMER profiles from Pfam/TIGRFAM put the description in ``DESC``; some
    GTDB-converted files use ``DESCRIPTION`` or a ``# DESCRIPTION`` comment
    instead. The first non-empty one is the text, and it is the profile's own
    annotation — nothing here is asserted about COG categories, because these
    profiles do not carry COG assignments.
    """
    wanted = ("DESC", "DESCRIPTION", "# DESCRIPTION", "NAME")
    seen: Dict[str, str] = {}
    patterns = [(key, re.compile(rf"^{re.escape(key)}\s+(.*)$")) for key in wanted]
    for line in hmm_path.read_text(encoding="utf-8", errors="replace").splitlines():
        for key, pattern in patterns:
            if key in seen:
                continue
            match = pattern.match(line)
            if match and match.group(1).strip():
                seen[key] = match.group(1).strip()
        if len(seen) == len(wanted):
            break
    for key in wanted:
        if key in seen:
            return seen[key]
    return ""


def build_cog_map() -> None:
    out = DATA / "cog"
    out.mkdir(exist_ok=True)
    hmms = sorted((REPO_DB).glob("*.HMM")) + sorted(REPO_DB.glob("*.hmm"))
    core_dir = DATA / "markers" / "small8_core"
    if not core_dir.exists():
        print("[skip] cog map: run 02 first (needs the core marker subset)")
        return
    core_ids = sorted(p.stem for p in core_dir.glob("*.faa"))
    by_stem = {p.stem: p for p in hmms}
    lines = ["marker_id\tfunctional_category"]
    for marker_id in core_ids:
        desc = hmm_description(by_stem[marker_id]) if marker_id in by_stem else ""
        lines.append(f"{marker_id}\t{desc or 'NA'}")
    path = out / "marker_category_map.tsv"
    filled = sum(1 for line in lines[1:] if not line.endswith("\tNA"))
    if filled == 0:
        raise SystemExit(
            f"{path}: no profile yielded a description line — the extractor is "
            "reading a field these profiles do not use, and the map would be "
            "unable to demonstrate --cog-category-map at all"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[write] {path} ({filled}/{len(lines) - 1} markers annotated; text = "
          f"the profile's own DESC/DESCRIPTION line)")


def build_mustpass(sets: Dict[str, List[str]], rows: List[Dict[str, str]]) -> None:
    """Taxonomy must-pass baselines plus its twin table.

    Two fixtures, one in each direction, because a gate that is only ever seen
    passing has not been shown to be able to fail:

    ``mustpass_domain.yaml`` — monophyly at ``domain``, min_taxa 2. Every
    shipped genome is Bacteria, so this is satisfiable and the case asserts the
    gate RAN (a non-zero number of groups tested) and passed.

    ``mustpass_splitgenus.yaml`` + ``taxonomy_small8_splitgenus.tsv`` — the
    same requirement at ``genus``, against a table in which three accessions
    from different orders and one different phylum are deliberately given one
    fake genus. Their gene trees cannot make that group exclusive of the other
    five tips, so the gate must abort the run with exit code 4.
    """
    out = DATA / "mustpass"
    out.mkdir(exist_ok=True)
    (out / "mustpass_domain.yaml").write_text(
        "# Taxonomy must-pass baseline. The requirement it encodes is\n"
        "# domain, where the shipped set has >= 2 representatives and every\n"
        "# member shares one value, so the gate can be shown to run and pass.\n"
        "# The failure direction is mustpass_splitgenus.yaml.\n"
        "markers:\n"
        "  note: \"applies to every gene tree the run produces\"\n"
        "  ids: []\n"
        "must_pass:\n"
        "  - relation: \"monophyletic\"\n"
        "    level: \"domain\"\n"
        "    min_taxa: 2\n"
        "    note: \"all shipped genomes are Bacteria\"\n",
        encoding="utf-8")
    (out / "mustpass_splitgenus.yaml").write_text(
        "# Failure control: monophyly at genus against\n"
        "# taxonomy/taxonomy_small8_splitgenus.tsv, which assigns three\n"
        "# accessions from different orders (and one from another phylum) to a\n"
        "# single fake genus. No marker tree can make that group exclusive of\n"
        "# the remaining tips, so the gate must abort with exit code 4.\n"
        "markers:\n"
        "  ids: []\n"
        "must_pass:\n"
        "  - relation: \"monophyletic\"\n"
        "    level: \"genus\"\n"
        "    min_taxa: 2\n"
        "    note: \"g__Splitgenus spans Bacillota orders and a Pseudomonadota\"\n",
        encoding="utf-8")
    print(f"[write] {out / 'mustpass_domain.yaml'}, mustpass_splitgenus.yaml")

    # The deliberately mislabelled twin table, written next to the honest ones.
    by_acc = {r["accession"]: r for r in rows}
    split_members = {
        "GCF_000009045.1",   # Bacillus subtilis (Caryophanales)
        "GCF_000006785.1",   # Streptococcus pyogenes (Lactobacillales)
        "GCF_000008105.1",   # Salmonella enterica (Pseudomonadota)
    }
    lines = ["genome_id\ttaxonomy"]
    for acc in sets["small8"]:
        row = dict(by_acc[acc])
        if acc in split_members:
            row["genus"] = "Splitgenus"
        lines.append(f"{acc}\t{lineage_b(row)}")
    path = DATA / "taxonomy" / "taxonomy_small8_splitgenus.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[write] {path} (3 tips share the fake genus Splitgenus)")


def build_configs(sets: Dict[str, List[str]]) -> None:
    """Config files for ``--config``, in all three accepted formats.

    Deliberately no ``input:`` key: the materialised input directory is a
    per-worker runtime path, so a config that pinned it would either point at a
    directory that does not exist or force the suite to write into one fixed
    place. The cases pass ``-i``/``-o``; the file supplies everything else.

    The set is ``quad4`` — the fast default — so a config-file run costs the
    same as any other feature-level case.
    """
    out = DATA / "configs"
    out.mkdir(exist_ok=True)
    yaml_lines = [
        "# YAML config for the validation quad4 cases. Keys are CLI option",
        "# names with dashes written as underscores; paths are relative to the",
        "# repository root, which is where the suite starts every run.",
        "threads: 2",
        "mode: standard",
        "marker_mode: gtdb_tk",
        "gtdb_markers_dir: validation/data/markers/quad4_core",
        "taxonomy_table: validation/data/taxonomy/taxonomy_quad4.tsv",
        "monophyly_rank: order",
        "max_markers: 4",
        "gene_tree_builder: fasttree",
        'coalescent_mode: "off"',
        "skip_checkm: true",
        "verbose: 1",
    ]
    yaml_text = "\n".join(yaml_lines) + "\n"
    (out / "quad4.yaml").write_text(yaml_text, encoding="utf-8")
    # TOML twin: same keys, TOML quoting.
    body = []
    for line in yaml_lines:
        if line.startswith("#"):
            continue
        key, _, value = line.partition(":")
        value = value.strip().strip('"')
        if value in ("true", "false"):
            body.append(f"{key.strip()} = {value}")
        elif value.isdigit():
            body.append(f"{key.strip()} = {value}")
        else:
            body.append(f'{key.strip()} = "{value}"')
    (out / "quad4.toml").write_text(
        "# TOML twin of quad4.yaml\n" + "\n".join(body) + "\n", encoding="utf-8")

    # JSON twin shaped like a run_config.json written by a previous run: its
    # Top-level "parameters" mapping must be flattened by --config.
    params = {}
    for line in body:
        key, _, value = line.partition(" = ")
        params[key] = value
    json_lines = ['{', '  "parameters": {']
    items = list(params.items())
    for i, (key, value) in enumerate(items):
        comma = "," if i < len(items) - 1 else ""
        json_lines.append(f'    "{key}": {value}{comma}')
    json_lines += ['  }', '}']
    (out / "quad4_flat_params.json").write_text(
        "\n".join(json_lines) + "\n", encoding="utf-8")

    for f in ("quad4.yaml", "quad4.toml", "quad4_flat_params.json"):
        print(f"[write] {out / f}")
    assert sets["quad4"], "set list is empty"


def build_sequence_fixtures(rows: List[Dict[str, str]],
                            sets: Dict[str, List[str]]) -> None:
    """Protein and nucleotide side files for --sequences / --mol-type paths."""
    out = DATA / "sequences"
    out.mkdir(exist_ok=True)
    by_acc = {r["accession"]: r for r in rows}

    # One protein sequence per tip, taken from that genome's own core markers.
    marker_dir = DATA / "markers" / "small8_core"
    if not marker_dir.exists():
        print("[skip] sequence fixtures: run 02 first")
        return
    picked: Dict[str, str] = {}
    for faa in sorted(marker_dir.glob("*.faa")):
        for tip, seq in read_fasta(faa).items():
            picked.setdefault(tip, seq)
    protein_path = out / "tip_sequences.faa"
    lines = []
    for tip, seq in sorted(picked.items()):
        lines.append(f">{tip}\n{seq[:120]}")
    protein_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[write] {protein_path} ({len(lines) // 2} tips)")

    # A short nucleotide FASTA per tip: only used to check the.fna sidecar
    # Handling and the --mol-type DNA alphabet validation.
    codon_table = "TCAG"
    nuc_lines = []
    for acc in sets["small8"]:
        dna = "".join(codon_table[(i * 7 + len(acc)) % 4] for i in range(300))
        nuc_lines.append(f">{acc}\n{dna}")
    (out / "tips_partial.fna").write_text("\n".join(nuc_lines) + "\n",
                                          encoding="utf-8")
    print(f"[write] {out / 'tips_partial.fna'} "
          f"(synthetic 300 nt per tip: alphabet checks only, no biology)")
    assert by_acc or True


def build_mag_named_variant(sets: Dict[str, List[str]], rows: List[Dict[str, str]]) -> None:
    """An identifier-scheme variant of quad4, for the MAG pathway.

    MarkerFinder classifies a sample by FILE NAME (``mag``/``bin`` -> MAG,
    ``sag`` -> SAG), and the sequences themselves are irrelevant to that rule.
    So the variant renames the four quad4 genomes to ``mag_01..mag_04`` and
    rewrites the marker FASTA headers and the taxonomy table to match. No
    sequence is altered: this is a labelling experiment, and it is recorded as
    one.
    """
    out = DATA / "variants" / "mag_named"
    genomes = out / "genomes"
    markers = out / "markers"
    taxonomy = out / "taxonomy"
    for d in (genomes, markers, taxonomy):
        d.mkdir(parents=True, exist_ok=True)
    members = sets["quad4"]
    rename = {acc: f"mag_{i + 1:02d}" for i, acc in enumerate(members)}

    for acc, new in rename.items():
        target = genomes / f"{new}.faa"
        if target.is_symlink() or target.exists():
            target.unlink()
        target.symlink_to(Path("..") / ".." / ".." / "genomes" / f"{acc}.faa")

    src_markers = DATA / "markers" / "quad4_core"
    for faa in sorted(src_markers.glob("*.faa")):
        seqs = read_fasta(faa)
        lines = [f">{rename[acc]}\n{seqs[acc]}" for acc in members if acc in seqs]
        (markers / faa.name).write_text("\n".join(lines) + "\n", encoding="utf-8")

    by_acc = {r["accession"]: r for r in rows}
    lines = ["genome_id\ttaxonomy"]
    for acc, new in rename.items():
        lines.append(f"{new}\t{lineage_b(by_acc[acc])}")
    (taxonomy / "taxonomy_mag_named.tsv").write_text("\n".join(lines) + "\n",
                                                    encoding="utf-8")
    (out / "README.txt").write_text(
        "Identifier-scheme variant of the quad4 set: the same four genomes and "
        "the same sequences, renamed mag_01..mag_04 so that MarkerFinder's "
        "filename-based MAG classification fires. Generated by "
        "scripts/03_build_fixtures.py; nothing here is a different measurement.\n"
        f"rename map: {rename}\n", encoding="utf-8")
    print(f"[write] {out} (4 MAG-named genomes, "
          f"{len(list(markers.glob('*.faa')))} marker files)")


# --------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------
def write_manifest() -> None:
    root = DATA
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root).as_posix()
        if rel.startswith((".download_tmp/", "markers/.hit_cache")):
            continue
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 16), b""):
                h.update(chunk)
        rows.append(f"{h.hexdigest()}  {rel}")
    path = root / "MANIFEST.sha256"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"[write] {path} ({len(rows)} files)")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true", help="(re)build every fixture")
    ap.add_argument("--skip-manifest", action="store_true")
    args = ap.parse_args(argv)

    if shutil.which("datasets") is None and not args.force:
        pass  # Metadata fetching is step 1's job; this step needs no network

    rows = read_provenance()
    sets = read_sets()
    build_taxonomy_tables(rows, sets)
    build_trees(rows, sets)
    build_checkm(rows, sets)
    build_cog_map()
    build_mustpass(sets, rows)
    build_configs(sets)
    build_sequence_fixtures(rows, sets)
    build_mag_named_variant(sets, rows)
    if not args.skip_manifest:
        write_manifest()
    return 0


if __name__ == "__main__":
    sys.exit(main())
