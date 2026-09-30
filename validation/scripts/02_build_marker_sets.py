#!/usr/bin/env python3
"""Step 2: build the per-marker inputs the pipeline consumes.

Produced artefacts (all under ``validation/data/``)
--------------------------------------------------
``markers/<set>/<marker_id>.faa``
    One file per marker, one sequence per genome: the layout GTDB-Tk writes
    into ``<gtdbtk_out>/align/marker_genes/*.faa`` and which MarkerFinder's
    ``--marker-mode gtdb_tk`` consumes.

    How they are obtained, stated plainly: the sequences are the best HMM hit
    of each genome's own proteome against the TIGRFAM/Pfam profiles bundled in
    this repository (``db/gtdb_markers/bac120``). They are a **reconstruction
    of the GTDB-Tk file layout**, not a GTDB-Tk release output -- GTDB-Tk needs
    a ~100 GB reference-tree package, which cannot ship inside a test bundle.
    What the reconstruction preserves is everything the pipeline reads: the
    one-file-per-marker layout, the filename stem as marker id, unaligned
    protein sequences, and per-marker occupancy that spans the range a real
    dataset spans. What it does not reproduce is GTDB-Tk's own taxonomy-aware
    best-hit selection.

``markers/OCCUPANCY_ALL.tsv``
    Occupancy of every profile that was scanned. This is the audit trail for
    which markers entered each directory and why.

``markers/<set>_core/``
    The eight complete (occupancy = 1.0) markers. Feature-level cases run on
    this subset: the pipeline prepares every marker in the directory (align,
    trim, gene tree, PIS and composition sidecars), so wall-clock time is
    almost proportional to the number of marker files -- measured: 7m19s for
    one 4-genome run over 30 markers, 1m34s over 6.

``markers/<set>/``
    Twenty markers chosen to SPAN the occupancy range, so that an occupancy
    threshold can be shown to reject as well as to accept.

``markers/degenerate/``
    The lowest-occupancy markers (2 and 3 tips). Robustness cases only.

``hmms/full/``, ``hmms/core/``
    The matching profile directories for ``--marker-mode hmm``.

``markers/chimera_bench12/`` + ``markers/chimera_ground_truth_bench12.json``
    Positive control: in three markers the *Bacillus anthracis* Ames recipient
    sequence is replaced by the *Lactococcus lactis* donor sequence -- a
    cross-order transfer of a position-known gene, so detection sensitivity can
    be measured rather than asserted. Every other sequence in the directory is
    byte-identical to the clean set, so the clean run is the negative control.

Usage
-----
    python 02_build_marker_sets.py # build what is missing
    python 02_build_marker_sets.py --force # rebuild everything
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent
DATA = ROOT / "data"
GENOMES = DATA / "genomes"
SETS_FILE = DATA / "sets"
MARKERS = DATA / "markers"
HMMS = DATA / "hmms"
CACHE = MARKERS / ".hit_cache"          # Regenerable, git-ignored

REPO_DB = ROOT.parent / "db" / "gtdb_markers" / "bac120"

MIN_HMM_SCORE = 20.0
N_SPANNING = 20     # Markers in the occupancy-spanning directory
N_CORE = 8          # Complete markers in the fast directory
MIN_TIPS_PER_MARKER = 3   # A 2-tip marker cannot be rooted or tested for clades
RANDOM_SEED = 20260929

CHIMERA_SET = "bench12"
CHIMERA_RECIPIENT = "GCF_000007845.1"   # Bacillus anthracis Ames (Caryophanales)
CHIMERA_DONOR = "GCF_000006865.1"       # Lactococcus lactis Il1403 (Lactobacillales)
N_CHIMERAS = 3


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def read_sets() -> Dict[str, List[str]]:
    if not SETS_FILE.exists():
        raise SystemExit(f"{SETS_FILE} missing — run 01_fetch_genomes.py first.")
    return {
        p.stem: [x.strip() for x in p.read_text(encoding="utf-8").splitlines()
                 if x.strip()]
        for p in sorted(SETS_FILE.glob("*.txt"))
    }


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hmm_files(source: Path, limit: int) -> List[Path]:
    files = sorted(source.glob("*.HMM")) + sorted(source.glob("*.hmm"))
    if not files:
        raise SystemExit(f"no HMM profiles found under {source}")
    return files[:limit] if limit else files


def best_hit(domtblout: Path, min_score: float) -> Optional[Tuple[str, float]]:
    best: Optional[Tuple[str, float]] = None
    for line in domtblout.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 14:
            continue
        try:
            score = float(parts[7])
        except ValueError:
            continue
        if score < min_score:
            continue
        if best is None or score > best[1]:
            best = (parts[0], score)
    return best


def scan_all(hmms: List[Path], sets: Dict[str, List[str]], force: bool) -> None:
    """Fill the (marker, genome) hmmsearch cache for every set."""
    accessions = sorted({a for members in sets.values() for a in members})
    CACHE.mkdir(parents=True, exist_ok=True)
    todo = [
        (hmm, acc) for hmm in hmms for acc in accessions
        if force or not (CACHE / f"{hmm.stem}_{acc}.domtblout").exists()
    ]
    print(f"[scan] {len(todo)} hmmsearch runs to run "
          f"({len(hmms)} profiles x {len(accessions)} proteomes)")
    with tempfile.TemporaryDirectory(prefix="mf_hmmtmp_") as tmp:
        tmp_dir = Path(tmp)
        for i, (hmm, acc) in enumerate(todo, 1):
            dom = CACHE / f"{hmm.stem}_{acc}.domtblout"
            subprocess.run(
                ["hmmsearch", "--noali", "--cpu", "1", "--domtblout", str(dom),
                 str(hmm), str(GENOMES / f"{acc}.faa")],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            if i % 200 == 0:
                print(f"[scan] {i}/{len(todo)}")
    print(f"[scan] cache complete at {CACHE}")


def occupancy(marker_id: str, members: List[str]) -> Dict[str, str]:
    """Accession -> sequence, for one marker within one set."""
    found: Dict[str, str] = {}
    for acc in members:
        dom = CACHE / f"{marker_id}_{acc}.domtblout"
        if not dom.exists():
            continue
        hit = best_hit(dom, MIN_HMM_SCORE)
        if hit is None:
            continue
        seq = _proteomes().get(acc, {}).get(hit[0])
        if seq:
            found[acc] = seq
    return found


_PROTEOME_CACHE: Dict[str, Dict[str, str]] = {}


def _proteomes() -> Dict[str, Dict[str, str]]:
    """Every shipped proteome, parsed once (the scan asks for them ~3000x)."""
    if not _PROTEOME_CACHE:
        for faa in sorted(GENOMES.glob("*.faa")):
            _PROTEOME_CACHE[faa.stem] = read_fasta(faa)
    return _PROTEOME_CACHE


def full_occ_table(hmms: List[Path], sets: Dict[str, List[str]]) -> List[dict]:
    rows = []
    for hmm in hmms:
        entry = {"marker_id": hmm.stem}
        for set_name, members in sets.items():
            n = len(occupancy(hmm.stem, members))
            entry[f"{set_name}_n"] = n
            entry[f"{set_name}_occ"] = round(n / len(members), 4) if members else 0
        rows.append(entry)
    header = ["marker_id"] + [f"{s}_{k}" for s in sets for k in ("n", "occ")]
    lines = ["\t".join(header)]
    for row in sorted(rows, key=lambda r: (-r["small8_occ"], r["marker_id"])):
        lines.append("\t".join(str(row[c]) for c in header))
    MARKERS.mkdir(parents=True, exist_ok=True)
    (MARKERS / "OCCUPANCY_ALL.tsv").write_text("\n".join(lines) + "\n",
                                               encoding="utf-8")
    print(f"[write] {MARKERS / 'OCCUPANCY_ALL.tsv'} ({len(rows)} profiles)")
    return rows


# --------------------------------------------------------------------------
# Directory materialisation
# --------------------------------------------------------------------------
def write_marker_dir(set_name: str, marker_ids: List[str],
                     members: List[str], label: str) -> List[Tuple[str, int]]:
    out_dir = MARKERS / set_name
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    counts: List[Tuple[str, int]] = []
    for marker_id in marker_ids:
        found = occupancy(marker_id, members)
        (out_dir / f"{marker_id}.faa").write_text(
            "\n".join(f">{acc}\n{found[acc]}" for acc in members if acc in found)
            + "\n", encoding="utf-8")
        counts.append((marker_id, len(found)))
    lines = ["marker_id\tsequences_in_file\tgenomes_in_set\toccupancy"]
    for marker_id, count in counts:
        lines.append(f"{marker_id}\t{count}\t{len(members)}\t{count / len(members):.4f}")
    (out_dir / "OCCUPANCY.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    spread = sorted({c for _m, c in counts})
    print(f"[write] {out_dir} ({label}: {len(counts)} markers, "
          f"tip counts {spread[0]}..{spread[-1]})")
    return counts


def choose_markers(rows: List[dict], members_key: str) -> Tuple[List[str], List[str]]:
    """Return (core_ids, spanning_ids) for one set.

    core = the N highest-occupancy markers, all complete: the fast input.
    spanning = a set that reaches down the occupancy range as far as the data
                allows while keeping >= MIN_TIPS_PER_MARKER tips, so thresholds
                have something to reject; padded up to N_SPANNING with complete
                markers.
    """
    usable = [r for r in rows if r[f"{members_key}_n"] >= MIN_TIPS_PER_MARKER]
    usable.sort(key=lambda r: (-r[f"{members_key}_occ"], r["marker_id"]))
    core = [r["marker_id"] for r in usable[:N_CORE]]

    partial = [r for r in usable if r[f"{members_key}_occ"] < 1.0]
    partial.sort(key=lambda r: (r[f"{members_key}_occ"], r["marker_id"]))
    spanning = [r["marker_id"] for r in partial[:N_SPANNING // 2]]
    fill = [r["marker_id"] for r in usable if r["marker_id"] not in spanning]
    spanning = sorted(set(spanning + fill[:N_SPANNING - len(spanning)]))
    return core, spanning


def materialize_hmms(hmms: List[Path], spanning_ids: List[str],
                     core_ids: List[str]) -> None:
    """Profile directories mirroring the two marker-set flavours.

    ``full`` matches ``markers/<set>`` (occupancy-spanning), ``core`` matches
    ``markers/<set>_core`` (complete), so a case can pair a FASTA directory and
    a profile directory without inventing a third combination.
    """
    for subdir, wanted in (("full", spanning_ids), ("core", core_ids)):
        target = HMMS / subdir
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
        for hmm in hmms:
            if hmm.stem in wanted:
                shutil.copy2(hmm, target / hmm.name)
        print(f"[write] {target} ({len(list(target.iterdir()))} profiles)")


# --------------------------------------------------------------------------
# Planted-HGT positive control
# --------------------------------------------------------------------------
def _chimera_dir_is_current(src: Path, dst: Path, truth: Path) -> bool:
    """Is the planted directory still the clean directory plus its swaps?

    Existence is not idempotency: the marker set a clean directory holds can
    change (a different profile-selection rule, a rebuilt genome pool) while the
    planted copy stays byte-for-byte as it was. A paired design built that way
    silently compares two different marker sets, and the detection numbers then
    mean nothing — so the pairing is checked, not assumed.
    """
    clean = {p.name for p in src.glob("*.faa")}
    planted_dir = {p.name for p in dst.glob("*.faa")}
    if clean != planted_dir:
        return False
    try:
        payload = json.loads(truth.read_text(encoding="utf-8"))
        swapped = {e["marker_id"] + ".faa" for e in payload["planted"]}
    except (ValueError, KeyError, TypeError):
        return False
    # Every non-planted file must be byte-identical to its clean counterpart.
    return all((dst / name).read_bytes() == (src / name).read_bytes()
               for name in clean - swapped)


def build_chimeras(members: List[str], force: bool) -> Path:
    src = MARKERS / CHIMERA_SET
    dst = MARKERS / f"chimera_{CHIMERA_SET}"
    truth = MARKERS / f"chimera_ground_truth_{CHIMERA_SET}.json"
    if dst.exists() and truth.exists() and not force:
        if _chimera_dir_is_current(src, dst, truth):
            print(f"[cache] {dst}")
            return truth
        print(f"[stale] {dst} no longer matches {src} — rebuilding")

    candidates = []
    for faa in sorted(src.glob("*.faa")):
        seqs = read_fasta(faa)
        if CHIMERA_RECIPIENT not in seqs or CHIMERA_DONOR not in seqs:
            continue
        if len(seqs[CHIMERA_RECIPIENT]) < 80 or len(seqs[CHIMERA_DONOR]) < 80:
            continue
        if seqs[CHIMERA_RECIPIENT] == seqs[CHIMERA_DONOR]:
            continue
        candidates.append(faa.stem)
    if len(candidates) < N_CHIMERAS:
        raise SystemExit(
            f"only {len(candidates)} markers carry both recipient and donor "
            f"sequences; need >= {N_CHIMERAS}"
        )

    rng = random.Random(RANDOM_SEED)
    selected = sorted(rng.sample(candidates, N_CHIMERAS))

    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

    planted = []
    for marker_id in selected:
        seqs = read_fasta(src / f"{marker_id}.faa")
        recipient_len = len(seqs[CHIMERA_RECIPIENT])
        donor_len = len(seqs[CHIMERA_DONOR])
        seqs[CHIMERA_RECIPIENT] = seqs[CHIMERA_DONOR]
        (dst / f"{marker_id}.faa").write_text(
            "\n".join(f">{g}\n{s}" for g, s in sorted(seqs.items())) + "\n",
            encoding="utf-8",
        )
        planted.append({
            "marker_id": marker_id,
            "recipient_genome": CHIMERA_RECIPIENT,
            "donor_genome": CHIMERA_DONOR,
            "recipient_sequence_length_before": recipient_len,
            "donor_sequence_length": donor_len,
            "replacement": "recipient sequence replaced by the donor ortholog",
            "sha256_of_written_file": sha256_of(dst / f"{marker_id}.faa"),
        })

    payload = {
        "description": (
            "Planted position-known HGT positives. Each listed marker's "
            "recipient sequence was overwritten with the donor ortholog; every "
            "other sequence in the directory is byte-identical to the clean "
            "run, so the clean run doubles as the negative control."
        ),
        "random_seed": RANDOM_SEED,
        "selection_rule": (
            f"{N_CHIMERAS} markers drawn with random.Random({RANDOM_SEED}) from "
            f"the {len(candidates)} markers of {CHIMERA_SET} that contain both "
            "recipient and donor and whose two sequences differ"
        ),
        "candidate_marker_pool": candidates,
        "recipient": {
            "accession": CHIMERA_RECIPIENT,
            "organism": "Bacillus anthracis str. Ames",
            "order": "Caryophanales",
        },
        "donor": {
            "accession": CHIMERA_DONOR,
            "organism": "Lactococcus lactis subsp. lactis Il1403",
            "order": "Lactobacillales",
        },
        "n_markers_in_directory": len(list(dst.glob("*.faa"))),
        "planted": planted,
    }
    truth.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"[write] {truth} ({len(planted)} planted of {len(candidates)} candidates)")
    return truth


# --------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true", help="rebuild every artefact")
    ap.add_argument("--hmm-dir", default=str(REPO_DB),
                    help="directory of individual TIGRFAM/Pfam profiles "
                         "(default: the repository's db/gtdb_markers/bac120)")
    ap.add_argument("--max-profiles", type=int, default=0,
                    help="limit the scanned profiles (0 = all found)")
    args = ap.parse_args(argv)

    sets = read_sets()
    hmms = hmm_files(Path(args.hmm_dir), args.max_profiles)
    marker_sets = {k: v for k, v in sets.items() if k != "single1"}
    print(f"[source] {len(hmms)} profiles from {args.hmm_dir}; "
          f"sets: {', '.join(sorted(marker_sets))}")

    scan_all(hmms, marker_sets, args.force)
    rows = full_occ_table(hmms, marker_sets)

    core_ids_all: List[str] = []
    spanning_all: List[str] = []
    for set_name, members in sorted(marker_sets.items()):
        core, spanning = choose_markers(rows, set_name)
        if set_name == "small8":
            core_ids_all, spanning_all = core, spanning
        write_marker_dir(f"{set_name}_core", core, members, "complete markers")
        write_marker_dir(set_name, spanning, members, "occupancy-spanning")

    # Degenerate inputs: the bundled library is a conserved-marker library, so
    # No profile is sparse on these closely related proteomes (the sparsest is
    # 7/8). A 2-tip and a 1-tip marker are therefore carved out of real
    # Sequences by subsetting a complete one -- recorded as such, and used only
    # By the robustness cases.
    degenerate = MARKERS / "degenerate"
    if degenerate.exists():
        shutil.rmtree(degenerate)
    degenerate.mkdir(parents=True)
    source = MARKERS / "small8_core" / f"{core_ids_all[0]}.faa"
    lines = source.read_text(encoding="utf-8").splitlines()
    records = [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1)
               if lines[i].startswith(">")]
    two_tip = "\n".join(f">{name}\n{seq}" for name, seq in records[:2]) + "\n"
    one_tip = "\n".join(f">{name}\n{seq}" for name, seq in records[:1]) + "\n"
    (degenerate / f"twotip_{core_ids_all[0]}.faa").write_text(two_tip,
                                                             encoding="utf-8")
    (degenerate / f"onetip_{core_ids_all[0]}.faa").write_text(one_tip,
                                                             encoding="utf-8")
    (degenerate / "README.txt").write_text(
        "Subsets carved out of the real "
        f"{core_ids_all[0]} marker file (small8_core) to hold exactly 2 and 1 "
        "tips. They are not observed absences: they exist to check that a "
        "below-minimally-informative marker file is reported rather than "
        "silently producing a tree.\n", encoding="utf-8")
    print(f"[write] {degenerate} (2-tip and 1-tip markers)")

    materialize_hmms(hmms, spanning_all, core_ids_all)
    build_chimeras(marker_sets[CHIMERA_SET], args.force)
    print("[done] core=", ",".join(core_ids_all),
          "\n       spanning=", ",".join(spanning_all))
    return 0


if __name__ == "__main__":
    sys.exit(main())
