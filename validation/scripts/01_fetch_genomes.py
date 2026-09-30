#!/usr/bin/env python3
"""Step 1 of the MarkerFinder validation data pipeline: real genomes + provenance.

What it produces
----------------
``data/genomes/<accession>.faa``
    Predicted-protein FASTA of each RefSeq assembly, taken from NCBI.
``data/PROVENANCE.tsv``
    One row per assembly: organism, tax_id, the seven NCBI lineage ranks,
    assembly statistics, NCBI's own CheckM completeness/contamination, the
    download URL, the retrieval date (UTC) and the SHA-256 of the shipped file.

Why the lineage is fetched and not typed
----------------------------------------
The taxonomy table drives the phylogenetic HGT screen (MAD rooting +
monophyly proportion), so a wrong rank is not a typo, it is a wrong
experiment. An earlier revision of this dataset hand-typed the lineages and
labelled ``GCF_000008105.1`` as ``o__Bacillales;f__Bacillaceae;g__Bacillus``;
NCBI reports it as *Salmonella enterica* serovar Choleraesuis (Pseudomonadota,
Enterobacterales, Enterobacteriaceae). Every rank in ``PROVENANCE.tsv`` and in
``data/taxonomy/*.tsv`` is therefore read from ``ncbi-datasets`` at build time
and carries its ``tax_id``, so a reader can re-derive it.

Usage
-----
    python 01_fetch_genomes.py # metadata + missing genomes
    python 01_fetch_genomes.py --metadata-only # never touch the FASTA files
    python 01_fetch_genomes.py --force # re-download everything

Requires the NCBI datasets CLI (``datasets``) on PATH, or set
``MARKERFINDER_DATASETS_BIN`` to its absolute path.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent                      # Validation/
DATA = ROOT / "data"
GENOMES = DATA / "genomes"
PROVENANCE = DATA / "PROVENANCE.tsv"

# --------------------------------------------------------------------------
# The dataset. Two nested sets, both spanning Bacillota (2 orders) plus one
# Pseudomonadota outgroup, so that cross-order incongruence is present in the
# Data rather than only asserted about it.
# --------------------------------------------------------------------------
ACCESSIONS: List[str] = [
    # --- order Caryophanales (formerly Bacillales), phylum Bacillota --------
    "GCF_000009045.1",   # Bacillus subtilis subsp. subtilis str. 168
    "GCF_000007825.1",   # Bacillus cereus ATCC 14579
    "GCF_000008445.1",   # Bacillus anthracis str. 'Ames Ancestor'
    "GCF_000007845.1",   # Bacillus anthracis str. Ames
    "GCF_000008505.1",   # Bacillus thuringiensis serovar konkukian str. 97-27
    "GCF_000009645.1",   # Staphylococcus aureus subsp. aureus N315
    "GCF_000196035.1",   # Listeria monocytogenes EGD-e
    # --- order Lactobacillales, phylum Bacillota ----------------------------
    "GCF_000007785.1",   # Enterococcus faecalis V583
    "GCF_000006785.1",   # Streptococcus pyogenes M1 GAS
    "GCF_000006865.1",   # Lactococcus lactis subsp. lactis Il1403
    "GCF_009913655.1",   # Lactiplantibacillus plantarum WCFS1
    # --- order Enterobacterales, phylum Pseudomonadota (distant outgroup) ---
    "GCF_000008105.1",   # Salmonella enterica serovar Choleraesuis str. SC-B67
]

# Membership of each named input set. ``quad4`` is the fast default for
# Feature-level cases: the phylogenetic step needs >= 4 tips per marker file
# (below that the pipeline refuses to build a gene tree and exits 3), and four
# Genomes keep one run to a couple of minutes.
SETS: Dict[str, List[str]] = {
    "quad4": [
        "GCF_000009045.1", "GCF_000007825.1", "GCF_000009645.1",
        "GCF_000196035.1",
    ],
    "small8": [
        "GCF_000009045.1", "GCF_000007825.1", "GCF_000008445.1",
        "GCF_000196035.1", "GCF_000009645.1", "GCF_000007785.1",
        "GCF_000006785.1", "GCF_000008105.1",
    ],
    "bench12": list(ACCESSIONS),
    # Below the 4-tip minimum: kept deliberately so the refusal is a tested
    # Behaviour rather than a rumour.
    "pair3": [
        "GCF_000009045.1", "GCF_000007825.1", "GCF_000009645.1",
    ],
    "single1": ["GCF_000009045.1"],
}

LINEAGE_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# RefSeq assembly accessions are strictly "GCF_" + 9 digits + "." + version.
# Enforcing the shape keeps every downstream path construction inside the
# data directory even if ACCESSIONS ever becomes file- or CLI-fed.
ACCESSION_RE = re.compile(r"GCF_[0-9]{9}\.[0-9]+")


def safe_genome_path(acc: str) -> Path:
    """Resolve ``GENOMES / <acc>.faa`` and refuse anything escaping GENOMES."""
    if not ACCESSION_RE.fullmatch(acc):
        raise SystemExit(f"unexpected accession format: {acc!r}")
    candidate = (GENOMES / f"{acc}.faa").resolve()
    if not candidate.is_relative_to(GENOMES.resolve()):
        raise SystemExit(f"path escapes the genome directory: {acc!r}")
    return candidate

PROVENANCE_COLUMNS = [
    "accession", "organism_name", "tax_id",
    "domain", "phylum", "class", "order", "family", "genus", "species",
    "assembly_level", "number_of_contigs", "total_sequence_length", "gc_percent",
    "protein_coding_genes", "checkm_completeness", "checkm_contamination",
    "refseq_release_date", "source_registry", "license",
    "download_url", "retrieved_utc", "file", "sha256",
]


def datasets_bin() -> str:
    override = os.environ.get("MARKERFINDER_DATASETS_BIN")
    if override:
        if not Path(override).exists():
            raise SystemExit(
                f"MARKERFINDER_DATASETS_BIN points at a non-existent file: {override}"
            )
        return override
    found = shutil.which("datasets")
    if not found:
        raise SystemExit(
            "NCBI datasets CLI not found on PATH. Install it "
            "(micromamba install -c bioconda ncbi-datasets-cli) or set "
            "MARKERFINDER_DATASETS_BIN=/abs/path/to/datasets."
        )
    return found


def _run_json(cmd: List[str]) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr[:500]}"
        )
    return json.loads(proc.stdout)


def fetch_metadata(acc: str, ds: str) -> Dict[str, str]:
    """Organism, lineage and assembly statistics for one accession."""
    genome = _run_json([ds, "summary", "genome", "accession", acc, "--pretty"])
    report = genome["reports"][0]
    organism = report["organism"]
    tax_id = organism["tax_id"]

    taxonomy = _run_json([ds, "summary", "taxonomy", "taxon", str(tax_id), "--pretty"])
    classification = taxonomy["reports"][0]["taxonomy"]["classification"]

    stats = report.get("assembly_stats", {})
    assembly = report.get("assembly_info", {})
    annotation = report.get("annotation_info", {}).get("stats", {}).get("gene_counts", {})
    checkm = report.get("checkm_info", {})

    row = {
        "accession": acc,
        "organism_name": organism.get("organism_name", ""),
        "tax_id": str(tax_id),
        "assembly_level": assembly.get("assembly_level", ""),
        "number_of_contigs": str(stats.get("number_of_contigs", "")),
        "total_sequence_length": str(stats.get("total_sequence_length", "")),
        "gc_percent": str(stats.get("gc_percent", "")),
        "protein_coding_genes": str(annotation.get("protein_coding", "")),
        "checkm_completeness": str(checkm.get("completeness", "")),
        "checkm_contamination": str(checkm.get("contamination", "")),
        "refseq_release_date": report.get("annotation_info", {}).get("release_date", ""),
        "source_registry": "NCBI RefSeq",
        "license": "NCBI RefSeq sequences are free to redistribute (no known "
                   "copyright restrictions)",
        "download_url": f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/"
                       f"{acc}/download",
    }
    for rank in LINEAGE_RANKS:
        node = classification.get(rank)
        row[rank] = node["name"] if node else ""
    return row


def download_protein_faa(acc: str, ds: str, dest: Path, workdir: Path) -> None:
    """Fetch ``protein.faa`` for one accession out of the datasets zip."""
    zip_path = workdir / f"{acc}.zip"
    subprocess.run(
        [ds, "download", "genome", "accession", acc,
         "--include", "protein", "--no-progressbar",
         "--filename", str(zip_path)],
        check=True, capture_output=True, text=True,
    )
    with zipfile.ZipFile(zip_path) as zf:
        candidates = [n for n in zf.namelist() if n.endswith("protein.faa")]
        if not candidates:
            raise RuntimeError(f"no protein.faa inside {zip_path}")
        chosen = next((n for n in candidates if f"/{acc}/" in n), candidates[0])
        with zf.open(chosen) as fh:
            payload = fh.read()
    dest.write_bytes(payload)
    zip_path.unlink(missing_ok=True)


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def read_existing_provenance(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    out: Dict[str, Dict[str, str]] = {}
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(header):
            continue
        out[cells[0]] = dict(zip(header, cells))
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true",
                    help="re-download every genome even if the file exists")
    ap.add_argument("--metadata-only", action="store_true",
                    help="refresh PROVENANCE.tsv only; never touch the FASTA files")
    ap.add_argument("--offline", action="store_true",
                    help="reuse the recorded PROVENANCE.tsv, contact nothing")
    args = ap.parse_args(argv)

    GENOMES.mkdir(parents=True, exist_ok=True)
    workdir = DATA / ".download_tmp"
    workdir.mkdir(exist_ok=True)

    if args.offline:
        print("[offline] PROVENANCE.tsv left as recorded; no network access.")
        return 0 if PROVENANCE.exists() else 1

    ds = datasets_bin()
    print(f"[datasets] using {ds}")

    existing = read_existing_provenance(PROVENANCE)
    today = _dt.date.today().isoformat()
    rows: List[Dict[str, str]] = []

    for acc in ACCESSIONS:
        faa = safe_genome_path(acc)
        if faa.exists() and faa.stat().st_size > 0 and not args.force:
            print(f"[cache] {faa.name}")
        elif args.metadata_only:
            raise SystemExit(
                f"{faa} is missing and --metadata-only was given: run without "
                "it first so the genome FASTA can be downloaded."
            )
        else:
            print(f"[download] {acc}")
            download_protein_faa(acc, ds, faa, workdir)

        try:
            row = fetch_metadata(acc, ds)
        except Exception as e:  # Noqa: BLE001 - keep the older record, say so
            print(f"[metadata] {acc}: refresh failed ({e}); reusing recorded row",
                  file=sys.stderr)
            row = dict(existing.get(acc) or {"accession": acc})
        row["file"] = f"genomes/{faa.name}"
        row["sha256"] = sha256_of(faa)
        row["retrieved_utc"] = existing.get(acc, {}).get("retrieved_utc") or today
        rows.append(row)
        print(f"[ok] {acc}  {row.get('organism_name', '?')}  "
              f"{row.get('order', '?')} / {row.get('genus', '?')}")

    shutil.rmtree(workdir, ignore_errors=True)

    provenance_path = (DATA / "PROVENANCE.tsv").resolve()
    if not provenance_path.is_relative_to(DATA.resolve()):
        raise SystemExit("PROVENANCE path escapes the validation data directory")
    with provenance_path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(PROVENANCE_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in PROVENANCE_COLUMNS) + "\n")
    print(f"[write] {PROVENANCE} ({len(rows)} assemblies)")

    sets_dir = DATA / "sets"
    sets_dir.mkdir(exist_ok=True)
    for name, members in SETS.items():
        (sets_dir / f"{name}.txt").write_text("\n".join(members) + "\n", encoding="utf-8")
        print(f"[write] {sets_dir / (name + '.txt')} ({len(members)} genomes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
