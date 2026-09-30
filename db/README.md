# Reference databases

This directory holds the metadata MarkerFinder verifies its databases against.
The profile data itself is **fetched, not committed**.

| Path | In git? | What it is |
|---|---|---|
| `expected_hashes.json` | yes | The SHA-256 each database is expected to have. `markerfinder --check` compares the directory on disk against it; an empty value is reported as `NOT RUN`, never as a pass. |
| `gtdb_markers/` | **no** (git-ignored) | One HMM profile per marker, in `bac120/` (bacteria) and `ar53/` (archaea). Used by `--marker-mode hmm` when `--marker-hmm-dir` is omitted. |

## Getting the HMM profiles

```bash
python3 scripts/fetch_marker_db.py --auto            # use a GTDB-Tk install here
python3 scripts/fetch_marker_db.py --set bac120 --from-hmm /path/gtdbtk_bac120.a.hmm
python3 scripts/fetch_marker_db.py --set ar53 --from-dir /path/to/ar53_marker_genes
```

`--marker-mode gtdb_tk` does not need any of this: it reads the marker FASTA
files that `gtdbtk align` already produced, passed with `--gtdb-markers-dir`.

## Why the profiles are not redistributed here

The models are third-party work — TIGRFAM profiles from NCBI and Pfam profiles
from InterPro — assembled into per-marker sets by GTDB-Tk. Committing them would
put their redistribution terms on this repository, so the tool ships the
assembler and the integrity check instead, and the user's own GTDB-Tk install
supplies the data.

Licensing of the upstream resources:

- Pfam / InterPro models: <https://www.interpro.org/terms/>
- NCBI TIGRFAM models: <https://www.ncbi.nlm.nih.gov/pubs/factsheet/>
- GTDB-Tk (which packages both): <https://github.com/ChanBurton/GTDB-Tk>

After fetching, verify what you got:

```bash
markerfinder --check                       # reports the computed and expected hashes
markerfinder --check --db-dir X            # checks the directory you actually use
```

`fetch_marker_db.py --record` updates `expected_hashes.json` with the hash of the
profiles on disk, which is how you pin a deliberate database change.
