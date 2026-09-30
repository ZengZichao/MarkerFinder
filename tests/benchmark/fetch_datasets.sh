#!/usr/bin/env bash
# fetch the benchmark genome set from NCBI.
#
# NOT auto-executed anywhere: GB-level download, user-triggered.
# Downloads land in tests/benchmark/downloads/ (gitignored); a
# SHA-256 manifest is written for reproducibility (pairs with ).
#
# Refuses to run while expected/accessions.yaml still carries the
# REPLACE_WITH_ORDER_NAME placeholder, or when disk headroom < 5 GB.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ACCESSIONS="$HERE/expected/accessions.yaml"
OUT="$HERE/downloads"

if grep -q "REPLACE_WITH_ORDER_NAME" "$ACCESSIONS"; then
  echo "REFUSING: expected/accessions.yaml still contains the placeholder order." >&2
  echo "Choose the concrete GTDB ar53 order, fill accessions, then re-run." >&2
  exit 2
fi

AVAIL_KB=$(df -k "$HERE" | awk 'NR==2 {print $4}')
if [ "$AVAIL_KB" -lt $((5 * 1024 * 1024)) ]; then
  echo "REFUSING: <5GB disk headroom (have ${AVAIL_KB}KB)." >&2
  exit 2
fi

mkdir -p "$OUT"
MANIFEST="$OUT/manifest.sha256"
: > "$MANIFEST"

# Extract GCF_* accessions from the YAML (no PyYAML dependency: plain grep).
ACCS=$(grep -oE 'GCF_[0-9]+\.[0-9]+' "$ACCESSIONS" | sort -u)
COUNT=$(printf '%s\n' "$ACCS" | wc -l | tr -d ' ')
echo "Fetching $COUNT genomes into $OUT"

for acc in $ACCS; do
  file="$OUT/${acc}_genomic.faa.gz"
  if [ -s "$file" ]; then
    echo " skip (exists): $acc"
  else
    url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${acc:0:4}/${acc:4:3}/${acc:7:3}/${acc:9:3}/${acc}/${acc}_protein.faa.gz"
    echo " downloading $acc"
    curl -fsSL --retry 3 -o "$file" "$url"
  fi
  shasum -a 256 "$file" >> "$MANIFEST"
done

echo "Done. Manifest: $MANIFEST"
echo "Next: python tests/benchmark/make_chimeras.py --downloads $OUT"
