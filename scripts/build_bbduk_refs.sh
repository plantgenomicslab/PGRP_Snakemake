#!/usr/bin/env bash
# Download 12 chloroplast + 12 mitochondrion + UniVec; concatenate into BBDuk reference set.
# Requires: efetch (entrez-direct), curl, gzip.
# Output: data/bbduk_refs/{cp_12sp.fa.gz,mt_12sp.fa.gz,univec.fa.gz}
# ribokmers.fa.gz is auto-located in the BBTools install at runtime by the Snakemake rule.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF_DIR="$REPO_ROOT/data/bbduk_refs"
LIST="$REF_DIR/SPECIES_LIST.tsv"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

if [ ! -f "$LIST" ]; then
  echo "ERROR: $LIST missing" >&2
  exit 1
fi

if ! command -v efetch >/dev/null 2>&1; then
  echo "ERROR: efetch not on PATH. Install entrez-direct: mamba install -c bioconda entrez-direct" >&2
  exit 1
fi

retry() {
  # retry CMD ARGS...  — up to 3 attempts, exponential backoff 5/10/15s
  local n=0 max=3 delay=5
  until "$@"; do
    n=$((n + 1))
    if [ "$n" -ge "$max" ]; then
      echo "ERROR: command failed after $max attempts: $*" >&2
      return 1
    fi
    echo "  retry $n/$max in $((delay * n))s: $*" >&2
    sleep $((delay * n))
  done
}

fetch_cp_set() {
  # NCBI recommends epost | efetch batching for >3 IDs.
  awk -F'\t' 'NR>1 && $3 != "" {print $3}' "$LIST" \
    | epost -db nucleotide \
    | efetch -format fasta > cp_12sp.fa.tmp
}

fetch_mt_set() {
  awk -F'\t' 'NR>1 && $4 != "" {print $4}' "$LIST" \
    | epost -db nucleotide \
    | efetch -format fasta > mt_12sp.fa.tmp
}

fetch_univec() {
  curl -fsSL https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec > univec.fa.tmp
}

echo "Fetching chloroplast set (cp_12sp.fa) via epost | efetch..."
retry fetch_cp_set
mv cp_12sp.fa.tmp cp_12sp.fa
gzip -f cp_12sp.fa
echo "cp_12sp.fa.gz built: $(zcat cp_12sp.fa.gz | grep -c '^>') sequences"

echo "Fetching mitochondrion set (mt_12sp.fa) via epost | efetch..."
retry fetch_mt_set
mv mt_12sp.fa.tmp mt_12sp.fa
gzip -f mt_12sp.fa
echo "mt_12sp.fa.gz built: $(zcat mt_12sp.fa.gz | grep -c '^>') sequences"

echo "Fetching UniVec..."
retry fetch_univec
mv univec.fa.tmp univec.fa
gzip -f univec.fa
echo "univec.fa.gz built: $(zcat univec.fa.gz | grep -c '^>') sequences"

echo
echo "Done. Reference set in $REF_DIR/"
ls -lh cp_12sp.fa.gz mt_12sp.fa.gz univec.fa.gz
