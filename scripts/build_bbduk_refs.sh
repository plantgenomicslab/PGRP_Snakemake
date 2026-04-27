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

fetch_set() {
  # $1 = column index (3 = cp_acc, 4 = mt_acc), $2 = output basename (without .gz)
  local col="$1" out="$2"
  : > "${out}.tmp"
  awk -F'\t' -v c="$col" 'NR>1 && $c != "" {print $c}' "$LIST" | while read -r acc; do
    echo "  fetch $acc -> $out"
    # 1-second courtesy spacing for NCBI Entrez (24 sequential calls in this script)
    efetch -db nucleotide -id "$acc" -format fasta >> "${out}.tmp"
    sleep 1
  done
  mv "${out}.tmp" "$out"
  gzip -f "$out"
  echo "${out}.gz built: $(zcat "${out}.gz" | grep -c '^>') sequences"
}

fetch_set 3 cp_12sp.fa
fetch_set 4 mt_12sp.fa

curl -fsSL https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec > univec.fa
gzip -f univec.fa
echo "univec.fa.gz built: $(zcat univec.fa.gz | grep -c '^>') sequences"

echo
echo "Done. Reference set in $REF_DIR/"
ls -lh cp_12sp.fa.gz mt_12sp.fa.gz univec.fa.gz
