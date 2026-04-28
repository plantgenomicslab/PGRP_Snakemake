#!/usr/bin/env bash
# CI helper: snakemake -n against a sandboxed config so the workflow
# can be validated without real genome refs / fastq / SRA access.
#
# Usage:  bash tests/ci_dry_run.sh {true|false}
#   $1 = bbduk_enable value to inject into the sandbox config

set -euo pipefail

BBDUK_ENABLE="${1:?Pass 'true' or 'false' as the first argument}"
case "$BBDUK_ENABLE" in
  true|false) ;;
  *) echo "ERROR: bbduk_enable must be 'true' or 'false'" >&2; exit 2 ;;
esac

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SANDBOX="$(mktemp -d)"
trap 'rm -rf "$SANDBOX"' EXIT

# Stub everything Snakefile parse-time checks expect
mkdir -p "$SANDBOX/refdir" "$SANDBOX/raw" "$SANDBOX/bbduk_refs"
touch "$SANDBOX/refdir/SA" "$SANDBOX/refdir/SAindex" \
      "$SANDBOX/refdir/exonGeTrInfo.tab" "$SANDBOX/refdir/genome.gtf" \
      "$SANDBOX/refdir/rsem_prep.seq"
touch "$SANDBOX/bbduk_refs/cp_12sp.fa.gz" \
      "$SANDBOX/bbduk_refs/mt_12sp.fa.gz" \
      "$SANDBOX/bbduk_refs/univec.fa.gz" \
      "$SANDBOX/bbduk_refs/ribokmers.fa.gz"

cat > "$SANDBOX/RunsByExperiment.tsv" <<'EOF'
Run	Experiment	Replicate	Treatment
S1_rep1	exp1	S1_rep1	A
S2_rep1	exp1	S2_rep1	B
EOF

# Build sandbox config from the canonical config.yml + overrides.
python3 - "$REPO/config.yml" "$SANDBOX" "$BBDUK_ENABLE" <<'PY'
import sys, yaml
src, sandbox, enable = sys.argv[1], sys.argv[2], sys.argv[3]
c = yaml.safe_load(open(src))
c.update({
    "source": "local",
    "rawInputDir": f"{sandbox}/raw",
    "libraryType": "PAIRED",
    "PAIR_LIST": ["_1", "_2"],
    "genomeDir": f"{sandbox}/refdir",
    "GTFname": f"{sandbox}/refdir/genome.gtf",
    "RSEM_prepared_genome": f"{sandbox}/refdir/rsem_prep",
    "readCounting": ["featureCounts"],
    "runDEG": False,
    "sample_contrast": f"{sandbox}/no_contrast.txt",
    "bbduk_enable": enable.lower() == "true",
    "bbduk_refs_dir": f"{sandbox}/bbduk_refs",
    "bbduk_ribokmers": f"{sandbox}/bbduk_refs/ribokmers.fa.gz",
    "bbduk_bin": "/bin/true",
})
yaml.dump(c, open(f"{sandbox}/config.yml", "w"), sort_keys=False)
PY

ln -sf "$REPO/Snakefile" "$SANDBOX/"
ln -sf "$REPO/scripts"   "$SANDBOX/"
ln -sf "$REPO/rules"     "$SANDBOX/"

echo "=== snakemake -n (bbduk_enable=$BBDUK_ENABLE) ==="
snakemake -n --scheduler greedy -s "$SANDBOX/Snakefile" -d "$SANDBOX" 2>&1 \
  | tee "$SANDBOX/dry_run.log"

# Assert: bbduk_filter_PAIRED appears iff bbduk_enable=true
if [ "$BBDUK_ENABLE" = "true" ]; then
  grep -q "^bbduk_filter_PAIRED" "$SANDBOX/dry_run.log" \
    || { echo "FAIL: expected bbduk_filter_PAIRED in dry-run output"; exit 1; }
else
  ! grep -q "^bbduk_filter_PAIRED" "$SANDBOX/dry_run.log" \
    || { echo "FAIL: bbduk_filter_PAIRED scheduled despite bbduk_enable=false"; exit 1; }
fi

# Always-required jobs (rule names match the post-Tier-2 unified layout)
for rule in align trim_PAIRED featureCounts; do
  grep -q "^$rule" "$SANDBOX/dry_run.log" \
    || { echo "FAIL: $rule missing from dry-run output"; exit 1; }
done

echo "PASS: dry-run for bbduk_enable=$BBDUK_ENABLE"
