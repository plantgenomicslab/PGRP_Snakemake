#!/usr/bin/env bash
# Minimal smoke test for the BBDuk pre-filter:
#   - Subsamples 100k reads from user-supplied paired FASTQs
#   - Runs bbduk.sh directly with the project reference set
#   - Asserts non-empty clean output and refstats
#
# Requires:
#   - Reference set already built (bash scripts/build_bbduk_refs.sh)
#   - bbmap conda env active (bbduk.sh on PATH)
#   - Env vars SMOKE_R1 and SMOKE_R2 pointing at paired-end FASTQ.GZ files

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"

if [ ! -f data/bbduk_refs/cp_12sp.fa.gz ] \
    || [ ! -f data/bbduk_refs/mt_12sp.fa.gz ] \
    || [ ! -f data/bbduk_refs/univec.fa.gz ]; then
  echo "Refs missing under data/bbduk_refs/. Run scripts/build_bbduk_refs.sh first." >&2
  exit 1
fi

if ! command -v bbduk.sh >/dev/null 2>&1; then
  echo "bbduk.sh not on PATH. Activate the bbmap conda env (envs/bbduk.yaml)." >&2
  exit 1
fi

: "${SMOKE_R1:?Set SMOKE_R1 to a paired R1 FASTQ.GZ}"
: "${SMOKE_R2:?Set SMOKE_R2 to a paired R2 FASTQ.GZ}"

WORK="tests/_workdir"
rm -rf "$WORK"
mkdir -p "$WORK"

# Take first 100k read pairs (4 lines/read => 400000 lines)
zcat "$SMOKE_R1" | head -400000 | gzip > "$WORK/R1.fq.gz"
zcat "$SMOKE_R2" | head -400000 | gzip > "$WORK/R2.fq.gz"

# Resolve ribokmers via the same logic the Snakemake rule uses
RIBO=$(python3 - <<'PY'
import shutil, glob, os, sys
b = shutil.which("bbduk.sh")
if not b:
    sys.exit("bbduk.sh not on PATH")
base = os.path.dirname(b)
for layout in ("opt/bbmap-*", "share/bbmap-*"):
    cands = glob.glob(os.path.join(base, "..", layout, "resources", "ribokmers.fa.gz"))
    if cands:
        print(cands[0])
        break
else:
    sys.exit("ribokmers.fa.gz not found")
PY
)
echo "ribokmers: $RIBO"

bbduk.sh \
    in1="$WORK/R1.fq.gz" in2="$WORK/R2.fq.gz" \
    out1="$WORK/clean_1.fq.gz" out2="$WORK/clean_2.fq.gz" \
    outm1="$WORK/contam_1.fq.gz" outm2="$WORK/contam_2.fq.gz" \
    ref="$RIBO,data/bbduk_refs/univec.fa.gz,data/bbduk_refs/cp_12sp.fa.gz,data/bbduk_refs/mt_12sp.fa.gz" \
    k=31 mcf=0.5 \
    stats="$WORK/stats.txt" refstats="$WORK/refstats.txt" \
    threads=4 -Xmx8g

echo "=== refstats ==="
cat "$WORK/refstats.txt"

[ -s "$WORK/clean_1.fq.gz" ] || { echo "FAIL: clean_1.fq.gz empty"; exit 1; }
[ -s "$WORK/clean_2.fq.gz" ] || { echo "FAIL: clean_2.fq.gz empty"; exit 1; }
[ -s "$WORK/refstats.txt" ]  || { echo "FAIL: refstats.txt empty"; exit 1; }
echo "PASS"
