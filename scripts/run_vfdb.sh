#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-8}"
AMR_ENV="${AMR_ENV:-}"

if [[ -n "$AMR_ENV" ]]; then
  export PATH="$AMR_ENV/bin:$PATH"
  export PERL5LIB="$AMR_ENV/lib/perl5/vendor_perl${PERL5LIB:+:$PERL5LIB}"
fi

OUTDIR="$ROOT/results/virulence/vfdb"
LOGDIR="$ROOT/logs/virulence"
mkdir -p "$OUTDIR" "$LOGDIR"

for SAMPLE in BC42 BC88; do
  ASM="$ROOT/results/polished/medaka/$SAMPLE/consensus.fasta"
  OUT="$OUTDIR/${SAMPLE}_vfdb.tsv"

  if [[ ! -s "$ASM" ]]; then
    echo "ERROR: assembly not found for $SAMPLE: $ASM" >&2
    exit 1
  fi

  abricate \
    --db vfdb \
    --minid 80 \
    --mincov 80 \
    --threads "$THREADS" \
    "$ASM" \
    > "$OUT" \
    2> "$LOGDIR/${SAMPLE}_vfdb.log"
done

abricate --summary "$OUTDIR/BC42_vfdb.tsv" "$OUTDIR/BC88_vfdb.tsv" \
  > "$OUTDIR/BC42_BC88_vfdb_summary.tsv" \
  2> "$LOGDIR/BC42_BC88_vfdb_summary.log"

awk '
  BEGIN { IGNORECASE = 1 }
  FNR == 1 && NR == 1 { print; next }
  FNR == 1 { next }
  $0 ~ /(exo[STUY]|toxA|type III|T3SS|psc|pcr|popB|popD|pcrV|pil|alg|las|rhl)/ { print }
' "$OUTDIR/BC42_vfdb.tsv" "$OUTDIR/BC88_vfdb.tsv" \
  > "$OUTDIR/BC42_BC88_vfdb_focused.tsv"

echo "VFDB summary: $OUTDIR/BC42_BC88_vfdb_summary.tsv"
echo "Focused VFDB hits: $OUTDIR/BC42_BC88_vfdb_focused.tsv"
