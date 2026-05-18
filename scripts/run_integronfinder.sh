#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

THREADS="${THREADS:-16}"
OUTDIR="$ROOT/results/mobile_elements/integronfinder"
INPUT_DIR="$OUTDIR/input"
LOGDIR="$ROOT/logs/mobile_elements"
INPUT="$INPUT_DIR/BC42_BC88_polished.prefixed.fasta"

BC42="$ROOT/results/polished/medaka/BC42/consensus.fasta"
BC88="$ROOT/results/polished/medaka/BC88/consensus.fasta"

mkdir -p "$INPUT_DIR" "$OUTDIR" "$LOGDIR"

if ! command -v integron_finder >/dev/null 2>&1; then
  echo "ERROR: integron_finder not found. Activate the pdr-pa-mobile environment first." >&2
  exit 1
fi

for fasta in "$BC42" "$BC88"; do
  if [[ ! -s "$fasta" ]]; then
    echo "ERROR: Missing input FASTA: $fasta" >&2
    exit 1
  fi
done

{
  awk '/^>/ {sub(/^>/, ">BC42|"); print; next} {print}' "$BC42"
  awk '/^>/ {sub(/^>/, ">BC88|"); print; next} {print}' "$BC88"
} > "$INPUT"

(
  cd "$OUTDIR"
  integron_finder --cpu "$THREADS" "$INPUT" 2>&1 | tee "$LOGDIR/integronfinder_BC42_BC88.log"
)

echo "Done. IntegronFinder output under: $OUTDIR"

