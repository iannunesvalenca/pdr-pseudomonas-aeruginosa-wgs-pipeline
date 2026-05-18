#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

THREADS="${THREADS:-16}"
OUTDIR="$ROOT/results/annotation/bakta"
LOGDIR="$ROOT/logs/annotation"

mkdir -p "$OUTDIR" "$LOGDIR"

if ! command -v bakta >/dev/null 2>&1; then
  echo "ERROR: bakta not found. Activate the pdr-pa-bakta environment first." >&2
  exit 1
fi

if [[ -z "${BAKTA_DB:-}" ]]; then
  if [[ -d "$ROOT/databases/bakta/db" ]] && find "$ROOT/databases/bakta/db" -mindepth 1 -print -quit | grep -q .; then
    BAKTA_DB="$ROOT/databases/bakta/db"
  elif [[ -d "$ROOT/databases/bakta" ]] && find "$ROOT/databases/bakta" -mindepth 1 -print -quit | grep -q .; then
    BAKTA_DB="$ROOT/databases/bakta"
  else
    BAKTA_DB="$ROOT/databases/bakta/db"
  fi
fi

if [[ ! -d "$BAKTA_DB" ]] || ! find "$BAKTA_DB" -mindepth 1 -print -quit | grep -q .; then
  echo "ERROR: Bakta database not found at: $BAKTA_DB" >&2
  echo "Download first, for example:" >&2
  echo "  bakta_db download --output $ROOT/databases/bakta" >&2
  echo "Then rerun this script. If the database is placed in another directory, set:" >&2
  echo "  export BAKTA_DB=/path/to/bakta/db" >&2
  exit 1
fi

rm -rf "$OUTDIR/BC42" "$OUTDIR/BC88"

bakta \
  --db "$BAKTA_DB" \
  --threads "$THREADS" \
  --keep-contig-headers \
  --output "$OUTDIR/BC42" \
  --prefix BC42 \
  "$ROOT/results/polished/medaka/BC42/consensus.fasta" \
  2>&1 | tee "$LOGDIR/BC42_bakta.log"

bakta \
  --db "$BAKTA_DB" \
  --threads "$THREADS" \
  --keep-contig-headers \
  --output "$OUTDIR/BC88" \
  --prefix BC88 \
  "$ROOT/results/polished/medaka/BC88/consensus.fasta" \
  2>&1 | tee "$LOGDIR/BC88_bakta.log"

echo "Done. Bakta output under: $OUTDIR"
