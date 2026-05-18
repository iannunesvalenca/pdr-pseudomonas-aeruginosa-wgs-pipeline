#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-8}"
AMR_ENV="${AMR_ENV:-}"

if [[ -n "$AMR_ENV" ]]; then
  export PATH="$AMR_ENV/bin:$PATH"
  export PERL5LIB="$AMR_ENV/lib/perl5/vendor_perl${PERL5LIB:+:$PERL5LIB}"
fi

QUERY="$ROOT/references/crpP_argannot_NG_062203_101_298.fasta"
ARGANNOT="${ARGANNOT_DB:-${AMR_ENV:+$AMR_ENV/db/argannot/sequences}}"
OUTDIR="$ROOT/results/blast/crpP"
DBDIR="$OUTDIR/db"
SUMMARY="$OUTDIR/BC42_BC88_crpP_summary.tsv"

mkdir -p "$OUTDIR" "$DBDIR" "$ROOT/logs/blast"

if [[ ! -s "$QUERY" ]]; then
  if [[ ! -s "$ARGANNOT" ]]; then
    echo "ERROR: ARG-ANNOT database not found. Set ARGANNOT_DB or provide $QUERY." >&2
    exit 1
  fi
  awk '
    /^>/ { keep = ($0 ~ /\(Flq\)crpP/) }
    keep { print }
  ' "$ARGANNOT" > "$QUERY"
fi

printf "sample\tcall\tquery\tsubject\tpident\tcoverage_pct\talignment_length\tquery_length\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore\n" > "$SUMMARY"

for SAMPLE in BC42 BC88; do
  ASM="$ROOT/results/polished/medaka/$SAMPLE/consensus.fasta"
  DB="$DBDIR/$SAMPLE"
  HITS="$OUTDIR/${SAMPLE}_crpP_blast.tsv"

  if [[ ! -s "$ASM" ]]; then
    echo "ERROR: assembly not found for $SAMPLE: $ASM" >&2
    exit 1
  fi

  makeblastdb -in "$ASM" -dbtype nucl -out "$DB" > "$ROOT/logs/blast/${SAMPLE}_crpP_makeblastdb.log" 2>&1

  blastn \
    -query "$QUERY" \
    -db "$DB" \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
    -perc_identity 70 \
    -qcov_hsp_perc 50 \
    -max_target_seqs 10 \
    -num_threads "$THREADS" \
    -out "$HITS" \
    2> "$ROOT/logs/blast/${SAMPLE}_crpP_blastn.log"

  sort -k11,11nr "$HITS" | awk -v sample="$SAMPLE" '
    BEGIN { OFS = "\t" }
    NR == 1 {
      cov = 100 * $4 / $5
      call = ($3 >= 80 && cov >= 80) ? "present_80id_80cov" : "below_80id_or_80cov"
      printf "%s\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", sample, call, $1, $2, $3, cov, $4, $5, $6, $7, $8, $9, $10, $11
      found = 1
    }
    END {
      if (!found) {
        printf "%s\tno_hit\tNA\tNA\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", sample
      }
    }
  ' >> "$SUMMARY"
done

echo "crpP BLAST summary: $SUMMARY"
