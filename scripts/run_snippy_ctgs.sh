#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-16}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

REF="$ROOT/references/AP012280.1_NCGM2S1.fasta"
BC42_ASM="$ROOT/results/polished/medaka/BC42/consensus.fasta"
BC88_ASM="$ROOT/results/polished/medaka/BC88/consensus.fasta"

OUTDIR="$ROOT/results/snp/snippy_ctgs"
LOGDIR="$ROOT/logs/snp/snippy_ctgs"

mkdir -p "$OUTDIR" "$LOGDIR"

for cmd in snippy snippy-core snp-dists; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: '$cmd' not found. Activate the pdr-pa-snippy environment first." >&2
    exit 1
  fi
done

for path in "$REF" "$BC42_ASM" "$BC88_ASM"; do
  if [[ ! -s "$path" ]]; then
    echo "ERROR: Missing required input: $path" >&2
    exit 1
  fi
done

snippy \
  --cpus "$THREADS" \
  --reference "$REF" \
  --ctgs "$BC42_ASM" \
  --outdir "$OUTDIR/snippy_BC42_ST235" \
  --prefix "snippy_BC42_ST235" \
  2>&1 | tee "$LOGDIR/BC42_snippy_ctgs.log"

snippy \
  --cpus "$THREADS" \
  --reference "$REF" \
  --ctgs "$BC88_ASM" \
  --outdir "$OUTDIR/snippy_BC88_ST235" \
  --prefix "snippy_BC88_ST235" \
  2>&1 | tee "$LOGDIR/BC88_snippy_ctgs.log"

(
  cd "$OUTDIR"

  snippy-core \
    --ref "$REF" \
    --prefix core_ST235 \
    snippy_BC42_ST235 snippy_BC88_ST235 \
    2>&1 | tee "$LOGDIR/snippy_core_ST235.log"

  snp-dists core_ST235.aln > matrix_snps.tab

  awk -F'\t' '
    BEGIN { OFS = "\t" }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "snippy_BC42_ST235") bc42 = i
        if ($i == "snippy_BC88_ST235") bc88 = i
      }
      if (!bc42 || !bc88) {
        print "ERROR: Could not find sample columns in core_ST235.tab" > "/dev/stderr"
        exit 1
      }
      print
      next
    }
    {
      b1 = toupper($bc42)
      b2 = toupper($bc88)
      if (b1 != b2 && b1 ~ /^[ATCG]$/ && b2 ~ /^[ATCG]$/) print
    }
  ' core_ST235.tab > only_diff_snps.tsv

  {
    echo "Snippy ctgs SNP summary"
    echo "Date: $(date -Iseconds)"
    echo
    echo "Reference: $REF"
    echo "BC42 assembly: $BC42_ASM"
    echo "BC88 assembly: $BC88_ASM"
    echo
    echo "snp-dists matrix:"
    cat matrix_snps.tab
    echo
    echo "Differential core SNP rows in only_diff_snps.tsv:"
    awk 'NR > 1 { n++ } END { print n + 0 }' only_diff_snps.tsv
  } > snippy_ctgs_summary.txt
)

echo "Done. Summary: $OUTDIR/snippy_ctgs_summary.txt"
