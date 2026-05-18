#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
export MAMBA_ROOT_PREFIX
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"

cd "${PROJECT_ROOT}"
mkdir -p results/amr/abricate logs/amr

run_amr() {
  "${MICROMAMBA}" run -n pdr-pa-amr bash -lc "$*"
}

samples=(BC42 BC88)
dbs=(card ncbi resfinder)

for db in "${dbs[@]}"; do
  for sample in "${samples[@]}"; do
    assembly="results/polished/medaka/${sample}/consensus.fasta"
    run_amr "abricate --db ${db} --threads 4 ${assembly}" \
      > "results/amr/abricate/${sample}_${db}.tsv" \
      2> "logs/amr/${sample}_${db}.log"
  done

  run_amr "abricate --summary results/amr/abricate/BC42_${db}.tsv results/amr/abricate/BC88_${db}.tsv" \
    > "results/amr/abricate/summary_${db}.tsv"
done

awk -F '\t' '
BEGIN {
  OFS = "\t"
  print "sample_id","contig","start","end","strand","gene","coverage_pct","identity_pct","resistance"
}
FNR > 1 {
  sample = ($1 ~ /BC42/ ? "BC42" : "BC88")
  print sample,$2,$3,$4,$5,$6,$10,$11,$15
}' results/amr/abricate/BC42_ncbi.tsv results/amr/abricate/BC88_ncbi.tsv \
  > results/amr/abricate/BC42_BC88_ncbi_focused.tsv

echo "Wrote:"
echo "  results/amr/abricate/summary_card.tsv"
echo "  results/amr/abricate/summary_ncbi.tsv"
echo "  results/amr/abricate/summary_resfinder.tsv"
echo "  results/amr/abricate/BC42_BC88_ncbi_focused.tsv"
