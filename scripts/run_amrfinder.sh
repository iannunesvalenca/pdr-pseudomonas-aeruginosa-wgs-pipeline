#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
export MAMBA_ROOT_PREFIX
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"

cd "${PROJECT_ROOT}"
mkdir -p results/amr/amrfinder logs/amr

run_amr() {
  "${MICROMAMBA}" run -n pdr-pa-amr bash -lc "$*"
}

for sample in BC42 BC88; do
  assembly="results/polished/medaka/${sample}/consensus.fasta"
  run_amr "amrfinder \
    --nucleotide ${assembly} \
    --organism Pseudomonas_aeruginosa \
    --plus \
    --threads 16 \
    --name ${sample} \
    --mutation_all results/amr/amrfinder/${sample}_all_mutations.tsv \
    --output results/amr/amrfinder/${sample}_amrfinder.tsv \
    --log logs/amr/${sample}_amrfinder.log"
done

awk -F '\t' '
BEGIN {
  OFS = "\t"
  print "sample_id","contig","start","stop","strand","symbol","scope","type","subtype","class","subclass","method","coverage_pct","identity_pct"
}
FNR > 1 && ($10 == "AMR" || $11 ~ /POINT/) {
  print $1,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14,$17,$18
}' results/amr/amrfinder/BC42_amrfinder.tsv results/amr/amrfinder/BC88_amrfinder.tsv \
  > results/amr/amrfinder/BC42_BC88_amrfinder_focused.tsv

for sample in BC42 BC88; do
  awk -F '\t' 'NR > 1 {print $7}' "results/amr/amrfinder/${sample}_amrfinder.tsv" | sort -u \
    > "results/amr/amrfinder/${sample}_symbols.txt"
done

{
  echo "Only BC42"
  comm -23 results/amr/amrfinder/BC42_symbols.txt results/amr/amrfinder/BC88_symbols.txt
  echo "Only BC88"
  comm -13 results/amr/amrfinder/BC42_symbols.txt results/amr/amrfinder/BC88_symbols.txt
  echo "Shared"
  comm -12 results/amr/amrfinder/BC42_symbols.txt results/amr/amrfinder/BC88_symbols.txt
} > results/amr/amrfinder/BC42_BC88_symbol_comparison.txt

echo "Wrote:"
echo "  results/amr/amrfinder/BC42_amrfinder.tsv"
echo "  results/amr/amrfinder/BC88_amrfinder.tsv"
echo "  results/amr/amrfinder/BC42_BC88_amrfinder_focused.tsv"
echo "  results/amr/amrfinder/BC42_BC88_symbol_comparison.txt"
