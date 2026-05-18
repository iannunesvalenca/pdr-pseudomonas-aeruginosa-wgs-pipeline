#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
export MAMBA_ROOT_PREFIX
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${PROJECT_ROOT}/.cache/matplotlib}"

cd "${PROJECT_ROOT}"
mkdir -p results/assembly results/qc/assembly logs/qc "${MPLCONFIGDIR}"

run_wgs() {
  "${MICROMAMBA}" run -n pdr-pa-wgs bash -lc "$*"
}

assemblies=(
  results/assembly/flye/BC42/assembly.fasta
  results/assembly/flye/BC88/assembly.fasta
  results/polished/medaka/BC42/consensus.fasta
  results/polished/medaka/BC88/consensus.fasta
)

run_wgs "seqkit stats -a -T ${assemblies[*]}" > results/assembly/assembly_seqkit_stats.tsv

run_wgs "quast.py --threads 16 --min-contig 1000 -o results/qc/assembly/quast_medaka results/polished/medaka/BC42/consensus.fasta results/polished/medaka/BC88/consensus.fasta" \
  2>&1 | tee logs/qc/quast_medaka.log

echo "Wrote:"
echo "  results/assembly/assembly_seqkit_stats.tsv"
echo "  results/qc/assembly/quast_medaka/report.tsv"

