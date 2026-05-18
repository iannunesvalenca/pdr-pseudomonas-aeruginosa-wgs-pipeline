#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
export MAMBA_ROOT_PREFIX
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"

cd "${PROJECT_ROOT}"
mkdir -p results/typing/mlst logs/typing

"${MICROMAMBA}" run -n pdr-pa-typing mlst \
  --scheme paeruginosa \
  results/polished/medaka/BC42/consensus.fasta \
  results/polished/medaka/BC88/consensus.fasta \
  > results/typing/mlst/BC42_BC88_mlst.tsv \
  2> logs/typing/BC42_BC88_mlst.log

echo "Wrote:"
echo "  results/typing/mlst/BC42_BC88_mlst.tsv"
