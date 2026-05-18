#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
export MAMBA_ROOT_PREFIX
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"

if [[ ! -x "${MICROMAMBA}" ]]; then
  echo "micromamba not found at ${MICROMAMBA}" >&2
  exit 1
fi

run_in_env() {
  local env_name="$1"
  shift
  "${MICROMAMBA}" run -n "${env_name}" bash -lc "$*"
}

echo "Mamba root: ${MAMBA_ROOT_PREFIX}"
echo

echo "[pdr-pa-wgs]"
run_in_env pdr-pa-wgs '
set -euo pipefail
python --version
seqkit version
nanoq --version
NanoPlot --version
chopper --version
filtlong --version 2>&1 | head -n 1
samtools --version | head -n 2
bcftools --version | head -n 2
minimap2 --version
flye --version
raven --version
racon --version 2>&1 | head -n 1
medaka --version
quast.py --version
fastANI --version 2>&1 | head -n 1
mash --version
if command -v mlst >/dev/null 2>&1; then
  mlst --version
else
  echo "mlst: not installed"
fi
'

echo
echo "[pdr-pa-amr]"
run_in_env pdr-pa-amr '
set -euo pipefail
abricate --version
amrfinder --version
blastn -version | head -n 2
'
