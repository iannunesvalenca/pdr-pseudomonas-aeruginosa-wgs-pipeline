#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="${PROJECT_ROOT}/.tools"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
MICROMAMBA="${TOOLS_DIR}/bin/micromamba"
ENV_FILE="${PROJECT_ROOT}/envs/pdr-pa-wgs.yml"

if [[ "${MAMBA_ROOT_PREFIX}" == *" "* ]]; then
  echo "MAMBA_ROOT_PREFIX contains spaces: ${MAMBA_ROOT_PREFIX}" >&2
  echo "Conda-style environments in paths with spaces break tools with shebang scripts." >&2
  echo "Set MAMBA_ROOT_PREFIX to a path without spaces, for example: ${HOME}/micromamba-pdr-pa-wgs" >&2
  exit 1
fi

mkdir -p "${TOOLS_DIR}" "${MAMBA_ROOT_PREFIX}" "${PROJECT_ROOT}/logs"

if [[ ! -x "${MICROMAMBA}" ]]; then
  echo "Downloading micromamba into ${TOOLS_DIR}"
  wget -qO "${TOOLS_DIR}/micromamba.tar.bz2" "https://micro.mamba.pm/api/micromamba/linux-64/latest"
  tar -xjf "${TOOLS_DIR}/micromamba.tar.bz2" -C "${TOOLS_DIR}" bin/micromamba
fi

echo "Creating/updating environment from ${ENV_FILE}"
"${MICROMAMBA}" env create \
  --yes \
  --root-prefix "${MAMBA_ROOT_PREFIX}" \
  --file "${ENV_FILE}" \
  2>&1 | tee "${PROJECT_ROOT}/logs/bootstrap_env.log"

cat <<EOF

Environment ready.

Activate it with:
  eval "\$(${MICROMAMBA} shell hook --shell bash --root-prefix ${MAMBA_ROOT_PREFIX})"
  micromamba activate pdr-pa-wgs

Quick check:
  snakemake --version
  flye --version
  medaka --version
EOF
