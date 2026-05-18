# Environment Setup

This project uses micromamba environments installed outside the project path:

```text
MAMBA_ROOT_PREFIX=${HOME}/micromamba-pdr-pa-wgs
```

Keeping the environment outside the project avoids broken shebang scripts when the project path contains spaces.

## Environments

`pdr-pa-wgs` is the core environment for read QC, assembly, polishing, mapping, assembly QC, and first-pass typing:

```bash
micromamba create -y -n pdr-pa-wgs \
  -c conda-forge -c bioconda \
  --strict-channel-priority \
  python=3.10 \
  pandas biopython matplotlib seaborn openpyxl python-docx \
  seqkit nanoq nanoplot chopper filtlong \
  "samtools>=1.20" "htslib>=1.20" "bcftools>=1.20" \
  minimap2 flye raven-assembler racon medaka \
  quast fastani mash mlst blast
```

`pdr-pa-amr` is a separate AMR screening environment:

```bash
micromamba create -y -n pdr-pa-amr \
  -c conda-forge -c bioconda \
  --strict-channel-priority \
  abricate ncbi-amrfinderplus blast
```

## Activation

From the project directory:

```bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba-pdr-pa-wgs"
eval "$("$PWD/.tools/bin/micromamba" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate pdr-pa-wgs
```

## Verification

```bash
scripts/check_envs.sh
```
