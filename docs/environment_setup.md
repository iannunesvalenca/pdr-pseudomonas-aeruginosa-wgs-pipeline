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

Additional tool-specific environments used by optional/downstream scripts:

```bash
# MLST typing (scripts/run_mlst.sh)
micromamba create -y -n pdr-pa-typing \
  -c conda-forge -c bioconda --strict-channel-priority \
  mlst

# IntegronFinder (scripts/run_integronfinder.sh)
micromamba create -y -n pdr-pa-mobile \
  -c conda-forge -c bioconda --strict-channel-priority \
  integron_finder

# Bakta genome annotation (scripts/run_bakta_annotation.sh)
micromamba create -y -n pdr-pa-bakta \
  -c conda-forge -c bioconda --strict-channel-priority \
  bakta

# Core SNP comparison (scripts/run_snippy_ctgs.sh)
micromamba create -y -n pdr-pa-snippy \
  -c conda-forge -c bioconda --strict-channel-priority \
  snippy snp-dists
```

You can also bootstrap `pdr-pa-wgs` directly from `envs/pdr-pa-wgs.yml`:

```bash
scripts/bootstrap_env.sh
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
