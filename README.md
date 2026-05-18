# PDR Pseudomonas aeruginosa Nanopore WGS Pipeline

Public-ready analysis pipeline for paired Oxford Nanopore bacterial whole-genome
sequencing data from *Pseudomonas aeruginosa* isolates.

The repository contains the workflow scripts and documentation used for:

- Nanopore read quality control
- Flye assembly and Medaka polishing
- QUAST assembly QC
- MLST typing
- AMR gene and mutation screening with AMRFinderPlus and Abricate
- targeted `crpP` confirmation
- VFDB virulence profiling
- Bakta annotation
- IntegronFinder analysis
- candidate SNP annotation and read-backed validation
- manuscript-quality visualization scripts

## Data Availability and Privacy

No raw sequencing reads, BAM files, assemblies, databases, patient documents,
draft manuscripts, or generated analysis results are included in this public
repository. The scripts expect users to provide their own input files locally.

Large or sensitive files should remain outside git or in ignored directories
such as `data/`, `results/`, `logs/`, `references/`, and `databases/`.

## Repository Layout

```text
scripts/              Pipeline and visualization scripts
docs/                 Setup and workflow documentation
config/               Example sample and path configuration files
envs/                 Environment notes
```

## Quick Start

1. Install micromamba.
2. Create the environments described in `docs/environment_setup.md`.
3. Place input FASTQ files in a local, ignored `data/` directory.
4. Adapt `config/samples.tsv.example` for the local sample names and paths.
5. Run the scripts in the order described in `docs/workflow.md`.

## Important Notes

The scripts were developed for a two-isolate comparative analysis and currently
use the sample labels `BC42` and `BC88` in several places. Treat these labels as
template values and adjust them for new datasets before reuse.

The repository is intended to make the analysis approach reproducible. It is not
a data release.
