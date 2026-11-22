# PDR Pseudomonas aeruginosa Nanopore WGS Pipeline

Analysis pipeline for paired Oxford Nanopore whole-genome sequencing data from *Pseudomonas aeruginosa* isolates.

The repository includes scripts and documentation for:

- Nanopore read quality control
- Flye assembly and Medaka polishing
- QUAST assembly QC
- MLST typing
- AMR gene and mutation screening with AMRFinderPlus and Abricate
- Targeted `crpP` confirmation
- VFDB virulence profiling
- Bakta annotation
- IntegronFinder analysis
- Candidate SNP annotation and read-backed validation
- Plotting scripts for figures

## Data Availability and Privacy

No raw reads, BAM files, assemblies, databases, patient documents, draft manuscripts, or analysis outputs are stored in this repository. Provide your own input files locally.

Keep large or sensitive files outside git or inside ignored directories (`data/`, `results/`, `logs/`, `references/`, and `databases/`).

## Repository Layout

```text
scripts/              Pipeline and plotting scripts
docs/                 Setup and workflow guides
config/               Example sample and path configs
envs/                 Environment specs and notes
```

## Quick Start

1. Install micromamba.
2. Create the environments listed in `docs/environment_setup.md` (or run `scripts/bootstrap_env.sh` for the core `pdr-pa-wgs` environment).
3. Put input FASTQ files in a local `data/` directory.
4. Copy and edit `config/samples.tsv.example` for your sample names and paths.
5. Run the scripts in the order shown in `docs/workflow.md`.

## Notes

These scripts were written for a two-isolate comparison and use `BC42` and `BC88` as sample labels in several files. Treat those names as placeholders and update them for new datasets.

This repository shares the analysis workflow for reproducibility. It does not release sequencing data.

## License

MIT. See [`LICENSE`](LICENSE).
