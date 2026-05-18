# Workflow

This document summarizes the intended order of execution. Commands assume the
repository root as the working directory and that the relevant micromamba
environment is active.

## 1. Environment setup

See `docs/environment_setup.md`.

```bash
scripts/check_envs.sh
```

## 2. Read QC

Inputs: filtered Nanopore FASTQ files.

```bash
scripts/run_read_qc.sh
```

## 3. Assembly and polishing

Use Flye followed by Medaka. The exact Flye/Medaka commands are documented in
the project reports and can be adapted for each dataset.

## 4. Assembly QC

```bash
scripts/run_assembly_qc.sh
```

## 5. Typing

```bash
scripts/run_mlst.sh
```

## 6. AMR screening

```bash
scripts/run_amrfinder.sh
scripts/run_abricate_amr.sh
scripts/check_crpp.sh
```

## 7. Annotation and mobile elements

```bash
scripts/run_bakta_annotation.sh
scripts/extract_bakta_focus.sh
scripts/run_integronfinder.sh
```

## 8. Virulence screening

```bash
scripts/run_vfdb.sh
```

## 9. SNP analysis

```bash
scripts/run_snippy_ctgs.sh
scripts/annotate_snps.py
```

Candidate SNPs from assembly comparison should be validated against the original
read alignments before being reported as fixed isolate-to-isolate differences.

## 10. Visualization

```bash
scripts/plot_genome_synteny.py
scripts/plot_integron_loci.py
scripts/plot_circular_genome_parallel.py
scripts/plot_snp_allele_support.py
scripts/plot_snp_coverage.py
scripts/plot_amr_landscape.py
scripts/plot_virulence_profile.py
```

Generated figures are intentionally ignored by git.
