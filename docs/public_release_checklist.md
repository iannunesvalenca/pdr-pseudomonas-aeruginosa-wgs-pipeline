# Public Release Checklist

Before pushing this repository publicly, verify that it does not contain:

- FASTQ/FAST5/POD5 reads
- BAM/CRAM/SAM alignments
- assembled FASTA/GenBank/GFF files
- local databases
- draft manuscripts or reviewer reports
- patient identifiers or clinical source documents
- generated result tables that should remain private until publication

Useful checks:

```bash
git status --short
git ls-files
find . -type f -size +5M -print
rg -n "patient|manuscript|docx|fastq|bam|barcode|/home/ian|sharing" .
```
