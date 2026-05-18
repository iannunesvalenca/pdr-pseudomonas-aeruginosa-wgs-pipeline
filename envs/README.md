# Software Environments

The analysis used multiple micromamba environments to avoid dependency
conflicts:

- `pdr-pa-wgs`: core QC, assembly, mapping, SNP work, and plotting
- `pdr-pa-amr`: AMRFinderPlus, Abricate, BLAST
- `pdr-pa-typing`: MLST typing
- `pdr-pa-mobile`: IntegronFinder
- `pdr-pa-bakta`: Bakta annotation

See `docs/environment_setup.md` for the concrete commands used during the
analysis.
