# AnalyzeSNP

This repository contains utilities for analyzing sequencing reads and extracting base substitution frequencies.

## Requirements

The scripts assume the following command line tools are installed and in your `PATH`:

- `bowtie2`
- `samtools`
- `varscan`
- `trim_galore`
- `pacbam`
- Python 3 with `pandas`

## Quick start

1. Place FASTQ files under a `data/` directory and a reference FASTA under `ref/`.
2. Run `process.sh` to align reads with Bowtie2, generate mpileup files and call SNPs with VarScan.
3. Run `process_pacbam.sh` to generate PacBio style pileups and calculate mutation frequencies.
4. Output files are written to `output/`.

## Scripts

- `calculate_mutation.py` – compute base substitution rates from an mpileup file.
- `parse_readcounts.py` – parse VarScan `readcounts` output and generate mutation matrices.
- `readcounts.py` – compute mutation frequencies from PacBAM pileup output.
- `plot_coverage.py` – show a simple text‑based coverage plot for an mpileup file.
- `process.sh` – an end‑to‑end pipeline using Bowtie2 and VarScan.
- `process_pacbam.sh` – similar pipeline using PacBAM.

See comments inside each script for more options and details.
