# Representative Sequences Workflow

## Cram File Inputs
The cram files used in the pipeline should be symlinked under the data/crams directory.
This aggregates crams from multiple sources in once spot.
It avoids having multiple or complex configurations.

## VCF(BCF) File Inputs
Uses the filtered genotype only bcf files.

## Preparation 
- Make samples mapping tsv
- Symlink required crams

## Running

Interactively on  SLURM:
```sh
nextflow run Sequences.nf -with-report seq_report.html -profile slurm
```

Run in background on SLURM:
```sh
nextflow run Sequences.nf -with-report seq_report.html -profile slurm -ansi-log false -bg > coverage.log
