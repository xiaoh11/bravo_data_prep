# Prepare VCF workflow

## Data Requirements
This pipeline requires both:
- full BCFs (those with DP & GQ fields)
- site only BCFs (those with qc_metrics specified in config but no genotypes)

## Preparation 
- Symlinks to reference: taken care of by ansible deployment.
- Symlinks to bcfs: symlink to full and sites directories.

The symlinks to bcfs and crams still done manually or by some script in `hacks/`

## Running
Interactively on  SLURM:
```sh
nextflow run vcf.nf -with-trace -profile slurm
```

Run in background on SLURM:
```sh
nextflow run vcf.nf -with-trace -profile slurm -ansi-log false -bg > coverage.log
```
For shared clusters, it would be better practice to have a wrapper to submit the nextflow run as a slurm job.
That would avoid putting load on the login node.
Since this workflow is written for a single purpose cluster, it isn't an issue.

## Gotchas
### Apriori determine *all* your samples ids are present.
If there are any sample ids in the sample list that are not present in the bcf files, ComputeAllelesAndHistograms will exit with error.

```sh
INPUT_BCF=/path/to/input.bcf
SAMPLES_FILE=/path/to/your/samples.txt

# Get sample ids from a bcf
BCF_SAMPLES_FILE=/path/where/ids/output.txt
bcftools query -l ${INPUT_BCF} > ${BCF_SAMPLES_FILE}

# Find which of your sample ids aren't in the BCF's sample ids
grep -F -v -f ${BCF_SAMPLE_IDS} ${SAMPLES_FILE}
```
