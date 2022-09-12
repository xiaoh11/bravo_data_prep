# Coverage Workflow
Calculate the read depth for various base pair widths (bin widths)

## Requires
Following executeables need to be findable in the PATH.  
On the slurm cluster, these should be provided by modules. 

- samtools
- bam
- bgzip
- tabix
- aggregate.py
- prune.py

## Cram File Inputs
The cram files used in the pipeline should be symlinked under the data/crams directory.
This facilitates aggregating crams from multiple sources or projects.
It avoids having multiple or complex globs in the nextflow.config.

## Running

Interactively on  SLURM:
```sh
nextflow run depth_statistics.nf -with-trace -profile slurm
```

Run in background on SLURM:
```sh
nextflow run depth_statistics.nf -with-trace -profile slurm -ansi-log false -bg > coverage.log
```
For shared clusters, it would be better practice to have a wrapper to submit the nextflow run as a slurm job.
That would avoid putting load on the login node.
Since this is a single purpose cluster, isn't an issue.
