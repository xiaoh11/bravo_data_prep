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
