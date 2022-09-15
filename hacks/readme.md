# One off scripts
Hacks to assist in local development and testing.
Dealing with moving directories of links or crams and trimming reference data to vignette size.

## Processing Development 
- `mlr_aggregate.sh`: Demonstration script for using Miller as part of pipeline to accumulate depth counts.
    Subesequent work should involve using mosdepth instead of samtools pileup.
- `tiered_demo_aggregation.nf`: Illustrate using an unrolled loop of nextflow processes to aggregate many files a few at a time.
- `casting_prune.py`: Short term hack to workaround a data formatting mistake that would take a long time to re-process data.
- `subset_filtered.sh`: Used for making small subset of chr11 calls for local testing. (from genotype only bcfs)
- `subset_full_bcf.sh`: Used for making small subset of chr11 calls for local testing. (from full bcfs)
- `make_link_rewrite_script`: Used to rewrite links when the target files were moved.

## Small Scale Development & Testing
### Prepare VCFs
`make_full_bcf_links.sh` used to create links to a subset of the full vcfs.
This keeps the processing time down to 20-ish minutes
Used 1000g_samples.txt as the sample ids list.

### Coverage
`make_1000g_cram_links.sh` used to create links to the crams 

Uses id list text files:
1000g_samples_ccdg13607.txt
1000g_samples_ccdg14151.txt

Makes links in coverage/data/crams directory

## Mid Scale Development & Testing

## Prepare VCFs
Use 198_samples.txt as samples list.

## Sequences
Use `198_samples.tsv` as samples list.  Same IDs as `.txt` file of same name, but with additonal fields for sequences process.

