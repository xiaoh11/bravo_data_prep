# One off scripts
Hacks to assist in local development and testing.
Dealing with moving directories of links or crams and trimming reference data to vignette size.

## Prepare VCFS
make_full_bcf_links.sh used to create links to a subset of the full vcfs.
This keeps the processing time down to 20-ish minutes
Used 1000g_samples.txt as the sample ids list.

## Coverage
make_1000g_cram_links.sh used to create links to the crams 

Uses id list text files:
1000g_samples_ccdg13607.txt
1000g_samples_ccdg14151.txt

Makes links in coverage/data/crams directory

