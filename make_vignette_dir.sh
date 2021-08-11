#!/usr/bin/env bash

# Create the runtime and basis data directories from completed workflow results.
#   The data prep workflows are not run by this script.
#
# Takes path to top level directory as an arg.
#   The basis and runtime dirs will be created and populated under that dir.
#
# Use:
# ./make_vignette_dir.sh /path/to/top/dir
#
#
# Basis data supports the Bravo API custom commands to populate mongo DB:
#  flask load-genes \
#          "${BASIS}/reference/canonical_transcripts.tsv.gz" \
#          "${BASIS}/reference/omim_ensembl_refs.tsv.gz" \
#          "${BASIS}/reference/hgcn_genenames.tsv.gz"\
#          "${BASIS}/reference/gencode_annotation.gtf.gz"
#
# flask load-snv 2 "${BASIS}"/vcfs/*.vcf.gz
#
# flask load-qc-metrics "${BASIS}"/qc_metrics/metrics.json.gz

# Default to making a data directory in working dir.
TOP_DIR="${1:-./data}"

# Assume workflows directory is in working dir.
WORKFLOWS="./workflows"

# Result directories
RES_COVERAGE="${WORKFLOWS}"/coverage/result 
RES_SEQUENCE="${WORKFLOWS}"/sequences/result
RES_PREP_VCF="${WORKFLOWS}"/prepare_vcf_teddy/result

# Create directory structure.
RUNTIME="${TOP_DIR}"/runtime
BASIS="${TOP_DIR}"/basis

mkdir -p "${RUNTIME}"
mkdir -p "${BASIS}"

#####################
# Fill Runtime Data #
#####################

mkdir -p "${RUNTIME}"/coverage
cp -Lru "${RES_COVERAGE}"/bin* "${RUNTIME}"/coverage/
cp -Lru "${RES_COVERAGE}"/full "${RUNTIME}"/coverage/

mkdir -p "${RUNTIME}"/crams
cp -Lru "${RES_SEQUENCE}"/sequences "${RUNTIME}"/crams/
cp -Lru "${RES_SEQUENCE}"/merged_variant_map.tsv.gz "${RUNTIME}"/crams/variant_map.tsv.gz
cp -Lru "${RES_SEQUENCE}"/merged_variant_map.tsv.gz.tbi "${RUNTIME}"/crams/variant_map.tsv.gz.tbi

# Create placeholder for external data
#  E.g hs38DH.fa and hs38DH.fa.fai
mkdir -p "${RUNTIME}"/reference
touch "${RUNTIME}"/reference/REF.fa.GOES_HERE
touch "${RUNTIME}"/reference/REF.fa.fai.GOES_HERE

# Create default API cache dir
mkdir -p "${RUNTIME}"/cache

###################
# Fill Basis Data #
###################

mkdir -p "${BASIS}"/vcfs
cp -Lru "${RES_PREP_VCF}"/final/vcfs/* "${BASIS}"/vcfs

mkdir -p "${BASIS}"/qc_metrics
cp -Lru "${RES_PREP_VCF}"/final/qc_metrics/metrics.json.gz "${BASIS}"/qc_metrics/metrics.json.gz

# Create placeholders for external data
mkdir -p "${BASIS}"/reference
touch "${BASIS}/reference/canonical_transcripts.tsv.gz.GOES_HERE"
touch "${BASIS}/reference/hgcn_genenames.tsv.gz.GOES_HERE"
touch "${BASIS}/reference/omim_ensembl_refs.tsv.gz.GOES_HERE"
touch "${BASIS}/reference/gencode_annotation.gtf.gz.GOES_HERE"
