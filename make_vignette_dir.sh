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

# Default to Reference locations on cluster for populationg ref data
BASIS_SRC_DIR="${BASIS_SRC_DIR:-/apps/reference/api}"
REF_SRC_DIR="${REF_SRC_DIR:-/apps/reference}"

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

# Create default API cache dir
mkdir -p "${RUNTIME}"/cache

mkdir -p "${RUNTIME}"/reference

# If REF_SRC_DIR is present, copy reference data under runtime directory
if [ -d ${REF_SRC_DIR} ]; then
  cp "${REF_SRC_DIR}"/hs38DH.fa     "${RUNTIME}"/reference/hs38DH.fa
  cp "${REF_SRC_DIR}"/hs38DH.fa.fai "${RUNTIME}"/reference/hs38DH.fa.fai
else
# Create placeholder for hs38DH reference
  touch "${RUNTIME}"/reference/hs38DH.fa.GOES_HERE
  touch "${RUNTIME}"/reference/hs38DH.fa.fai.GOES_HERE
fi

###################
# Fill Basis Data #
###################

mkdir -p "${BASIS}"/vcfs
cp -Lru "${RES_PREP_VCF}"/final/vcfs/* "${BASIS}"/vcfs

mkdir -p "${BASIS}"/qc_metrics
cp -Lru "${RES_PREP_VCF}"/final/qc_metrics/metrics.json.gz "${BASIS}"/qc_metrics/metrics.json.gz

mkdir -p "${BASIS}"/reference

# If BASIS_SRC_DIR is present, copy reference data under basis directory
if [ -d ${BASIS_SRC_DIR} ]; then
  cp "${BASIS_SRC_DIR}/canonical_transcripts.tsv.gz" "${BASIS}/reference/canonical_transcripts.tsv.gz"
  cp "${BASIS_SRC_DIR}/hgcn_genenames.tsv.gz"        "${BASIS}/reference/hgcn_genenames.tsv.gz"
  cp "${BASIS_SRC_DIR}/omim_ensembl_refs.tsv.gz"     "${BASIS}/reference/omim_ensembl_refs.tsv.gz"
  cp "${BASIS_SRC_DIR}/gencode.v38.annotation.gtf.gz"    "${BASIS}/reference/gencode.v38.annotation.gtf.gz"
else
# Otherwise Create placeholders for external data
  touch "${BASIS}/reference/canonical_transcripts.tsv.gz.GOES_HERE"
  touch "${BASIS}/reference/hgcn_genenames.tsv.gz.GOES_HERE"
  touch "${BASIS}/reference/omim_ensembl_refs.tsv.gz.GOES_HERE"
  touch "${BASIS}/reference/gencode_annotation.gtf.gz.GOES_HERE"
fi
