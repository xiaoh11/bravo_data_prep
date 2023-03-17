#!/usr/bin/env bash

# Copy result from process vcf workflow to results GCP bucket

SRC=/home/grosscol_umich_edu/data_prep/workflows/process_vcf/result
RUN_NAME=freeze10
BUCKET=bravo-results
DEST=gs://${BUCKET}/${RUN_NAME}/process_vcf

# Requires gsutil to be installed
command -v gsutil --version 1> /dev/null 2>&1 || \
{ echo >&2 "gsutil required but it's not installed.  Aborting."; exit 1; }

if [ ! -d ${SRC} ]; then
  echo "Dir not found: ${SRC}"
  exit 1
fi

if ! gsutil ls "gs://${BUCKET}" 1> /dev/null 2>&1; then
  echo "Bucket not found: gs://${BUCKET}"
  exit 1
fi

# Copy the full coverage result
echo "Copying qc_metrics" 
gsutil -m cp -r ${SRC}/final/qc_metrics ${DEST}

# Copy the de-identified sequences
echo "Copying vcfs" 
gsutil -m cp -r ${SRC}/final/vcfs ${DEST}
