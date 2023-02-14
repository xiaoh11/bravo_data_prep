#!/usr/bin/env bash

# Copy result from process vcf workflow to results GCP bucket

SRC=/home/grosscol_umich_edu/data_prep/workflows/process_vcf/result
RUN_NAME=fr10_2chrom
DEST=/mnt/results/${RUN_NAME}/process_vcf

if [ ! -d ${SRC} ]; then
  echo "Dir not found: ${SRC}"
  exit 1
fi

if [ ! -d ${DEST} ]; then
  mkdir -p ${DEST}
fi

# Copy the full coverage result
echo "Copying qc_metrics" 
cp -Lr ${SRC}/final/qc_metrics ${DEST}

# Copy the de-identified sequences
echo "Copying vcfs" 
cp -Lr ${SRC}/final/vcfs ${DEST}
