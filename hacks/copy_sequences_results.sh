#!/usr/bin/env bash

# Copy result from sequences workflow to results GCP bucket

SRC=/home/grosscol_umich_edu/data_prep/workflows/sequences/result
RUN_NAME=fr10_2chrom
DEST=/mnt/results/${RUN_NAME}/sequences

if [ ! -d ${SRC} ]; then
  echo "Dir not found: ${SRC}"
  exit 1
fi

if [ ! -d ${DEST} ]; then
  mkdir -p ${DEST}
fi

# Wrap sequences result into a tarball
echo "Tar'ing sequences"
cd ${SRC}
# tar -chf sequences_final.tar merged_variant_map.* sequences

# Copy the full coverage result
echo "Copying sequences final results tarball into bucket" 
gsutil cp sequences_final.tar gs://bravo-results/${RUN_NAME}/sequences/sequences_final.tar

echo "Done"
