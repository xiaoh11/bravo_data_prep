#!/usr/bin/env bash

# Copy result from coverage workflow to results GCP bucket

SRC=/home/grosscol_umich_edu/data_prep/workflows/coverage/result
RUN_NAME=fr10_2chrom
DEST=/mnt/results/${RUN_NAME}/coverage

if [ ! -d ${SRC} ]; then
  echo "Dir not found: ${SRC}"
  exit 1
fi

if [ ! -d ${DEST} ]; then
  mkdir -p ${DEST}
fi

# Copy the full coverage result
echo "Copying full coverage result" 
cp -Lr ${SRC}/full ${DEST}

# Copy the pruned coverage result
echo "Copying full coverage result" 
cp -Lr ${SRC}/bin_* ${DEST}
