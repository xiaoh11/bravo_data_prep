#!/usr/bin/env bash

SRC=/home/grosscol_umich_edu/data_prep/workflows/coverage/result
RUN_NAME=freeze10
BUCKET=bravo-results
DEST=gs://${BUCKET}/${RUN_NAME}/coverage

# Requires gsutil to be installed
command -v gsutil 1> /dev/null 2>&1 || \
{ echo >&2 "gsutil required but it's not installed.  Aborting."; exit 1; }

# Verify source and destination
if [ ! -d ${SRC} ]; then
  echo "Dir not found: ${SRC}"
  exit 1
fi

if ! gsutil ls "gs://${BUCKET}" 1> /dev/null 2>&1; then
  echo "Bucket not found: gs://${BUCKET}"
  exit 1
fi

# Copy result from coverage workflow to results GCP bucket

# Copy the full coverage result
echo "Copying full coverage result" 
# gsutil -m cp -r ${SRC}/full ${DEST}
 
# Copy the pruned coverage result
echo "Copying full coverage result" 
gsutil -m cp -r "${SRC}/bin_*" ${DEST}
