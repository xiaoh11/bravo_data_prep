#!/usr/bin/env bash
# Script to create links to a 198 samples from 1000 Genome crams in convenience data directory
#  for testing the sequence workflow on slurm cluster.
# All 198 samples are from the CCDG 14151 project

SRC14151=/mnt/crams/CCDG_14151/Project_CCDG_14151_B01_GRM_WGS.cram.2020-02-12
DEST_DIR=/home/grosscol_umich_edu/data_prep/workflows/sequences/data/crams
HACKS_DIR=/home/grosscol_umich_edu/data_prep/hacks

mkdir -p ${DEST_DIR}
readarray -t ARR_198 < ${HACKS_DIR}/198_samples.txt

for ID in ${ARR_198[@]}
do
  CRAM=${SRC14151}/Sample_${ID}/analysis/${ID}.final.cram
  CRAI=${SRC14151}/Sample_${ID}/analysis/${ID}.final.cram.crai

  ln -sf ${CRAM} ${DEST_DIR}/${ID}.final.cram
  ln -sf ${CRAI} ${DEST_DIR}/${ID}.final.cram.crai
done
