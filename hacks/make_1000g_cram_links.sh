#!/usr/bin/env bash

# Script to create links to a 1000 Genome crams in convenience data directory
#  for testing the coverage workflow on slurm cluster.

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DEST_DIR=/home/grosscol_umich_edu/data_prep/workflows/coverage/data/crams

SRC13607=/mnt/crams/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06
SRC14151=/mnt/crams/CCDG_14151/Project_CCDG_14151_B01_GRM_WGS.cram.2020-02-12

mkdir -p ${DEST_DIR}

readarray -t ARR13607 < ${SCRIPT_DIR}/1000g_samples_ccdg13607.txt
readarray -t ARR14151 < ${SCRIPT_DIR}/1000g_samples_ccdg14151.txt

for ID in ${ARR13607[@]}
do
  CRAM=${SRC13607}/Sample_${ID}/analysis/${ID}.final.cram
  CRAI=${SRC13607}/Sample_${ID}/analysis/${ID}.final.cram.crai

  ln -sf ${CRAM} ${DEST_DIR}/${ID}.final.cram
  ln -sf ${CRAI} ${DEST_DIR}/${ID}.final.cram.crai
done

for ID in ${ARR14151[@]}
do
  CRAM=${SRC14151}/Sample_${ID}/analysis/${ID}.final.cram
  CRAI=${SRC14151}/Sample_${ID}/analysis/${ID}.final.cram.crai

  ln -sf ${CRAM} ${DEST_DIR}/${ID}.final.cram
  ln -sf ${CRAI} ${DEST_DIR}/${ID}.final.cram.crai
done
