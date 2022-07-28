#!/usr/bin/env bash

# Script to create links to a subset of full bcfs in convenience data directory
#  for testing the process_vcf workflow on slurm cluster.

SRC_DIR=/mnt/vcfs/freeze10/genotypes/merged
DEST_DIR=/home/grosscol_umich_edu/data_prep/workflows/process_vcf/data/bcfs/subset_full

mkdir -p ${DEST_DIR}/chr11

ln -sf ${SRC_DIR}/chr11/merged.chr11_5100001_5200000.genotypes.bcf ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5100001_5200000.genotypes.bcf.csi ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5200001_5300000.genotypes.bcf ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5200001_5300000.genotypes.bcf.csi ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5300001_5400000.genotypes.bcf ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5300001_5400000.genotypes.bcf.csi ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5400001_5500000.genotypes.bcf ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5400001_5500000.genotypes.bcf.csi ${DEST_DIR}/chr11
ln -sf ${SRC_DIR}/chr11/merged.chr11_5500001_5600000.genotypes.bcf ${DEST_DIR}/chr11

 mkdir -p ${DEST_DIR}/chr12

ln -sf ${SRC_DIR}/chr12/merged.chr12_5100001_5200000.genotypes.bcf ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5100001_5200000.genotypes.bcf.csi ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5200001_5300000.genotypes.bcf ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5200001_5300000.genotypes.bcf.csi ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5300001_5400000.genotypes.bcf ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5300001_5400000.genotypes.bcf.csi ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5400001_5500000.genotypes.bcf ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5400001_5500000.genotypes.bcf.csi ${DEST_DIR}/chr12
ln -sf ${SRC_DIR}/chr12/merged.chr12_5500001_5600000.genotypes.bcf ${DEST_DIR}/chr12
