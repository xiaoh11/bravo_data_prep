#!/usr/bin/env bash

# Subset full bcf to create testing subset of chr11
# Use baked in paths from GCP cluster provisioning and freeze 10 bucket.

OUT_DIR=~/bravo-data-subset

SRC_BCF=/mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_5200001_5300000.genotypes.bcf
OUT_BCF=${OUT_DIR}/chr11.subset.full.bcf
INTERMEDIATE_BASE=${OUT_DIR}/chr11

# Make output directory
mkdir -p ${OUT_DIR}

ARR=(\
/mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_5200001_5300000.genotypes.bcf \
/mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_5300001_5400000.genotypes.bcf \
/mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_5400001_5500000.genotypes.bcf \
/mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_5500001_5600000.genotypes.bcf )

# Trim BCFS to samples and region
COUNTER=0
INTERMEDIATES=()
for INFILE in "${ARR[@]}"
do
  echo "${INFILE}"

  OUTPATH="${INTERMEDIATE_BASE}.part${COUNTER}.bcf"

  # Subset the bcf file
  bcftools view ${INFILE}\
          -o ${OUTPATH} \
          --output-type b \
          -r chr11:5220000-5520000 \
          --samples-file sample_ids.txt \
          --threads 1

  INTERMEDIATES[${COUNTER}]=${OUTPATH}
  COUNTER=$((COUNTER+1))
done

# Consolidate intermediates
bcftools concat -n -o ${OUT_BCF} --output-type b \
  "${INTERMEDIATE_BASE}.part0.bcf" \
  "${INTERMEDIATE_BASE}.part1.bcf" \
  "${INTERMEDIATE_BASE}.part2.bcf" \
  "${INTERMEDIATE_BASE}.part3.bcf" 

# Index the subset file
bcftools index ${OUT_BCF} --threads 1
