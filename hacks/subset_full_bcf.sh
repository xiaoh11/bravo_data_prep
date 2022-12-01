#!/usr/bin/env bash

# Subset full bcf to create testing subset of chr11
# Use baked in paths from GCP cluster provisioning and freeze 10 bucket.

OUT_DIR=~/subset/sampled
INFILE_LIST=chr11.list.txt

# Make output directory
mkdir -p ${OUT_DIR}

# Read list of bcf paths 
# Example contents:
# /mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_2600001_2700000.genotypes.bcf
# /mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_2700001_2800000.genotypes.bcf
# /mnt/vcfs/freeze10/genotypes/merged/chr11/merged.chr11_2800001_2900000.genotypes.bcf

readarray -t ARR < ${INFILE_LIST}

# Trim BCFS to samples and region
for INFILE in "${ARR[@]}"
do
  BNAME=$(basename ${INFILE})
  OUTPATH="${OUT_DIR}/${BNAME}"

  # Subset the bcf file
  bcftools view ${INFILE}\
          -o ${OUTPATH} \
          --output-type b \
          -r chr11:5220000-5520000 \
          --samples-file sample_ids.txt \
          --threads 1

  bcftools index ${OUTPATH} --threads 1
done
