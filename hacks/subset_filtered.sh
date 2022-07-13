#!/usr/bin/env bash

OUT_DIR=~/bravo-data-subset

SRC_BCF=/mnt/vcfs/freeze10/filtered_genotypes/minDP0/merged.chr11.gtonly.minDP0.filtered.bcf

OUT_BCF=${OUT_DIR}/chr11_filtered_subset.bcf

# Subset the bcf file
bcftools view ${SRC_BCF}\
  -o ${OUT_BCF} \
  --output-type b \
  -r chr11:5220000-5520000 \
  --samples-file sample_ids.txt \
  --threads 1

# Index the subset file
bcftools index ${OUT_BCF} --threads 1
