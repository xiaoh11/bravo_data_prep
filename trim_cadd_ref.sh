#!/usr/bin/env bash

# Trim reference cadd tsv file to position range.
#  Assumes the cadd.tsv.gz file is single chromosome per file and in order.
#
# Use: ./trim_cadd_ref.sh START END IN_FILE > OUT_FILE
#
# The output file will subsequently have to be tabix indexed.
#  tabix -b 2 -e 2 OUT_FILE

START_POS=$1
END_POS=$2
IN_FILE=$3

zcat "${IN_FILE}" |\
  awk -F "\t" -v lower=${START_POS} -v upper=${END_POS} \
    'NR < 3 {print}; $2 >= lower && $2 <= upper {print}' |\
  bgzip

