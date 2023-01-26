#!/usr/bin/env bash

# Process the samples data (csv file with info about the included samples) to sample ids.
# Check that the user ids are present in the bcf prior to running process_vcf workflow.
# If there are any sample ids in the sample list that are not present in the bcf files,
# ComputeAllelesAndHistograms will exit with error.
# Output corrected sample ids to include only ids that are present in the bcf.

# Requires bcftools to be in PATH
command -v bcftools 1> /dev/null 2>&1 || \
  { echo >&2 "bcftools required but not found.  Try loading the bcftools module."; exit 1; }

###############
# INPUT FILES #
###############
BCF=/mnt/vcfs/freeze10/genotypes/merged/chr22/merged.chr22_9300001_9400000.genotypes.bcf
SAMPLES_FILE=/home/grosscol_umich_edu/BRAVO_Samples_Freeze10_20230117.txt

######################
# INTERMEDIATE FILES #
######################
SAMPLE_IDS_FILE=~/freeze10_raw_sample_ids.txt
BCF_SAMPLES_FILE=/tmp/bcf_ids.txt
SORTED_SAMPLES_FILE=/tmp/sorted_samples.txt
SAMPLE_IDS_NOT_IN_BCF_FILE=/tmp/sample_ids_NOT_IN_bcf.txt

################
# OUTPUT FILES #
################
CORRECTED_SAMPLES_FILE=~/freeze10_corrected_sample_ids.txt

################
# PROCESS DATA #
################

# Get sample ids from samples data. Use ccdg_id where present.
tail -n +2 ${SAMPLES_FILE} |\
  awk -F ',' '{ if ( $11 =="NA" ) {print $1} else {print $11} }' > ${SAMPLE_IDS_FILE}

# Get sample ids from a bcf
bcftools query -l ${BCF} | sort > ${BCF_SAMPLES_FILE}

echo -e "\n----------"
echo "lines in input samples file:"
wc -l ${SAMPLE_IDS_FILE}

echo -e "\n----------"
echo "lines bcf ids:"
wc -l ${BCF_SAMPLES_FILE}

# Find which of your sample ids aren't in the BCF's sample ids
sort ${SAMPLE_IDS_FILE} > ${SORTED_SAMPLES_FILE}
comm -13 "${BCF_SAMPLES_FILE}" "${SORTED_SAMPLES_FILE}" > ${SAMPLE_IDS_NOT_IN_BCF_FILE}

echo -e "\n----------"
echo "Number of extraneous ids:"
wc -l ${SAMPLE_IDS_NOT_IN_BCF_FILE}

echo -e "\n----------"
echo "List of extraneous recorded in file:"
echo ${SAMPLE_IDS_NOT_IN_BCF_FILE}

# Write corrected list of ids
comm -12 "${BCF_SAMPLES_FILE}" "${SORTED_SAMPLES_FILE}" > ${CORRECTED_SAMPLES_FILE}

echo -e "\n----------"
echo "Number of corrected ids:"
wc -l ${CORRECTED_SAMPLES_FILE}
echo ""
echo "Corrected IDs file written:"
echo ${CORRECTED_SAMPLES_FILE}

