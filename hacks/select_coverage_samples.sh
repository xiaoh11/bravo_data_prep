#!/usr/bin/env bash

# This is a preparation script for running the coverage estimation workflow. It handles putting the
#  cram and crai files in a convenient place for the workflow.  This is accomplished with symlinks 
#  in the data/ directory of the coverage workflow.
#
# - Collect cram object paths in crams buckets which have corresponding file mounts.
# - Create mapping tsv file for crams with id and mounted path columns. 
# - Create symlinks to a random subset of cram and crai in the coverage data directory.

##########
# INPUTS #
##########
SAMPLE_FILE=~/freeze10_corrected_sample_ids.txt
COVERAGE_CRAM_DIR=~/data_prep/workflows/coverage/data/crams

# Buckets where crams are located
CRAM_BUCKETS=("fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e"
              "ccdg-east"
              "ccdg-firecloud-crams-tmp"
              "topmed-irc-share")

# Corresponding mounts for buckets
CRAM_DIRS=("/mnt/crams_1000g" 
           "/mnt/crams_ccdg" 
					 "/mnt/crams_ccdg_tmp"
					 "/mnt/crams_topmed")

###########
# OUTPUTS #
###########
ALL_CRAMS_INDEX=/tmp/all_cram_paths.txt
ALL_CRAMS_MAP=/tmp/all_crams_map.txt

##############
# PROCESSING #
##############
# Sanity check that both arrays are same length
if [ ! "${#CRAM_BUCKETS[@]}" -eq "${#CRAM_DIRS[@]}" ]; then
	echo "Error: Arrays unequal length!"
	exit 1
fi

# Create indexes of crams in buckets
cat /dev/null > ${ALL_CRAMS_INDEX}
for I in "${!CRAM_BUCKETS[@]}"
do
	BUCKET=${CRAM_BUCKETS[$I]}
	DIR=${CRAM_DIRS[$I]}
	BUCKET_IDX="/tmp/${BUCKET}_idx.txt"
	CRAMS_IDX="/tmp/${BUCKET}_crams.txt"

	if [ ! -f "${BUCKET_IDX}" ]; then
		echo "writing ${BUCKET_IDX}"
		gsutil -u genome-variant-server ls -r "gs://${BUCKET}/**" > ${BUCKET_IDX}
	fi

	if [ ! -f "${CRAMS_IDX}" ]; then
		echo "writing ${CRAMS_IDX}"
		grep -e '\.cram$' "${BUCKET_IDX}" |\
      sed "s;gs://${BUCKET};${DIR};" > ${CRAMS_IDX}
	fi

	# Append to unsorted crams index
	#  Remove lines with public/ to avoid dupes
	grep -v -e 'public/' ${CRAMS_IDX} >> ${ALL_CRAMS_INDEX}
done


# Create mapping file of ids to path
if [ ! -f ${ALL_CRAMS_MAP} ]; then
	EXTRACT_ID_PATT='s;.*/\([[:alnum:]_-]*\).*;\1;'
	mkfifo /tmp/basenames.fifo
	(cat ${ALL_CRAMS_INDEX} | sed "${EXTRACT_ID_PATT}" > /tmp/basenames.fifo &)
	paste -d "\t" /tmp/basenames.fifo ${ALL_CRAMS_INDEX} | sort > ${ALL_CRAMS_MAP}
	rm /tmp/basenames.fifo
fi

# Create ALL the symlinks for both the cram and crai
OLD_IFS=$IFS
IFS=$'\t'
while read ID TARGET; do
	echo "$ID"
	echo "$TARGET"
	# ln -f -s "${TARGET}" "${COVERAGE_CRAM_DIR}/${ID}.cram"
	# ln -f -s "${TARGET}.crai" "${COVERAGE_CRAM_DIR}/${ID}.cram.crai"
done < <(shuf -n 10 ${ALL_CRAMS_MAP})
