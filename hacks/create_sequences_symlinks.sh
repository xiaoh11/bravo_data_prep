#!/usr/bin/env bash

# This is a preparation script for running the sequences processing workflow. It handles putting the
#  cram and crai files in a convenient place for the workflow.  This is accomplished with symlinks 
#  in the data/ directory of the sequences workflow.
# This is a direct re-use of the coverage workflow setup script: select_coverage_samples.sh 

##########
# INPUTS #
##########
SAMPLE_FILE=~/freeze10_corrected_sample_ids.txt
SEQUENCES_CRAM_DIR=~/data_prep/workflows/sequences/data/crams
RND_ID_CHARS=10

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
ALL_CRAIS_INDEX=/tmp/all_crai_paths.txt
ALL_CRAMS_MAP=/tmp/all_crams_map.txt
ALL_CRAMS_IDS=/tmp/all_crams_ids.txt
RND_CRAMS_IDS=/tmp/rnd_crams_ids.txt

TRIM_SAMPLES_MAP=/tmp/trim_crams_map.txt
TRIM_CRAMS_REMAP=~/sequences_samples_map.tsv

##############
# PROCESSING #
##############
# Sanity check that both arrays are same length
if [ ! "${#CRAM_BUCKETS[@]}" -eq "${#CRAM_DIRS[@]}" ]; then
	echo "Error: Arrays unequal length!"
	exit 1
fi

# Create indexes of crams in buckets
if [ ! -f "${ALL_CRAMS_INDEX}" ]; then
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
fi

# Create index of expected crai files
sed 's;$;.crai;' < ${ALL_CRAMS_INDEX} > ${ALL_CRAIS_INDEX} 

# Create list of ids from the crams files
EXTRACT_ID_PATT='s;.*/\([[:alnum:]_-]*\).*;\1;'
sed "${EXTRACT_ID_PATT}" < ${ALL_CRAMS_INDEX} > ${ALL_CRAMS_IDS}

# Create sorted mapping of id cram crai and restrict by samples list.
paste -d "\t" ${ALL_CRAMS_IDS} ${ALL_CRAMS_INDEX} ${ALL_CRAIS_INDEX} |\
  sort -u -t $'\t' -k 1,1 |\
 join -t $'\t' -j 1 ${SAMPLE_FILE} - > ${TRIM_SAMPLES_MAP}

# Create list of unique random ids
N_IDS=$(wc -l ${TRIM_SAMPLES_MAP} | cut -d' ' -f 1)

let "N_RND = N_IDS + 1000"
let "N_RAW_CHARS = N_RND * RND_ID_CHARS"

cat /dev/urandom |\
  LC_ALL=C tr -dc A-Z0-9 |\
  head -c ${N_RAW_CHARS} |\
  fold -w ${RND_ID_CHARS} |\
  awk '!x[$0]++' |\
  sort |\
  uniq |\
  shuf -n ${N_IDS} > ${RND_CRAMS_IDS}

# Create samples mapping file expected for sequences workflow
paste -d "\t" ${RND_CRAMS_IDS} ${TRIM_SAMPLES_MAP} > ${TRIM_CRAMS_REMAP}

# Create symlinks for both the cram and crai
rm -f ${SEQUENCES_CRAM_DIR}/*
OLD_IFS=$IFS
IFS=$'\t'
while read RND_ID SAMPLE_ID CRAM_PATH CRAI_PATH; do
	ln -f -s "${CRAM_PATH}" "${SEQUENCES_CRAM_DIR}/${SAMPLE_ID}.cram"
	ln -f -s "${CRAI_PATH}.crai" "${SEQUENCES_CRAM_DIR}/${SAMPLE_ID}.cram.crai"
done < ${TRIM_CRAMS_REMAP}
