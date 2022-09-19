#!/usr/bin/env bash
# Create samples.tsv for sequences workflow with randomly generated new IDs.
# Assume cram symlinks have been made in workflows/sequences/data/crams

CRAMS_DIR="data/crams"
HACKS_DIR="/home/grosscol_umich_edu/data_prep/hacks"
DEST_FILE="${HACKS_DIR}/198_samples.tsv"
# Number of characters per random ID
ID_SIZE=7

# Read original IDs
readarray -t ARR_198 < ${HACKS_DIR}/198_samples.txt
N_IDS=${#ARR_198[@]}

# Generate Random Anonymous IDs and a few extra to handle duplication.
let "N_RND = N_IDS + 100"
let "N_RAW_CHARS = N_RND * ID_SIZE"

echo "Samples list has ${N_IDS} ids."
echo "Generating ${N_RND} random ${ID_SIZE} character IDs."

# Awk one liner to ensure unique ids. https://stackoverflow.com/a/11532197
readarray -t RND_IDS < <( cat /dev/urandom |\
  LC_ALL=C tr -dc A-Z0-9 |\
  head -c ${N_RAW_CHARS} |\
  fold -w ${ID_SIZE} |\
  awk '!x[$0]++')

echo "Generated ${#RND_IDS[@]} unique random ids"

# Truncate destination file
cat /dev/null > ${DEST_FILE}

# Append to destination file
for (( i=0; i < ${N_IDS}; i++ ))
do
  ID=${ARR_198[$i]}
  RND_ID=${RND_IDS[$i]}
  CRAM_PATH="${CRAMS_DIR}/${ID}.final.cram"
  CRAI_PATH="${CRAMS_DIR}/${ID}.final.cram.crai"

  echo -e "${RND_ID}\t${ID}\t${CRAM_PATH}\t${CRAI_PATH}" >> ${DEST_FILE}
done


# for ID in ${ARR_198[@]}
# do
#   # Generate sequence random upper case and numeric characters
#   RND_ID=$(cat /dev/urandom | LC_ALL=C tr -dc A-Z0-9 | head -c7)
# 
# done

