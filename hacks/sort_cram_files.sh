#!/usr/bin/env bash

# Hack to solve the issue of having all cram files in a single directory.
#  API expects cram files will be sample_id.cram and reside in a two hex char subdirectory.

# Run in the directory containing crams.
# Assume .cram.crai sits in same directory as corresponding .cram

# Create all 256 two character hex subdirectories 00 to ff
for i in {0..255}; do printf '%.2x\n' $i; done | xargs mkdir -p

# Move each file to the subdirectory it's expected in
for f in *.cram; do
  TARGET_DIR=$(sed 's/.cram//' < "$f" | md5sum | cut -c1-2)
  mv "${f}" "${TARGET_DIR}/${f}" 
  mv "${f}.crai" "${TARGET_DIR}/${f}.crai" 
done

echo "Done."

