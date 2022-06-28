#!/usr/bin/env bash

# Find all links under LN_SRC that are currently targeting TARGET_DIR, and rewrite using sed. 
#
# Create a bespoke rewriting script for links in a directory (LN_SRC) that currently target
#   files in another diretory (TARGET_DIR). Rewrite links according to sed expression (SUBST)
#
# Use:
#  ./make_link_rewrite_script.sh > my_relink.sh
#  bash my_relink.sh
# Derived from this S.E. [answer](https://superuser.com/a/157832)

LN_SRC="/current/path/to/workflows/coverage/result/vep_ok"
TARGET_DIR="/old/path/to/workflows"
SUBST="s/sequences_freeze8/sequences/"

# Make array of all links in LN_SRC that point anywhere into problematic TARGET_DIR
OLDIFS=${IFS}
IFS=$'\n'
LINK_ARR=($(find ${LN_SRC} -type l -lname "${TARGET_DIR}/*"))
IFS=${OLDIFS}

# Generate ln commands to overwrite links
for LINK_NAME in "${LINK_ARR[@]}"
do
  NEW_LINK_TARGET=$(readlink "${LINK_NAME}" | sed "${SUBST}" )
  echo "ln -nsf ${NEW_LINK_TARGET} ${LINK_NAME}"
done
