#!/usr/bin/env sh

# Create a bespoke rewriting script for links in a directory (LN_SRC) that currently target
#   files in another diretory (TARGET_DIR). Rewrite links according to sed expression (SUBST)
#
# Use
#  ./make_link_rewrite_script.sh > my_rewrite.sh
#  ./my_rewrite.sh
# Derived from this S.E. [answer](https://superuser.com/a/157832)

LN_SRC="/current/path/to/workflows/coverage/result/vep_ok"
TARGET_DIR="/old/path/to/workflows"
SUBST="s/old/current/"

find ${LN_SRC} -type l \
  -lname "${TARGET_DIR}/*" -printf \
  'ln -nsf "$(readlink "%p"|sed ${SUBST})" "$(echo "%p"|sed ${SUBST})"\n'
