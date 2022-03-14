#!/usr/bin/env bash

# Start by assuming it was the path invoked.
THIS_SCRIPT="$0"

# Handle resolving this script location regardless of symlinking.
# Using ls instead of readlink, because bsd and gnu flavors have different behavior.
while [ -h "$THIS_SCRIPT" ] ; do
  ls=`ls -ld "$THIS_SCRIPT"`
  # Drop everything prior to ->
  link=`expr "$ls" : '.*-> \(.*\)$'`
  if expr "$link" : '/.*' > /dev/null; then
    THIS_SCRIPT="$link"
  else
    THIS_SCRIPT=`dirname "$THIS_SCRIPT"`/"$link"
  fi
done

# Get absolute path to the scripts directory.
SCRIPT_DIR=$(dirname "${THIS_SCRIPT}" | pwd -P)

echo "Making workflow bin directories"
mkdir -p ${SCRIPT_DIR}/workflows/{coverage,prepare_vcf,sequences}/bin

echo "Making bin links"
CPP_BIN=${SCRIPT_DIR}/tools/cpp_tools/cget/bin
ln -sf ${CPP_BIN}/* ${SCRIPT_DIR}/workflows/coverage/bin/
ln -sf ${CPP_BIN}/* ${SCRIPT_DIR}/workflows/prepare_vcf/bin/
ln -sf ${CPP_BIN}/* ${SCRIPT_DIR}/workflows/sequences/bin/

PY_BIN=${SCRIPT_DIR}/tools/py_tools
ln -sf ${PY_BIN}/* ${SCRIPT_DIR}/workflows/coverage/bin/
ln -sf ${PY_BIN}/* ${SCRIPT_DIR}/workflows/prepare_vcf/bin/
ln -sf ${PY_BIN}/* ${SCRIPT_DIR}/workflows/sequences/bin/

if [ -f /opt/ensembl-vep/vep ]; then
	ln -sf /opt/ensembl-vep/vep ${SCRIPT_DIR}/workflows/coverage/bin/vep
	ln -sf /opt/ensembl-vep/vep ${SCRIPT_DIR}/workflows/prepare_vcf/bin/vep
	ln -sf /opt/ensembl-vep/vep ${SCRIPT_DIR}/workflows/sequences/bin/vep
else
  echo "vep link NOT created"
  echo "vep not found at /opt/ensembl-vep/vep"
fi
