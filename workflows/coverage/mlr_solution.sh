#!/usr/bin/env bash

INFILES=(chr11/chr11.HG00337.depth.gz
         chr11/chr11.HG01069.depth.gz
         chr11/chr11.HG02008.depth.gz
         chr11/chr11.HG02729.depth.gz
         chr11/chr11.HG03113.depth.gz
         chr11/chr11.HG00956.depth.gz
         chr11/chr11.HG01178.depth.gz
         chr11/chr11.HG02390.depth.gz
         chr11/chr11.HG03065.depth.gz
         chr11/chr11.HG03249.depth.gz
)

MAX_POS=$1
MIN_POS=0 

DELTA=$(( ${MAX_POS} - ${MIN_POS} ))
STEP=1000000

if [ -z ${MAX_POS} ]
then
  echo "Must provide maxium position to run to." 1>&2
  exit 1
fi

echo "From ${MIN_POS} to ${MAX_POS}"
echo "${DELTA} positions ${STEP} at a time"

# Initialize result file
RESULT_FILE=result.tsv.gz
cat /dev/null > ${RESULT_FILE}

ITT=0
POS=${MIN_POS}
while [ ${POS} -lt ${MAX_POS} ]
do
  # Calculate end of block
  END=$((${POS} + ${STEP}))

  # Create pipes for tabix reading
  PIPES=()
  PIPE_COUNT=0
  for INFILE in ${INFILES[@]}
  do
    PIPE_NAME=pipe_${PIPE_COUNT}
    mkfifo $PIPE_NAME
    PIPES+=( $PIPE_NAME )

    # Start reading depth file into pipe
    (tabix ${INFILE} chr11:${POS}-${END} > ${PIPE_NAME}) &

    let "PIPE_COUNT++"
  done

  # Aggregate depths from depth file chunks
  mlr -N --tsv 'nest' --ivar ";" -f 3 ${PIPES[@]} |\
    sort --numeric-sort --key=2 |\
    bgzip >> ${RESULT_FILE}

  # Clean up pipes
  for P in ${PIPES[@]}
  do
    rm ${P}
  done

  POS=$((${END} + 1))
  let "ITT++"
done

echo "Done in ${ITT} iterations."

