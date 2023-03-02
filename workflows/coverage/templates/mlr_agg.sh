# Enumerate inputs
mapfile -d '' INFILES < <(find . -name '*.gz' -print0)

# Create file to append aggregate results
cat /dev/null > agg.gz 

ITT=0
POS=${min_pos}
while [ \${POS} -lt ${max_pos} ]
do
  echo "\${POS}"
  # Calculate end of block
  let "END = POS + ${step}"


  # Create pipes for tabix reading
  PIPES=()
  PIPE_COUNT=0
  for INFILE in \${INFILES[@]}
  do
    PIPE_NAME=pipe_\${PIPE_COUNT}
    mkfifo \${PIPE_NAME}
    PIPES+=( \${PIPE_NAME} )

    # Start reading depth file into pipe
    (tabix \${INFILE} ${chromosome}:\${POS}-\${END} > \${PIPE_NAME}) &

    let "PIPE_COUNT = PIPE_COUNT + 1"
  done

  # Aggregate depths from depth file chunks
  mlr -N --tsv 'nest' --ivar ";" -f 3 \${PIPES[@]} |\
    sort --numeric-sort --key=2 |\
    bgzip >> agg.gz

  # Clean up pipes
  for P in \${PIPES[@]}
  do
    rm \${P}
  done

  # Advance the position beyond current end
  let "POS = END + 1"
  let "ITT = ITT + 1"
done

# Work around tabix not working for some concatenated bgzip files.
zcat agg.gz | bgzip --threads 2 > ${result_file}

# Index the results
tabix -s 1 -b 2 -e 2 ${result_file}
