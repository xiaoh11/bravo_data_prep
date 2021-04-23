process random {

  input:
  file vcf from Channel.fromPath(params.vcfs_path)

  // Use a value channel for the samples file so it can be read unlimited times.
  file samples_file from Channel.value( file(params.samples_path) )

  output:
  tuple stdout, file("${vcf.baseName}.random.vcf.gz"), file("${vcf.baseName}.random.vcf.gz.tbi") into random_vcf

  publishDir "result/random", pattern: "*.random.vcf.gz"

  script:
  // Handle optional file parameter
  def samples_opt = samples_file.name == 'NO_FILE' ?  "" : "-s $samples_file"
  """
  ${params.random.exec} -i ${vcf} ${samples_opt} -k ${params.random.max_hethom} -e ${params.random.seed} -o ${vcf.baseName}.random.vcf.gz > ${vcf.baseName}.random.log
  tabix ${vcf.baseName}.random.vcf.gz > tabix.log
  tabix -l ${vcf.baseName}.random.vcf.gz > chromosomes.txt
  if [[ \$(wc -l <chromosomes.txt) -gt 1 ]]; then
     echo "Multiple chromosomes within one input file are not allowed." 1>&2
     exit 1
  elif [[ \$(wc -l <chromosomes.txt) -eq 0  ]]; then
     printf "empty"
  fi
  printf "`head -n1 chromosomes.txt`"
  """
}


process sequences {

  input:
  // Include index file as input so it's symlinked in the work dir
  tuple val(chromosome), file(vcf), file(vcf_idx) from random_vcf

  // Use value channel to allow reference to be used repeatedly
  file reference from Channel.value( file(params.reference_path) )
  file py_script from Channel.value( file(params.sequences.script_path) )

  // Make text file containing: id /path/to/id.cram
  // Use Reduce to return a value Channel.
  file crams_list from Channel.fromPath(params.cram_files)
                                  .map{ "${it.getSimpleName()} ${it}" }
                                  .collectFile(name: "cram_files_list.txt", newLine: true)
                                  .reduce { result -> result }

  when:
  chromosome != 'empty'

  output:
  tuple val(chromosome), file("${vcf.baseName}.sorted.cram"), file("${vcf.baseName}.sorted.cram.crai") into sequences
  
  publishDir "result/crams/${chromosome}", pattern: "*.sorted.cram*"

  script:
  """
  echo "vcf: ${vcf}"
  echo "idx: ${vcf_idx}"
  python ${py_script} cram -i ${vcf} -r ${reference} -c ${crams_list} -w ${params.sequences.window} -o ${vcf.baseName}.cram > ${vcf.baseName}.log
  ${params.samtools.exec} sort -m ${params.samtools.max_mem} -o ${vcf.baseName}.sorted.cram ${vcf.baseName}.cram
  ${params.samtools.exec} index ${vcf.baseName}.sorted.cram
  """
}
