process select_files {
  input:
  path data_dir from Channel.fromPath("data")

  output:
  stdout paths_subset

  publishDir "result/dummy_select", pattern: "*.txt"

  script:
  // TODO: List file paths, select 100, output appropriate lists or maps to next process
  """
  ls ${params.cram_files} |\
     shuf -n 10 |\
     tee selected_crams.txt
  """

}

/*
process pileup {
  input:
  // get ID, ID.cram, and ID.cram.crai from glob of crams
  tuple val(name), file(cram), file(crai) from Channel 
                                                 .fromPath(params.cram_files)
                                                 .map{ ["${it.getSimpleName()}",
                                                        file("${it}"),
                                                        file("${it}.crai")] }
                                                 .randomSample(10, 1234)
  each chromosome from Channel.from(params.chromosomes)

  // Use a value channel allowing reuse
  file ref from Channel.value( file(params.reference_path) )
  file ref_idx from Channel.value( file(params.reference_path + ".fai") )

  output:
  tuple val(chromosome), 
        file("${chromosome}.${name}.depth.gz"), 
        file("${chromosome}.${name}.depth.gz.tbi") into pileups

  publishDir "result/dummy_pileup/${chromosome}", pattern: "*.depth.gz*"

  script:
  """

  # Dummy sample data
  echo -e "chr11\t10001\t10\nchr11\t10002\t11\nchr11\t10020\t12\nchr11\t10021\t13\n" |\
  bgzip > ${chromosome}.${name}.depth.gz
  tabix -s1 -b2 -e2 ${chromosome}.${name}.depth.gz
  """
}

process aggregate {
  input:
  tuple val(chromosome), 
        file(depth_files), 
        file(depth_tbis) from pileups.groupTuple(size: 5)

  publishDir "result/dummy_full", pattern: "*.txt*"

  """
  echo "foo" > bar.txt
  """
}
*/
