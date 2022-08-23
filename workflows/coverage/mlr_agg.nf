process pileup {
  // Debug: Do not run
  when false

  label "highcpu"

  input:
  // Get ID, ID.cram, and ID.cram.crai from glob of crams
  // Random sample 10 of the crams for processing
  tuple val(name), file(cram), file(crai) from Channel 
                                                 .fromPath(params.cram_files)
                                                 .map{ ["${it.getSimpleName()}",
                                                        file("${it}"),
                                                        file("${it}.crai")] }
                                                 .randomSample(10)
  each chromosome from Channel.from(params.chromosomes)

  // Use a value channel allowing reuse
  file ref from Channel.value( file(params.reference_path) )
  file ref_idx from Channel.value( file(params.reference_path + ".fai") )

  output:
  tuple val(chromosome), 
        file("${chromosome}.${name}.depth.gz"), 
        file("${chromosome}.${name}.depth.gz.tbi") into pileups

  publishDir "result/test_pileup/${chromosome}", pattern: "*.depth.gz*"

  script:
  """
  samtools view -q 20 -F 0x0704 -h ${cram} ${chromosome} |\
    samtools calmd -AEr - ${ref} |\
    bam clipOverlap --in - --out - |\
    samtools mpileup -f ${ref} -Q 20 - |\
    cut -f1,2,4 | bgzip > ${chromosome}.${name}.depth.gz
  tabix -s1 -b2 -e2 ${chromosome}.${name}.depth.gz
  """
}

// Mock up using cached output from previous step.
process aggregate {
   label "highmem"

   // input:
   // tuple val(chromosome), file(depth_files), file(depth_tbis) from pileups.groupTuple()

   tuple file(depth_files), file(depth_tbis) 
      fromFilePairs("/home/grosscol_umich_edu/data_prep/workflows/coverage/result/test_pileup/chr11/*.{depth.gz,depth.gz.tbi}", flat: true)

   val chromsome from "chr11"

   output:
   tuple val(chromosome), file("${chromosome}.full.json.gz"), file("${chromosome}.full.json.gz.tbi") into aggregated

   publishDir "result/test_full", pattern: "*.full.json.gz*"

   """
   MD5=\$(ls *.gz | md5sum | cut -f 1 -d ' ')
   find . -name "${chromosome}.*.depth.gz" > files_list.txt
   aggregate.py -i files_list.txt -o ${chromosome}.full.json.gz
   """
}
