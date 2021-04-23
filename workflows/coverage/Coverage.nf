process pileup {

  input:
  // get ID, ID.cram, and ID.cram.crai from glob of crams
  tuple val(name), file(cram), file(crai) from Channel 
                                                 .fromPath(params.cram_files)
                                                 .map{ ["${it.getSimpleName()}",
                                                        "${it}",
                                                        "${it}.crai"] }
  each chromosome from Channel.from(params.chromosomes)

  // Use a value channel allowing reuse
  file ref from Channel.value( file(params.reference_path) )
  file ref_idx from Channel.value( file(params.reference_path + ".fai") )

  output:
  tuple val(chromosome), 
        file("${chromosome}.${name}.depth.gz"), 
        file("${chromosome}.${name}.depth.gz.tbi") into pileups

  publishDir "result/pileup/${chromosome}", pattern: "*.depth.gz*"

  script:
  """
  ${params.samtools} view -q 20 -F 0x0704 -h ${cram} ${chromosome} |\
    ${params.samtools} calmd -AEr - ${ref} |\
    bam clipOverlap --in - --out - |\
    ${params.samtools} mpileup -f ${ref} -Q 20 - |\
    cut -f1,2,4 | bgzip > ${chromosome}.${name}.depth.gz
  tabix -s1 -b2 -e3 ${chromosome}.${name}.depth.gz
  """
}


// Were these processes being worked on or vestigial? 
/*process aggregate {
   errorStrategy "ignore"

   input:
   set val(chromosome), file(depth_files) from pileups.groupTuple()

   output:
   set val(chromosome), file("${chromosome}.full.json.gz"), file("${chromosome}.full.json.gz.tbi") into aggregated

   publishDir "results/full", pattern: "*.full.json.gz*"

   """
   find . -name "${chromosome}.*.depth.gz" > files_list.txt
   ${params.aggregate} -i files_list.txt -o ${chromosome}.full.json.gz
   """
}


process prune {
   errorStrategy "ignore"

   input:
   set val(chromosome), file(full_json), file(full_json_tbi) from aggregated
   each limit from Channel.from(params.prune_limits)

   output:
   set file("${chromosome}.bin_${limit}.json.gz"), file("${chromosome}.bin_${limit}.json.gz.tbi")

   publishDir "results/bin_${limit}", pattern: "*.bin_*.json.gz*"

   """
   ${params.prune} -i ${full_json} -l ${limit} -o ${chromosome}.bin_${limit}.json.gz
   """
}*/
