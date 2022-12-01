process pileup {
  memory "0.75 GB"

  input:
  // get ID, ID.cram, and ID.cram.crai from glob of crams
  tuple val(name), file(cram), file(crai) from Channel 
                                                 .fromPath(params.cram_files)
                                                 .map{ ["${it.getSimpleName()}",
                                                        file("${it}"),
                                                        file("${it}.crai")] }
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
  samtools view -q 20 -F 0x0704 -h ${cram} ${chromosome} |\
    samtools calmd -AEr - ${ref} |\
    bam clipOverlap --in - --out - |\
    samtools mpileup -f ${ref} -Q 20 - |\
    cut -f1,2,4 | bgzip > ${chromosome}.${name}.depth.gz
  tabix -s1 -b2 -e2 ${chromosome}.${name}.depth.gz
  """
}

process aggregate {
   memory "6 GB"

   input:
   tuple val(chromosome), file(depth_files), file(depth_tbis) from pileups.groupTuple()

   output:
   tuple val(chromosome), file("${chromosome}.full.json.gz"), file("${chromosome}.full.json.gz.tbi") into aggregated

   publishDir "result/full", pattern: "*.full.json.gz*"

   """
   find . -name "${chromosome}.*.depth.gz" > files_list.txt
   aggregate.py -i files_list.txt -o ${chromosome}.full.json.gz
   """
}

process prune {
   memory "4 GB"

   input:
   set val(chromosome), file(full_json), file(full_json_tbi) from aggregated
   each limit from Channel.from(params.prune_limits)

   output:
   tuple file("${chromosome}.bin_${limit}.json.gz"), file("${chromosome}.bin_${limit}.tsv.gz.tbi")

   publishDir "result/bin_${limit}", pattern: "*.bin_*.json.gz*"

   """
   prune.py -i ${full_json} -l ${limit} -o ${chromosome}.bin_${limit}.tsv.gz
   """
}
