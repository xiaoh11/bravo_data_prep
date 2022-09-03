import java.nio.file.Paths

grp_size = 5

process select_files {
  input:
  path data_dir from Channel.fromPath(params.cram_base_dir)

  output:
  stdout paths_subset
  path "selected_crams.txt"

  publishDir "result/subset", pattern: "*.txt"

  script:
  """
  ls ${params.cram_files} |\
     shuf -n 10 |\
     tee selected_crams.txt
  """

}

process pileup {
  label "highcpu"

  input:
  // Use chromosome as processing grouping.
  // Output will be also grouped as a result.
  each chromosome from Channel.from(params.chromosomes)
  
  // Split filenames input by line and remove trailing newlines. 
  tuple val(name), file(cram), file(crai) from paths_subset
    .splitText(){ it.trim() }
    .map{ [Paths.get(it).getSimpleName(), file(it), file("${it}.crai")] } 

  // Use a value channel allowing indefinite reuse
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

process aggregate_depths_rnd_1 {
  label "highmem"

  input:
  tuple val(chromosome), 
        file(depth_files), 
        file(depth_tbis) from pileups.groupTuple(size: grp_size)

  output:
  tuple val(chromosome), 
        file("${chromosome}_rnd_1_${task.index}.tsv.gz"),
        file("${chromosome}_rnd_1_${task.index}.tsv.gz.tbi") into agg_rnd_1

  script:

  min_pos = 0
  max_pos = params.chrom_lengths[chromosome]
  step = 1000000
  result_file="${chromosome}_rnd_1_${task.index}.tsv.gz"

  template 'mlr_agg.sh'
}

process aggregate_depths_rnd_2 {
  label "highmem"
  input:
  tuple val(chromosome), file(dat_file_list), file(idx_file_list) from agg_rnd_1
    .toSortedList({ lhs, rhs -> lhs[0] <=> rhs[0] })
    .flatMap()
    .groupTuple(size: grp_size, remainder: true)

  output:
  tuple val(chromosome), 
        file("${chromosome}_rnd_2_*.tsv.gz"),
        file("${chromosome}_rnd_2_*.tsv.gz.tbi") into agg_rnd_2
  
  script:

  min_pos = 0
  max_pos = params.chrom_lengths[chromosome]
  step = 1000000
  result_file="${chromosome}_rnd_2_${task.index}.tsv.gz"

  template 'mlr_agg.sh'
}

process aggregate_depths_rnd_3 {
  label "highmem"

  input:
  tuple val(chromosome), file(dat_file_list), file(idx_file_list) from agg_rnd_2
    .toSortedList({ lhs, rhs -> lhs[0] <=> rhs[0] })
    .flatMap()
    .groupTuple(size: grp_size, remainder: true)

  output:
  tuple val(chromosome), 
        file("${chromosome}_rnd_3_*.tsv.gz"),
        file("${chromosome}_rnd_3_*.tsv.gz.tbi") into agg_rnd_3

  publishDir "result/agg"
  
  script:

  min_pos = 0
  max_pos = params.chrom_lengths[chromosome]
  step = 1000000
  result_file="${chromosome}_rnd_3_${task.index}.tsv.gz"

  template 'mlr_agg.sh'
}
