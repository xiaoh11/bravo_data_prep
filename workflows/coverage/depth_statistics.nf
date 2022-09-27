import java.nio.file.Paths

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
    shuf -n $params.n_indiv |\
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

/****************************************************************************** 
 BEGIN Aggregation Rounds 
  Read pileups $grp_size at a time using Miller (mlr) to compile depths at
  each position to semicolon delimted list of depths.
******************************************************************************/

process aggregate_depths_rnd_1 {
  label "highmem"

  input:
  tuple val(chromosome), 
        file(depth_files), 
        file(depth_tbis) from pileups.groupTuple(size: params.grp_size)

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
    .groupTuple(size: params.grp_size, remainder: true)

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
    .groupTuple(size: params.grp_size, remainder: true)

  output:
  tuple val(chromosome), 
        file("${chromosome}_rnd_3_*.tsv.gz"),
        file("${chromosome}_rnd_3_*.tsv.gz.tbi") into agg_rnd_3

  script:

  min_pos = 0
  max_pos = params.chrom_lengths[chromosome]
  step = 1000000
  result_file="${chromosome}_rnd_3_${task.index}.tsv.gz"

  template 'mlr_agg.sh'
}

process aggregate_depths_rnd_4 {
  label "highmem"

  input:
  tuple val(chromosome), file(dat_file_list), file(idx_file_list) from agg_rnd_3
    .toSortedList({ lhs, rhs -> lhs[0] <=> rhs[0] })
    .flatMap()
    .groupTuple(size: params.grp_size, remainder: true)

  output:
  tuple val(chromosome), 
        file("${chromosome}_rnd_4_*.tsv.gz"),
        file("${chromosome}_rnd_4_*.tsv.gz.tbi") into agg_rnd_4

  publishDir "result/depth_aggregation"
  
  script:

  min_pos = 0
  max_pos = params.chrom_lengths[chromosome]
  step = 1000000
  result_file="${chromosome}_rnd_4_${task.index}.tsv.gz"

  template 'mlr_agg.sh'
}

/*** END Aggregation Rounds ***************************************************/

/*******************************************************************************
 Summarize Depths
  From aggregated depths, caluclate the start and end positions for each input 
  file such that each will be handled in $num_chunks sequential parts.
  Depth staticstics are tabulated on each smaller part.  Finally, part summary
  data is concatenated back into depth statistics file covering one chromosome.
*******************************************************************************/

process prep_summarization {
  label "anyqueue"

  input:
  tuple val(chrom), file(data_path), file(idx_path) from agg_rnd_4  

  each part_num from Channel.from(1..params.num_chunks)

  output:
  tuple val(chrom), env(START), env(END), file(data_path), file(idx_path) into depth_chunking

  script:
  max_pos = params.chrom_lengths[chrom]
  """
  MIN_POS=\$(zcat $data_path | head -n 1 | cut -f2)

  let "INC = 1 + (${max_pos} - MIN_POS) / ${params.num_chunks}"
  let "END = MIN_POS + (INC * $part_num)"
  let "START = 1 + END - INC"
  """
}

process summarize_depths {
  label "anyqueue"

  input:
  tuple val(chrom), val(start_pos), val(end_pos), file(data_file), file(idx_file) from depth_chunking

  output:
  tuple val(chrom), file("${summ_path}") into summary_parts 
  
  script:
  summ_path = "${chrom}.${start_pos.padLeft(9,'0')}.summary.tsv.gz"
  pos_range = "${chrom}:${start_pos}-${end_pos}"

  """
  tabix ${data_file} ${pos_range} |\
   calc_agg_pileups.py -n ${params.n_indiv} |\
   bgzip > ${summ_path} 
  """
}

process combine_summaries {
  label "anyqueue"

  input:
  tuple val(chrom), file(dat_file_list) from summary_parts
    .groupTuple()

  output:
  tuple val(chrom), file(out_path), file(out_tbi) into depth_summaries

  publishDir "result/depth_summary"

  script:
  out_path = "${chrom}_depth_summary.tsv.gz"
  out_tbi = "${chrom}_depth_summary.tsv.gz.tbi"
  """
  find . -name '*.tsv.gz' | sort | xargs -I {} cat {} >> tmp.out
  mv tmp.out ${out_path}
  tabix -s 1 -b 2 -e 3 ${out_path}
  """
}

process prune {
  input:
  set val(chromosome), file(full_depth) from depth_summaries
  each limit from Channel.from(params.prune_limits)

  output:
  tuple file("${chromosome}.bin_${limit}.tsv.gz")

  publishDir "result/bin_${limit}", pattern: "*.bin_*.tsv.gz*"

  """
  prune.py -i ${full_depth} -l ${limit} -o ${chromosome}.bin_${limit}.tsv.gz
  """
}
