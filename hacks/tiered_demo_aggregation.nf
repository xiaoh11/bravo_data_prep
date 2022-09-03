/* Illustrate sorting and re-grouping depth files for subsequent rounds of aggregation.
  In practice, only the 1st aggregation will get the files already grouped by chromosome.
  Each aggregation process emits: [chrom, date_file, index_file]
  The ordering isn't garunteed so they will have to be collected, sorted, and grouped.
*/

agg_grp_size = 5

dummy_in = Channel.from( "chr11", "chr12", "chr11", "chr12", "chr11", "chr12", "chr11", "chr12", "chr11", "chr12")

// First aggregation step generates initial data and index files
//   In practice the data will come into the first aggregation step already grouped by chromosome.
process scratch_rnd_1 {
  input:
  val chromosome from dummy_in

  output:
  tuple val(chromosome), 
        file("*.dat"),
        file("*.dat.idx") into agg_rnd_1

  script:
  out_dat = "${chromosome}.rnd_1.${task.index}.dat"
  out_idx = "${out_dat}.idx"
  """
  echo "foo" > ${out_dat}
  echo "bar" > ${out_idx}
  """
}

process scratch_rnd_2 {
  input:
  tuple val(chromosome), val(dat_file_list), val(idx_file_list) from agg_rnd_1
    .toSortedList({ lhs, rhs -> lhs[0] <=> rhs[0] })
    .flatMap()
    .groupTuple(size: agg_grp_size, remainder: true)

  output:
  tuple val(chromosome), 
        file("*.dat"),
        file("*dat.idx") into agg_rnd_2
  
  publishDir "result/scratch"

  script:
  out_dat = "${chromosome}.rnd_2.${task.index}.dat"
  out_idx = "${out_dat}.idx"
  """
  find . -name '*.dat' > ${out_dat}
  echo "duq" > ${out_idx}
  """
}

process scratch_rnd_3 {
  input:
  tuple val(chromosome), val(dat_file_list), val(idx_file_list) from agg_rnd_2
    .toSortedList({ lhs, rhs -> lhs[0] <=> rhs[0] })
    .flatMap()
    .groupTuple(size: 3, remainder: true)

  output:
  tuple val(chromosome), 
        file("*.dat"),
        file("*.dat.idx") into agg_rnd_3
  
  publishDir "result/scratch"

  script:
  out_dat = "${chromosome}.rnd_3.${task.index}.dat"
  out_idx = "${out_dat}.idx"
  """
  find . -name '*.dat' > ${out_dat}
  echo "duq" > ${out_idx}
  """
}

agg_rnd_3.subscribe onNext: { println it }, onComplete: { println 'Done' }
