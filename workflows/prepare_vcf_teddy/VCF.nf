/*
  Generate samples file when not provided
*/
if(params.samples_path != 'NO_FILE') {
  samples_file_chan = Channel.fromPath(params.samples_file)
} 

process generate_samples_file {
  input:
  // Assume all samples are present in the first VCF
  tuple val(id), file(vcf), file(idx) from Channel.fromFilePairs(params.vcfs_glob, flat: true).first()

  output:
  file("samples.txt") into samples_file_chan

  when: 
  params.samples_path == 'NO_FILE'

  script:
  """
  bcftools query -l $vcf > samples.txt
  """
}

process aggregate {
  input:
  tuple val(id), file(vcf), file(idx) from Channel.fromFilePairs(params.vcfs_glob, flat: true)

  // Use a value channel for the samples file so it can be read unlimited times.
  file samples_file from samples_file_chan

  output:
  tuple stdout, file("${vcf.baseName}.aggr.gz") into aggregated

  publishDir "results/aggr", pattern: "*.aggr.gz"

  script:
  """
  ${params.aggregate} -i ${vcf} -s ${samples} -f ${params.qc_metrics} -o ${vcf.baseName}.aggr.gz > ${vcf.baseName}.log
  tabix ${vcf.baseName}.aggr.gz > tabix.log
  tabix -l ${vcf.baseName}.aggr.gz > chromosomes.txt
  if [[ \$(wc -l <chromosomes.txt) -gt 1 ]]; then
     echo "Multiple chromosomes within one input file are not allowed." 1>&2
     exit 1
  elif [[ \$(wc -l <chromosomes.txt) -eq 0  ]]; then
     printf "empty"
  fi
  printf "`head -n1 chromosomes.txt`"
  """
}


process vep {
   input:
   tuple val(chromosome), file(vcf) from aggregated

   when:
   chromosome != 'empty'

   output:
   tuple val(chromosome), file("${vcf.baseName}.vep.gz"), file("${vcf.baseName}.vep.gz.tbi") into vep
   
   publishDir "results/vep/${chromosome}", pattern: "*.vep.gz*"

   """
   ${params.vep.vep} -i ${vcf} -o ${vcf.baseName}.vep.gz --format vcf --cache --offline --vcf --compress_output bgzip --no_stats --dir_plugins ${params.vep.loftee_path} --plugin LoF,loftee_path:${params.vep.loftee_path},human_ancestor_fa:${params.vep.loftee_human_ancestor_fa},conservation_file:${params.vep.loftee_conservation_file},gerp_bigwig:${params.vep.loftee_gerp_bigwig} ${params.vep.flags}
   tabix ${vcf.baseName}.vep.gz
   """
}


process concat {
   input:
   tuple val(chromosome), file(vcfs), file(indices) from vep.groupTuple()

   output:
   tuple val(chromosome), file("${chromosome}.aggr.vep.concat.gz"), file("${chromosome}.aggr.vep.concat.gz.tbi") into concat

   publishDir "results/concat", pattern: "*.aggr.vep.concat.gz*"
 
   """
   find . -name "*.vep.gz" > files_list.txt
   bcftools concat -a -f files_list.txt -Oz -o ${chromosome}.aggr.vep.concat.gz
   tabix ${chromosome}.aggr.vep.concat.gz
   """
}


process cadd {
  input:
  tuple val(chromosome), file(vcf), file(index) from concat
  file cadds from Channel.fromPath(params.cadd.tsv_path).collect()
  file cadds_indices from Channel.fromPath(params.cadd.tsv_path + ".tbi").collect()

  output:
  tuple val(chromosome), file("${vcf.baseName}.cadds.gz"), file("${vcf.baseName}.cadds.gz.tbi") into cadds

  publishDir "results/cadds", pattern: "*.cadds.gz*"

  script:
  def script_file = file(params.cadd.script)
  """
  python ${script_file} -i ${vcf} -c ${cadds} -o ${vcf.baseName}.cadds.gz
  tabix ${vcf.baseName}.cadds.gz
  """
}


(cadds1, cadds2) = cadds.into(2)


process percentiles {
   input:
   val qc_metric from Channel.from(params.percentiles.qc_metrics)
   file vcfs from cadds1.map { it[1] }.collect()

   output:
   file "${qc_metric}.variant_percentile.vcf.gz" into variant_percentiles
   file "${qc_metric}.variant_percentile.vcf.gz.tbi" into variant_percentiles_indices
   file "${qc_metric}.all_percentiles.json.gz" into metric_summaries

   publishDir "results/percentiles", pattern: "*.variant_percentile.vcf.gz*"
   publishDir "results/percentiles", pattern: "*.all_percentiles.json.gz"

   """
   if [[ "${qc_metric}" = "QUAL" ]]; then
      ${params.percentiles.exec} -t 4 -i ${vcfs} -p 10 -o ${qc_metric}
   else
      message=`bcftools view -h ${vcfs[0]} | grep 'INFO=<ID=${qc_metric}' | grep -oP 'Description="\\K.+(?=")'`
      ${params.percentiles.exec} -t 4 -i ${vcfs} -p 10 -m ${qc_metric} -d \"\${message}\" -o ${qc_metric}
   fi
   tabix ${qc_metric}.variant_percentile.vcf.gz
   """
}


process merge_vcf {
   input:
   tuple val(chromosome), file(vcf), file(vcf_index) from cadds2
   val qc_metrics from Channel.from(params.percentiles.qc_metrics).toList().map { it.join(" ") }
   file variant_percentiles from variant_percentiles.collect()
   file variant_percentiles_indices from variant_percentiles_indices.collect()

   output:
   tuple file("${chromosome}.bravo.vcf.gz"), file("${chromosome}.bravo.vcf.gz") into merged_vcf

   publishDir "results/final/vcfs", pattern: "*.bravo.vcf.gz*"

   """
   cp ${vcf} temp.vcf.gz
   tabix temp.vcf.gz
   for qc_metric in ${qc_metrics}; do
      bcftools annotate -a \${qc_metric}.variant_percentile.vcf.gz -c INFO/\${qc_metric}_PCTL temp.vcf.gz -Oz -o ${chromosome}.bravo.vcf.gz
      cp ${chromosome}.bravo.vcf.gz temp.vcf.gz
      tabix temp.vcf.gz
   done
   tabix ${chromosome}.bravo.vcf.gz
   rm temp.vcf.gz*
   """
}


process merge_metrics {
   input:
   file metrics from metric_summaries.collect()

   output:
   file "metrics.json.gz" into merged_metrics

   publishDir "results/final/qc_metrics", pattern: "metrics.json.gz"

   """
   for f in ${metrics}; do gzip -dc "\${f}"; done | gzip > metrics.json.gz
   """
}
