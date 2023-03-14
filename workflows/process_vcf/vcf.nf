/*
  Generate samples file when NO_FILE is specified.
*/
process generate_samples_file {
  input:
  // Assume all samples are present in the first VCF
  tuple val(id), file(bcf), file(idx) from Channel.fromFilePairs(params.bcfs_full, flat: true).first()

  output:
  file("samples.txt") into samples_file_chan

  when:
  params.samples_path == 'NO_FILE'

  script:
  """
  bcftools query -l $bcf > samples.txt
  """
}


/*
  Use samples file if provided.
*/
if(params.samples_path != 'NO_FILE') {
  samples_file_chan = Channel.fromPath(params.samples_path)
}


/*
  Calculate depth histograms from full BCFs.
    Apply filter to input file pairs to exclude M and Y chromosomes.
*/
process calc_histograms {
  label "small"

  input:
  tuple val(id), file(bcf), file(idx) from Channel
    .fromFilePairs(params.bcfs_full, flat: true)

  // Use a value channel for the samples file so it can be read unlimited times.
  file samples_file from samples_file_chan.first()

  output:
  tuple chr_str, file("${bcf.baseName}.hist.bcf"), file("${bcf.baseName}.hist.bcf.csi") into histogramed

  script:
  chr_str = (bcf.baseName =~ /chr[0-9X]{1,2}/)[0]
  outfile = "${bcf.baseName}.hist.bcf"

  """
  ComputeAlleleCountsAndHistograms -i ${bcf} -s ${samples_file} \
    -o ${outfile} > ${bcf.baseName}.hist.log

  tabix --csi -f ${outfile} > tabix.log
  """
}


/*
  Annotate INFO from sites only vcfs into histogram vcfs
*/
process annotate_bcfs {
  label "annotate"

  input:
  tuple val(chr), file(bcf), file(idx) from histogramed
  path(sites_dir) from Channel.fromPath(params.bcfs_sites_dir).first()

  output:
  tuple val(chr), file("${bcf.baseName}.anno.bcf"), file("${bcf.baseName}.anno.bcf.csi") into annotated

  script:
  out_bcf = "${bcf.baseName}.anno.bcf"
  out_idx = "${bcf.baseName}.anno.bcf.csi"
  sites_bcf = "${sites_dir}/freeze10.topmed.${chr}.filtered.sites.bcf"
  """
  bcftools annotate -a ${sites_bcf} \
    --output-type b \
    --columns ${params.anno_fields} \
    --output ${out_bcf} \
    ${bcf}

  # Index resulting file
  tabix --csi -f ${out_bcf}
  """
}


process vep {
  label "vep"
  memory { task.attempt > 1 ? 50.GB : 10.GB }

  publishDir "result/vep/${chromosome}", pattern: "*.vep.gz*"

  input:
  tuple val(chromosome), file(bcf), file(idx) from annotated

  when:
  chromosome != 'empty'

  output:
  tuple val(chromosome), file("${bcf.baseName}.vep.gz"), file("${bcf.baseName}.vep.gz.csi") into vep

  script:
  outfile = "${bcf.baseName}.vep.gz"
  plugin_opts = ["LoF",
                 "loftee_path:${params.vep.loftee_path}",
                 "human_ancestor_fa:${params.vep.loftee_human_ancestor_fa}",
                 "conservation_file:${params.vep.loftee_conservation_file}",
                 "gerp_bigwig:${params.vep.loftee_gerp_bigwig}"
                ].join(",")
  """
  bcftools view ${bcf} |\
  vep -o ${outfile} \
    --fasta ${params.vep.ref_fasta} \
    --dir_cache ${params.vep.cache} \
    --dir_plugins ${params.vep.loftee_path} \
    --plugin ${plugin_opts} \
    ${params.vep.static_flags}

  tabix -f --csi ${outfile}
  """
}


/* Depends on:
  - incoming vcf split by position and ordered.
  - underscore delime position in filename:  chr_start_end
  demo.chr12_132700001_132800000.genotypes.hist.anno.vep.gz
  demo.chr12_132800001_132900000.genotypes.hist.anno.vep.gz
  demo.chr12_132900001_133000000.genotypes.hist.anno.vep.gz
*/
process concat {
  label "small"

  input:
  tuple val(chromosome), file(vcfs), file(indices) from vep.groupTuple()

  output:
  tuple val(chromosome), file("${chromosome}.aggr.vep.concat.gz"), file("${chromosome}.aggr.vep.concat.gz.csi") into concat

  publishDir "result/concat", pattern: "*.aggr.vep.concat.gz*"

  """
  # Sort file names by first position in file name
  mapfile -t ALL_FILE_NAMES < <(find . -name "*.vep.gz" | sort -n -k 2 -t '_')

  cat /dev/null > sorted_files_list.txt
  cat /dev/null > omitted_files_list.txt

  # Omit VCF files with no calls as they bork concat.
  for F_NAME in "\${ALL_FILE_NAMES[@]}"
  do
    N_LINES=\$(bcftools view --no-header "\${F_NAME}" | head | wc -l)

    if [ \${N_LINES} -gt 0 ]
    then
      echo "\${F_NAME}" >> sorted_files_list.txt
    else
      echo "\${F_NAME}" >> omitted_files_list.txt
    fi
  done

  bcftools concat -f sorted_files_list.txt --naive -O z -o ${chromosome}.aggr.vep.concat.gz
  tabix -f --csi ${chromosome}.aggr.vep.concat.gz
  """
}


process cadd {
  label "small"

  input:
  tuple val(chromosome), file(vcf), file(index) from concat
  file cadds from Channel.fromPath(params.cadd.tsv_path).collect()
  file cadds_indices from Channel.fromPath(params.cadd.tsv_path + ".tbi").collect()

  output:
  tuple val(chromosome),
        file("${vcf.baseName}.cadds.gz"),
        file("${vcf.baseName}.cadds.gz.csi") into(cadds1, cadds2)

  publishDir "result/cadds", pattern: "*.cadds.gz*"

  script:
  def script_file = file(params.cadd.script)
  """
  python ${script_file} -i ${vcf} -c ${cadds} -o ${vcf.baseName}.cadds.gz
  tabix -f --csi ${vcf.baseName}.cadds.gz
  """
}

process percentiles {
  label "highmem"

  input:
  val qc_metric from Channel.from(params.percentiles.qc_metrics)
  file vcfs from cadds1.map { it[1] }.collect()

  output:
  file "${qc_metric}.variant_percentile.vcf.gz" into variant_percentiles
  file "${qc_metric}.variant_percentile.vcf.gz.csi" into variant_percentiles_indices
  file "${qc_metric}.all_percentiles.json.gz" into metric_summaries

  publishDir "result/percentiles", pattern: "*.variant_percentile.vcf.gz*"
  publishDir "result/percentiles", pattern: "*.all_percentiles.json.gz"

  """
  if [[ "${qc_metric}" = "QUAL" ]]; then
     ComputePercentiles -t 4 -i ${vcfs} -p 10 -o ${qc_metric}
  else
     message=`bcftools view -h ${vcfs[0]} | grep 'INFO=<ID=${qc_metric}' | grep -oP 'Description="\\K.+(?=")'`
     ComputePercentiles -t 4 -i ${vcfs} -p 10 -m ${qc_metric} -d \"\${message}\" -o ${qc_metric}
  fi
  tabix -f --csi ${qc_metric}.variant_percentile.vcf.gz
  """
}


process merge_vcf {
  label "merge"

  input:
  tuple val(chromosome), file(vcf), file(vcf_index) from cadds2
  val qc_metrics from Channel.from(params.percentiles.qc_metrics).toList().map { it.join(" ") }
  file variant_percentiles from variant_percentiles.collect()
  file variant_percentiles_indices from variant_percentiles_indices.collect()

  output:
  tuple file("${chromosome}.bravo.vcf.gz"), file("${chromosome}.bravo.vcf.gz.csi") into merged_vcf

  publishDir "result/final/vcfs", pattern: "*.bravo.vcf.gz*"

  """
  cp ${vcf} temp.vcf.gz
  cp ${vcf_index} temp.vcf.gz.csi
  for qc_metric in ${qc_metrics}; do
     bcftools annotate -a \${qc_metric}.variant_percentile.vcf.gz \
       -c INFO/\${qc_metric}_PCTL \
       -Oz -o ${chromosome}.bravo.vcf.gz \
       temp.vcf.gz
     cp ${chromosome}.bravo.vcf.gz temp.vcf.gz
     tabix -f --csi temp.vcf.gz
  done
  tabix -f --csi ${chromosome}.bravo.vcf.gz
  rm temp.vcf.gz*
  """
}


process merge_metrics {
  label "small"

  input:
  file metrics from metric_summaries.collect()

  output:
  file "metrics.json.gz" into merged_metrics

  publishDir "result/final/qc_metrics", pattern: "metrics.json.gz"

  """
  for f in ${metrics}; do gzip -dc "\${f}"; done | gzip > metrics.json.gz
  """
}
