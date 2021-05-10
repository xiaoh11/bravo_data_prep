#!/usr/bin/env nextlow

chromosome               = params.chromosome
vcfs                     = params.vcfs
samples                  = params.samples
window_size              = params.window_size_bp
counts_exec              = params.counts_exec
add_cadd_script          = params.add_cadd_script
vep                      = params.vep
loftee_path              = params.loftee_path
loftee_flags             = params.loftee_flags
loftee_human_ancestor_fa = params.loftee_human_ancestor_fa
loftee_conservation_file = params.loftee_conservation_file
loftee_gerp_bigwig       = params.loftee_gerp_bigwig
cadd_tsvs                = params.cadd_tsvs

vcfs      = Channel.fromPath(vcfs, relative: true)
samples   = Channel.fromPath(samples, relative: true)
windows   = Channel.from((1..300000000).by(window_size)).map { w -> [ w, w + window_size - 1 ] }
cadds     = Channel.fromPath(cadd_tsvs, relative: true).collect()
cadds_tbi = Channel.fromPath(cadd_tsvs + ".tbi", relative: true).collect()

process compute {

  label "small_mem"

  input:
  set val(vcf), val(window_start), val(window_stop) from vcfs.combine(windows)

  // Use a value channel for the samples file so it can be read unlimited times.
  file samples_file from Channel.value( file(params.samples) )


  output:
  set stdout, file("region_summary.vcf.gz"), file("region_summary.vcf.gz.tbi") into regions


  script:
  // Handle optional file parameter
  def samples_opt = samples_file.name == 'NO_FILE' ?  "" : "-s $samples_file"

  """
  ${counts_exec} -i ${vcf} ${samples_opt} -r ${chromosome}:${window_start}-${window_stop} -o region_summary.vcf.gz > log.txt
  tabix region_summary.vcf.gz
  n_variants=`gzip -dc region_summary.vcf.gz | grep -v "#" | wc -l`
  region_start=`gzip -dc region_summary.vcf.gz | grep -v "#" | head -n1 | cut -f2`
  region_stop=`gzip -dc region_summary.vcf.gz | grep -v "#" | tail -n1 | cut -f2`
  printf "\${n_variants},${chromosome},\${region_start},\${region_stop}"
  """
}


process vep {

   label "small_mem"

   input:
   set val(region_chromosome), val(region_start), val(region_stop), file(vcf), file(tbi) from regions.filter{ it[0].split(",")[0].toInteger() > 0 }.map{ r, f, t -> e = r.split(","); return [e[1], e[2], e[3] , f, t];}

   publishDir "result/counts", pattern: "*[0-9].vcf.gz*", mode: "move"
   publishDir "result/vep_ok", pattern: "*[0-9].vep_ok.vcf.gz*"

   output:
   file "${region_chromosome}_${region_start}_${region_stop}.vcf.gz"
   file "${region_chromosome}_${region_start}_${region_stop}.vcf.gz.tbi"
   file "${region_chromosome}_${region_start}_${region_stop}.vep_ok.vcf.gz" into regions_vep
   file "${region_chromosome}_${region_start}_${region_stop}.vep_ok.vcf.gz.tbi"

   """
   cp ${vcf} ${region_chromosome}_${region_start}_${region_stop}.vcf.gz
   cp ${tbi} ${region_chromosome}_${region_start}_${region_stop}.vcf.gz.tbi
   ${vep} -i ${vcf} -o ${region_chromosome}_${region_start}_${region_stop}.vep_ok.vcf.gz --format vcf --cache --offline --vcf --compress_output bgzip --no_stats --dir_plugins ${loftee_path} --plugin LoF,loftee_path:${loftee_path},human_ancestor_fa:${loftee_human_ancestor_fa},conservation_file:${loftee_conservation_file},gerp_bigwig:${loftee_gerp_bigwig} ${loftee_flags}
   tabix ${region_chromosome}_${region_start}_${region_stop}.vep_ok.vcf.gz
   """
}


process cadd {

  label "small_mem"

  input:
  file vcf from regions_vep
  file cadd from cadds
  file cadd_tbi from cadds_tbi

  publishDir "result/cadd_ok", pattern: "*.cadd_ok.vcf.gz*"

  output:
  file "${vcf}.cadd_ok.vcf.gz"
  file "${vcf}.cadd_ok.vcf.gz.tbi"  
 
  script:
  def script_file = file(params.add_cadd_script)
  """
  python ${script_file} -i ${vcf} -c ${cadd} -o ${vcf}.cadd_ok.vcf.gz
  tabix ${vcf}.cadd_ok.vcf.gz
  """
}

