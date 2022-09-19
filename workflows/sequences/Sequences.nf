/* Handle Samples File
  
  Samples file channel: file with New_ID\tOld_ID 
    Should be a value channel to permit indefinite re-use.
    Either generated from list of crams or provided as text file.

  Samples info channel:  new_id, old_id, path/to/cram, path/to/crai
    Either generated from list of crams or parsed from provided samples file.
*/
if(params.samples_path == 'NO_FILE') {
  samples_file_chan = Channel.fromPath(params.cram_files)
                             .map{ "${it.getSimpleName()}\t${it.getSimpleName()}" }
                             .collectFile(name: "generated_samples_list.tsv", newLine: true)
                             .first()

  samples_info_chan = Channel.fromPath(params.cram_files)
                             .map{ 
                               ["${it.getSimpleName()}","${it.getSimpleName()}",file("${it}"),file("${it}.crai")]
                             }
} else {
  samples_file_chan = Channel.fromPath(params.samples_path)
                             .first()

  samples_info_chan = Channel.from(file(params.samples_path).readLines())
                             .map { 
                               line -> def fields = line.split();
                               [fields[0], fields[1], file(fields[2]), file(fields[3])]
                             }
}


process variants_by_sample {
  label "highcpu"

  input:
  file vcf from Channel.fromPath("${params.bcfs_path}")
  //file csi from Channel.fromPath("${params.bcfs_path}.csi")
  file samples_file from samples_file_chan

  output:
  file "*.variant_map.tsv.gz" optional true into variant_maps
  file "*.sample_map.tsv.gz" optional true into sample_maps

  script:
  """
  VCF_HASH=`basename ${vcf} | md5sum | cut -f1 -d" "`
  ${params.random.exec} -i ${vcf} -s ${samples_file} -k ${params.random.max_hethom} -e ${params.random.seed} -o \${VCF_HASH}
  """
}


process merge_sample_maps {
  label "highmem"

  input:
  file sample_map_list from sample_maps.collect()

  output:
  tuple file("merged_sample_map.tsv.gz"), 
        file("merged_sample_map.tsv.gz.gzi"), 
        file("merged_sample_map.tsv.gz.offsets") into merged_sample_map 

  script:
  first_sample_map = sample_map_list.first()
  """
  {
     gzip -dc $first_sample_map | grep "^#";
     for f in *.sample_map.tsv.gz; do
        gzip -dc \${f} | grep -v "^#" || true;
     done | sort -k 1,1 -k2,2 -k3n,3n --compress-program=gzip --temporary-directory=`pwd`;
  } | bgzip -c -i -I merged_sample_map.tsv.gz.gzi > merged_sample_map.tsv.gz;

  # Create tab delim index of: sample ID, offset, length
  bgzip -dc merged_sample_map.tsv.gz |\
    grep -bPo "^[^\t#]*" |\
    awk -F":" '{if (!(\$2 in beg)) beg[\$2] = \$1; end[\$2] = \$1;}\
               END { for (k in beg) {print k, beg[k], end[k] - beg[k]}}' > merged_sample_map.tsv.gz.offsets
  """
}

//  This generates the backing variant map for BRAVO API
process merge_variant_maps {
  label "highmem"

  input:
  file variant_map_list from variant_maps.collect()

  output:
  file "merged_variant_map.tsv.gz*"

  publishDir "result/", pattern: "merged_variant_map.tsv.gz*"

  script:
  first_variant_map = variant_map_list.first()
  """
  {
     gzip -dc $first_variant_map | grep "^#";
     for f in *.variant_map.tsv.gz; do
        gzip -dc \${f} | grep -v "^#" || true;
     done | sort -k 1,1 -k2n,2n --compress-program=gzip --temporary-directory=`pwd`; 
  } | bgzip -c > merged_variant_map.tsv.gz;
  tabix -s1 -b2 -e2 merged_variant_map.tsv.gz;
  """
}

//  This generates the backing crams for BRAVO API
process extract_sequences {
  label "highmem"

  scratch true
  errorStrategy "ignore"

  input:
  tuple val(new_sample_name), val(sample_name), val(cram), val(crai) from samples_info_chan 

  tuple file(sample_map), 
        file(sample_map_gzi), 
        file(sample_map_offsets) from merged_sample_map

  file reference from Channel.value( file(params.reference_path) )

  output:
  file "${new_sample_name}.cram*" optional true into samples_crams

  publishDir "result/sequences/", pattern: "${new_sample_name}.cram*",
                                  saveAs: { filename -> "${new_sample_name.md5().take(2)}/${filename}" }

  script:
  """
  header_bytes=`sort -k2n,2n ${sample_map_offsets} | head -n1 | cut -f2 -d" "`
  header_bytes=\$((header_bytes * 2))
  offset=`grep "^${sample_name}\\s" ${sample_map_offsets} | cut -f2 -d" "`
  bytes=`grep "^${sample_name}\\s" ${sample_map_offsets} | cut -f3 -d" "`
  bytes=\$((bytes * 2))

  {
     bgzip -b 0 -s \${header_bytes} ${sample_map} | grep "^#";
     bgzip -b \${offset} -s \${bytes} ${sample_map} | grep "^${sample_name}\\s";
  } | bgzip -c > variants.tsv.gz

  ${params.sequences.exec} -c ${cram} -i ${crai} -v variants.tsv.gz \
    -r ${reference} -w ${params.sequences.window} -o ${new_sample_name}.no_sort.cram

  ${params.samtools.exec} sort -m ${params.samtools.max_mem} -O cram \
    --reference ${reference} -o ${new_sample_name}.cram ${new_sample_name}.no_sort.cram

  ${params.samtools.exec} index ${new_sample_name}.cram
  """
}

