crams = Channel.from(file(params.crams_list_path).readLines()).map { line -> def fields = line.split();  [fields[0], fields[1], fields[2], fields[3]] }

/*
  Make samples file list of crams optional by generating on the fly.
    Stop-gap until samples file is optional for RandomHetHom3
  Use first() to create a value channel that can be read unlimited times.
*/
if(params.samples_path == 'NO_FILE') {
  // Debugging using take(5)
  samples_file_chan = Channel.fromPath(cram_files)
                             .take(5)
                             .map{ "${it.getSimpleName()}\t${it.getSimpleName()}" }
                             .collectFile(name: "generated_samples_list.tsv", newLine: true)
                             .first()
} else {
  samples_file_chan = Channel.fromPath(samples_file)
                             .first()
}

process variants_by_sample {
  scratch true
  errorStrategy "retry"
  maxRetries 3

  input:
  file vcf from Channel.fromPath(params.vcfs_path)
  file samples_file from samples_file_chan

  output:
  file "*.variant_map.tsv.gz" optional true into variant_maps
  file "*.sample_map.tsv.gz" optional true into sample_maps

  """
  vcf_hash=`basename ${vcf} | md5sum | cut -f1 -d" "`
  ${params.random.exec} -i ${vcf} -s ${samples_file} -k ${params.random.max_hethom} -e ${params.random.seed} -o \${vcf_hash}
  """
}


/* Two processes: Concat sample maps together, Create sample map offsets file

    Get header from first sample map
    Strip headers from all sample maps
    concat all into sort
      sort, compress
    bgzip and index.

    grep non-comment lines, emit byte offset
    awk splitting on ':' ie. byte_offset:content_matched
    track first and last byte of instance of content_matched

    output:
    content\tbegin\tlength
*/
process merge_sample_maps {
  input:
  file "sample_map_*" from sample_maps.collect()

  output:
  tuple file("sample_map.tsv.gz"), file("sample_map.tsv.gz.gzi"), file("sample_map.tsv.gz.offsets") into sample_map 

  script:
  """
  {
     gzip -dc sample_map_1 | grep "^#";
     for f in sample_map_*; do
        gzip -dc \${f} | grep -v "^#" || true;
     done | sort -k 1,1 -k2,2 -k3n,3n --compress-program=gzip --temporary-directory=`pwd`;
  } | bgzip -c -i -I sample_map.tsv.gz.gzi > sample_map.tsv.gz;

  bgzip -dc sample_map.tsv.gz |\
    grep -bPo "^[^\t#]*" |\
    awk -F":" '{if (!(\$2 in beg)) beg[\$2] = \$1; end[\$2] = \$1;}\
               END { for (k in beg) {print k, beg[k], end[k] - beg[k]}}' > sample_map.tsv.gz.offsets
  """
}


process merge_variant_maps {
  echo true   

  input:
  file "variant_map_*" from variant_maps.collect()

  output:
  file "variant_map.tsv.gz*"

  publishDir "result/", pattern: "variant_map.tsv.gz*"

  script:
  """
  {
     gzip -dc variant_map_1 | grep "^#";
     for f in variant_map_*; do
        gzip -dc \${f} | grep -v "^#" || true;
     done | sort -k 1,1 -k2n,2n --compress-program=gzip --temporary-directory=`pwd`; 
  } | bgzip -c > variant_map.tsv.gz;
  tabix -s1 -b2 -e2 variant_map.tsv.gz;
  """
}


process extract_sequences {
   scratch true
   errorStrategy "ignore"

   input:
   set val(new_sample_name), val(sample_name), val(cram), val(crai), file(sample_map), file(sample_map_gzi), file(sample_map_offsets) from crams.combine(sample_map)
   each file(reference) from Channel.fromPath(params.reference_path)

   output:
   file "${new_sample_name}.cram*" optional true into samples_crams

   publishDir "result/sequences/", pattern: "${new_sample_name}.cram*", saveAs: { filename -> "${new_sample_name.md5().take(2)}/${filename}" }

   """
   header_bytes=`sort -k2n,2n ${sample_map_offsets} | head -n1 | cut -f2 -d" "`
   header_bytes=\$((header_bytes * 2))
   offset=`grep "^${sample_name}\\s" ${sample_map_offsets} | cut -f2 -d" "`
   bytes=`grep "^${sample_name}\\s" ${sample_map_offsets} | cut -f3 -d" "`
   bytes=\$((bytes * 2))
   {
      bgzip -b 0 -s \${header_bytes} sample_map.tsv.gz | grep "^#";
      bgzip -b \${offset} -s \${bytes} ${sample_map} | grep "^${sample_name}\\s";
   } | bgzip -c > variants.tsv.gz
   ${params.sequences.exec} -c ${cram} -i ${crai} -v variants.tsv.gz -r ${reference} -w ${params.sequences.window} -o ${new_sample_name}.no_sort.cram
   ${params.samtools.exec} sort -m ${params.samtools.max_mem} -O cram --reference ${reference} -o ${new_sample_name}.cram ${new_sample_name}.no_sort.cram
   ${params.samtools.exec} index ${new_sample_name}.cram
   """
}

