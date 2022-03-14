# Details of BRAVO Data Preparation Processes

## Prepare VCF

We recommend to process each chromosome separately in parallel. You can further parallelize the process by specifying chromosomal regions in steps (2) and (3).

1. Compile data preparation tools:
   ```
   cd tools/cpp_tools/
   cget install .
   ```
   After successful compilation, the executables will be installed in `tools/cpp_tools/cget/bin`.
   
2. Preprare VCF with the INFO fields NS, AN, AC, AF, Hom, Het, DP, AVGDP, AVGDP_R, AVGGQ, AVGGQ_R, DP_HIST, DP_HIST_R, GQ_HIST, and GQ_HIST_R:
   ```
   ./cget/bin/ComputeAlleleCountsAndHistograms -i [input bcf/vcf] -s [samples file] -r [CHR:START-END] -o [output.vcf.gz]
   ```
   Input BCF/VCF must have DP and GQ FORMAT fields. Input BCF/VCF can be accessed both from local and from Google bucket storages. The input samples file (one sample ID per line) and chromosomal region CHR:START-END are optional.

3. Run [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep) on the VCF created in step (2):
   ```
   ./vep -i [input vcf.gz] --plugin LoF[,options]  --assembly [GRCh37/GRCh38] --cache --offline --vcf --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --pubmed --shift_hgvs 0 --allele_number --format vcf --force --buffer_size 100000 --compress_output gzip --no_stats -o [output vcf.gz]
   ```
   Specify [LoF plugin](https://github.com/konradjk/loftee) configuration options as you need.

4. (Optional) Obtain CADD scores from [https://cadd.gs.washington.edu](https://cadd.gs.washington.edu) and annotate VCF from step (3):
   ```
   python add_cadd_scores.py -i [input vcf.gz] -c [cadd_file1.tsv.gz] [cadd_file2.tsv.gz] ...  -o [output vcf.gz]
   ```
   CADD score files must by accompanied by the corresponding index files. If multiple CADD score files are specified, then the maximal CADD score across all files will be used.
<!-- 5. Now you are ready to import VCF's from step (4) into Mongo database. Index all input VCF files with `tabix` and run the following command:
   ```
   python manage.py variants -t [threads] -v [input chr1 vcf.gz] [input chr2 vcf.gz] ...
   ``` -->

## Prepare percentiles
Percentiles must be computed separately for each INFO field.

1. For each VCF INFO field (AVGDP, BQZ, CYZ, DP, FIBC_I, FIBC_P, HWE_SLP_I, HWE_SLP_P, IOR, NM0, NM1, NMZ, QUAL, STZ, SVM, ABE, ABZ) run:
   ```
   ./cget/bin/ComputePercentiles -i [input vcf.gz] -m [INFO field] -t [threads] -f [min MAF] -F [max MAF] -a [allele count] -p [number of perceniles] -d [description] -o [prefix for output files]
   ```
   Examples:
   ```
   ./cget/bin/ComputePercentiles -i /mymachine/myhome/mydata/chr*.mystudy.vcf.gz -t 10 -p 10 -o QUAL
   ./cget/bin/ComputePercentiles -i /mymachine/myhome/mydata/chr*.mystudy.vcf.gz -m ABE -t 10 -p 10 -d "Expected allele Balance towards Reference Allele on Heterozygous Sites" -o ABE
   ```

2. For each INFO field `X` in step (1), you will have two files `X.all_percentiles.json.gz` and `X.variant_percentile.vcf.gz`. The first is a compressed text file with INFO field description and percentiles in JSON format. The second is a compressed VCF file with `X_PCTL` INFO field which stores the corresponding percentile for every variant.
   
3. Index `X.variant_percentile.vcf.gz` using `tabix` and annotate your VCF files from previous step:
   ```
   find . -maxdepth 1 -name "*.variant_percentile..gz" -exec tabix {} \;
   python add_percentiles.py -i [input vcf.gz] - p QUAL.variant_percentile.vcf.gz ABE.variant_percentile.vcf.gz ... -o [output vcf.gz]
   ```

<!-- 3. Import `ALL.all_percentiles.gz` from step (2) into Mongo database:
    ```
    python manage.py metrics -m ALL.all_percentiles.gz
    ```
 4. Update Mongo database with variant percentiles from `*.variant_percentiles.gz` files:
    ```
    [will be added soon]
    ``` -->

## Prepare coverage

To prepare a coverage data for each base-pair position, you can use all your BAM/CRAM files or only a random subset of them (e.g. 1,000) if you need to reduce computational time.

1. For each chromosome and for each BAM/CRAM file extract depth per base-pair:
   ```
   samtools view -q 20 -F 0x0704 -uh [CRAM/BAM file] [chromosome] | samtools calmd -uAEr - [reference FASTA] | bam clipOverlap --in -.ubam --out -.ubam | samtools mpileup -f [reference FASTA] -Q 20 -t DP - | cut -f1-4 | bgzip > [chromosome].[sample].depth.gz
   ```
   In this step we use [clipOverlap from BamUtil](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap).
   
2. For each chromosome, create tabix index files for `[chromosome].[sample].depth.gz` e.g.:
   ```
   for f in 10.*.depth.gz; do tabix $f; done
   ```
 
3. For each chromosome, aggregate base-pair coverage acrosss output files `[chromosome].[sample].depth.gz` from step (1):
   ```
   python base_coverage/create_coverage.py -i [files list] aggregate -c [chromosome] -s [start bp] -e [end bp] | bgzip -c > [chromosome].[start].[end].json.gz
   ```
   The `files list` is a text file which lists all output files for a single chromosome from step (1). For example, you can create such list for chromosome 10 with `find . -name "10.*.depth.gz" -printf "%f\n" > depth_files.txt`.
   
   To generate chunks for each chromosome you can use the following command:
   ```
   python base_coverage/create_coverage.py -i [files list] chunk -c [chromosome] -s [chunk size in bp]
   ```
   A typical chunk size is 250,000 bp or 500,000 bp.
   
4. For each chromosome, merge files `[chromosome].[start].[end].json.bgz` from step (3):
   ```
   python base_coverage/merge_coverage.py -i [files list] -o [chromosome].full.json.gz
   ```
   The `files list` is a text file which lists all output files for a single chromosome from step (3).
 
5. After step (4), you should have coverage summary across your samples for each base pair in files `1.full.json.gz`, `2.full.json.gz`, ..., `22.full.json.gz`. For faster web-based visualization, you should prepare several pruned version of the coverage summary e.g.:
   ```
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 0.25 -o 22.bin_0.25.json.gz
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 0.50 -o 22.bin_0.50.json.gz
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 0.75 -o 22.bin_0.75.json.gz
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 1.00 -o 22.bin_1.00.json.gz
   ```
6. Tabix all coverage summary files.
7. Reference all of the coverage files in `BASE_COVERAGE` in `default.py`.

## Prepare CRAM

BRAVO uses IGV.js to visualize raw sequenced from up to 10 random alternate allele carriers. To enable this visualization BRAVO uses a pre-computed combined CRAM file with all reads from the carriers. We recommend to prepare a separate combined CRAM for each chromosome. The following steps describe hot to prepare the combined CRAM file:

1. In this step, for each variant you will extract IDs for 5 random heterozygous and 5 random homozygous alternate allele carriers from your original VCF file with genotypes. (optional) To speed-up the process, split each chromosome into chunks.
   ```
   RandomHetHom -k 5 -e 1987 --i [input vcf.gz] -s [samples file] -r [CHR:START-END] -o [output carriers_ids.vcf.gz]
   ```
   The value `1987` is a random seed, which you may want to change.
   
2. Prepare a text file (e.g. `samples.txt`) which lists BAM/CRAM file path for each sample:
   ```
   SAMPLEID1    /drive1/batch1/sampleid1.bam
   SAMPLEID2    /drive1/batch1/sampleid2.bam
   SAMPLEID3    /drive1/batch2/sampleid3.bam
   ...
   ```
3. In this step, you will create a combined CRAM file, by extracting all carriers' reads that overlap +/-100bp window around each variant.
   ```
   python prepare_sequences.py cram -i [carriers_ids.vcf.gz] -c samples.txt -w 100 -o [output combined.cram]
   ```

Note: sample IDs and read names in the combined CRAM files will be anonymized.
