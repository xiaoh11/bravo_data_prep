# Prepare VCF Pipeline Dependencies
Install dependencies and take note of where you put them.
The executables either need to be in the PATH of the executing environment,
or the full path to them will need to be supplied to the nextflow configs.

## [Nextflow](https://www.nextflow.io/)
Data processing workflows are written in Nextflow's DSL.

## External Data Sources

### Basis Data
Basis data is loaded in to backing database directly.
- Genenames: HUGO Gene Nomenclature Commitee (HGNC)
- Canonical Transcripts: Ensembl ID to Ensembl Transcript ID mapping
- OMIM File: Ensemble ID to OMIM ID mapping
- Gencode File: Source list of genes
See [basis\_data.md](basis_data.md) for where to obtain this data.

### CADD Scores
CADD scores that get added to the VCF backing data during `prepare_VCF` pipeline originate from [here](https://cadd.gs.washington.edu/)
The pregenerated "All possible SNVs of GRCh38/hg38" provided the backing data.

A pipeline to generate CADD scores locally is not in the scope of this project.

### Reference Sequence Data
```
curl ftp://share.sph.umich.edu/vt/grch38/hs38DH.fa -o hs38DH.fa
curl ftp://share.sph.umich.edu/vt/grch38/hs38DH.fa.fai -o hs38DH.fa.fai
```

## BamUtil
Clone and install from [bamUtil repo](https://github.com/statgen/bamUtil)

## VEP
Installation [instructions](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html) on ensembl.org.

## Loftee
Master branch of [LoF plugin](https://github.com/konradjk/loftee) doesn't work for GRCh38
per this [issue](https://github.com/konradjk/loftee/issues/73#issuecomment-733109901)
```
git clone --depth 1 --branch grch38 --single-branch git@github.com:konradjk/loftee.git
```

### Getting supporting data for loftee options
Running loftee in the prepare vcf workflow

```
  loftee_human_ancestor_fa = "/path/to/VEP/Plugins/loftee_data/human_ancestor.fa.gz"
  loftee_conservation_file = "/path/to/VEP/Plugins/loftee_data/loftee.sql"
  loftee_gerp_bigwig       = "/path/to/VEP/Plugins/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
```

#### Human ancestor sequence data at the time of writing.
```
mkdir loftee/data
curl https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz -o loftee/data/human_ancestor.fa.gz
curl https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai -o loftee/data/human_ancestor.fa.gz.fai
curl https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi -o loftee/data/human_ancestor.fa.gz.gzi

echo "f8c79d45c8fdffb52ef6926d540f2dd3  loftee/data/human_ancestor.fa.gz" | md5sum -c
echo "205a31051be5f1a312c31abf8a298ed7  loftee/data/human_ancestor.fa.gz.fai" | md5sum -c
echo "121343b868d5da87cc04646d55e806c3  loftee/data/human_ancestor.fa.gz.gzi" | md5sum -c
```

#### Conservation file
```
curl https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz -o loftee/data/phylocsf_gerp.sql.gz
gunzip loftee/data/phylocsf_gerp.sql.gz
```

GERP scores for GRCh38 were lifted over from GRCh37 see [tweet thread](https://twitter.com/konradjk/status/1093324906773786624)
Information about conservation scores on the [Ensembl site](https://ensembl.org/info/genome/compara/conservation_and_constrained.html)
```
curl ftp://ftp.ensembl.org/pub/current_compara/conservation_scores/90_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw -o loftee/data/gerp_conservation_scores.homo_sapiens.GRCh38.bw
```



