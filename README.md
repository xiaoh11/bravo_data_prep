# BRAVO Data Pipeline
Processing data to power the BRowse All Variants Online (BRAVO) API

1. Build, download, or install dependencies.
    1. Compile custom tools
    1. Install external tools
    1. Download external data
1. Collect data to be processed into convenient location.
1. Modify nextflow configs to match paths on your system or cluster.
1. Run nextflow workflows

## Input Data
**Naming:** The pipeline depends on the names of the input cram files having the sample ID as the first part of the filename.
Specifically, the expectation that the ID preceeds the first `.` such that a call to `getSimpleName()` yields the ID.

### Sequence Data
Source cram files.  Original sequences from which the variant calls were made.

### Variant calls
Source bcf files. Generated running the [topmed variant calling pipeline](https://github.com/statgen/topmed_variant_calling) 

## Data Preparation Tools

### Compile Custom Tools
In the `tools/` directory you will find tools/scripts to prepare your data for importing into Mongo database and using in BRAVO browser.

```sh
cd tools/cpp_tools
cget install .
```
This build executables in `tools/cpp_tools/cget/bin`

### External Tools
BamUtil, VEP, and Loftee tools required are described in [dependencies.md](dependencies.md)

### External Data
Gencode, Ensembl, dbSNP, and HUGO data required are described in [basis\_data.md](basis_data.md)

## Nextflow Scripts
In the `workflows/` directory are three Nextflow configs and scripts used to prepare the runtime data for the BRAVO API.

Details about the steps of the pipeline are detailed in [data\_prep\_steps.md](data_prep_steps.md).

## Data Directory Contents Setup

There are two data sets that Bravo API needs to run:
- *Backing Data* are flat files on disk read at runtime.
- *Basis Data* files processed and loaded into mongo db.

We consolidate the results from the nextflow scripts into a single data directory.
It powers the BRAVO API.
It contains all the runtime data on disk as well as the basis data that gets processed into the mongo db.

- `reference/` holds the refercence fasta files for the genome
- API's `SEQUENCE_DIR` config val is asking for directory that contains the 'sequences' directory.
  - sequences dirname is hardcoded
  - `variant_map.tsv.gz` file name is hardcoded.
  - `variant_map.tsv.gz.tbi` file name is hardcoded.
- Under sequence/, directory structure and filenames are perscribed.
  - All two hex character directories 00 to ff should exist as subdirectories.
  - cram files must have the filename in the exact form of `sample_id.cram`
  - The sub dir a cram belongs in is the first two characters of the md5 hexdigest of the sample_id.
    - E.g. foobar123.cram would be in directory "ae"
        ```python
        hashlib.md5("foobar123".encode()).hexdigest()[:2]
        ```
- coverage directory contents are taken from result/ dir of coverage workflow
- `variant_map.tsv.gz` is an output of `RandomHetHom3`
    

```sh
data/
├── cache
├── coverage
│   ├── bin_1
│   ├── bin_25e-2
│   ├── bin_50e-2
│   ├── bin_75e-2
│   └── full
├── crams
│   ├── sequences
│   ├── variant_map.tsv.gz
│   └── variant_map.tsv.gz.tbi
└── reference
    ├── hs38DH.fa
    └── hs38DH.fa.fai
```
