# BRAVO Data Pipeline
Processing data to power the BRowse All Variants Online (BRAVO) API

1. Build, download, or install dependencies.
    1. Compile custom tools
    1. Install external tools
    1. Download external data
1. Collect data to be processed into convenient location.
1. Modify nextflow configs to match paths on your system.
1. Run nextflow workflows

## Input Data
**Naming:** The pipeline depends on the names of the input cram files having the sample ID as the first part of the filename.
Specifically, the expectation that the ID preceeds the first `.` such that a call to `getSimpleName()` yields the ID.

## Data Preparation Tools

### Compile Custom Tools
In the `data_prep/` directory you will find tools/scripts to prepare your data for importing into Mongo database and using in BRAVO browser.

```sh
cd data_prep/cpp_tools
cget install .
```
This build executables in `data_prep/cpp_tools/cget/bin`

### External Tools and Data
BamUtil, VEP, Loftee, and refernce data required is described in [dependencies.md](dependencies.md)

## Nextflow Scripts
In the `workflows/` directory are three Nextflow configs and scripts used to prepare the backing data for the BRAVO API.

Details about the steps of the pipeline are detailed in [data\_prep\_steps.md](data_prep_steps.md).

### Data Directory Contents Setup
Consolidating the results from the nextflow scripts into a single data directory as follows to power the BRAVO API.

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
- `variant_map.tsv.gz` is the result of
    

```sh
data/
├── cache
├── coverage
│   ├── bin_1
│   ├── bin_25e-2
│   ├── bin_50e-2
│   ├── bin_75e-2
│   └── full
├── external
├── crams
│   ├── sequences
│   ├── variant_map.tsv.gz
│   └── variant_map.tsv.gz.tbi
└── reference
    ├── hs38DH.fa
    └── hs38DH.fa.fai
```
