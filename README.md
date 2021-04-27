# BRAVO Data Pipeline
Processing data to power the BRowse All Variants Online (BRAVO) API

1. Build, download, or install dependencies.
    1. Compile custom tools
    1. Install external tools
    1. Download external data
1. Collect data to be processed into convenient location.
1. Modify nextflow configs to match paths on your system.
1. Run nextflow workflows

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

### Recommended Data Setup
Consolidating the results from the nextflow scripts into a single data directory as follows to power the BRAVO API.

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
