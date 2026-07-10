# rz_discovery_supplements

## Files

- `code` contains scripts used in the manuscript
  - `01_processing` contains scripts for processing the raw sequencing data
  - `02_analysis` contains scripts for performing the analysis
    - `script_hsa` for human
    - `script_mmu_osa_dre` for mouse, rice and zebrafish
    - `script_dhy_ani_ath_tgu` for *Aspergillus nidulans*, *Trichoderma guizhouense*, *Arabidopsis thaliana* and *Drosophila hydei*
- `demo` contains demo data to test the code
  - `01_processing` contains demo sequencing data to test the processing code
  - `02_analysis` contains demo data to test the analysis code

## System requirements

The scripts have been tested on CentOS Linux 7, with the following dependencies:

```
snakemake 6.8.0
R 4.1
python 3.11.3
biopython 1.81
bedtools 2.30.0
cutadapt 3.4
FastQC 0.11.9
MultiQC 1.14
bowtie2 2.4.4
samtools 1.12
```

When running on real data, a high-performance computing machine is recommended.

No non-standard hardware is required.

## Installation guide

Install the above dependencies from repository of your Linux system or using conda.

The install time should be very short with a good internet connection.

## Demo

Use the demo data in folder `demo` to test the code:

1. download the human reference genome (GRCh38) and snakemake wrappers.
1. modify the `config.yaml` and `run_pipelines.sh` in `code/01_processing`, to change the path to the genome and wrappers to where you put them, and only run pipelines for human (hsa).
1. run `code/01_processing/run_pipelines.sh`.
1. modify the scripts in `02_analysis/scripts_hsa`, changing the path to the genome to where you put them.
1. run `code/02_analysis/scripts_hsa/find_candidates.r`

The expected run time should be short, since the demo data is very small.

The `code/01_processing/run_pipelines.sh` should generate counts of mapped fragments to the genome and spike-ins.

The `code/02_analysis/scripts_hsa/find_candidates.r` should not generate results for genome since the demo data is too small, while it should generate results for the spike-ins.

## Instructions for use

To run on your real data:

1. download the corresponding reference genome and snakemake wrappers.
1. modify the `config.yaml` and `run_pipelines.sh` in `code/01_processing`, to change the corresponding paths; modify the Snakefiles and adapter sequences if necessary.
1. run `code/01_processing/run_pipelines.sh`.
1. make a `sample_info.xlsx` according to your sample naming scheme and path (use `demo/02_analysis/sample_info.xlsx` as a template).
1. make fasta files for the last 10 and 4 bases for your adapters (use those `demo/02_analysis` if no change).
1. modify the scripts for the corresponding species in `02_analysis`, changing the corresponding paths.
1. run `find_candidates.r` in `code/02_analysis` for corresponding species.
