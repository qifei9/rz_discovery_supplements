#! /bin/sh
#
# run_pipelines.sh
# Copyright (C) 2022 qifei9 <qifei9@gmail.com>
#
# Distributed under terms of the MIT license.
#

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_ani_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_ani_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_ani_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_ath_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_ath_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_ath_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dhy_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dhy_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dhy_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_tgu_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_tgu_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_tgu_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_hsa_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_hsa_L2
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_hsa_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_hsa_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_mmu_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_mmu_L2
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_mmu_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_mmu_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_osa_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_osa_L2
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_osa_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_osa_TC

snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dre_L1
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dre_L2
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dre_TA
snakemake -j 100 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_library_dre_TC

snakemake -j 20 --use-conda --wrapper-prefix 'git+file://path/to/snakemake-wrappers' -s ./Snakefile_QC
