#! /bin/sh
#
# calculate_seq_identity_for_peaks.sh
# Copyright (C) 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.
#


script="../script/calculate_seq_identity.py"


genome="../../../biodata/genomes/GRCh38/GRCh38_no_alt_plus_hs38d1_analysis_set/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.genome"
fa="../../../biodata/genomes/GRCh38/GRCh38_no_alt_plus_hs38d1_analysis_set/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"

lms=(L1 L2 TA TC)

lengths=(4 10)

for len in ${lengths[@]}
do
  for lm in ${lms[@]}
  do
    python $script "../data/oligos_${lm}_${len}.fa" "../analysis/peaks_${lm}_slopped_${len}.fa" > "../analysis/peaks_${lm}_slopped_${len}_iden.tsv"
  done
done
