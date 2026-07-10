#! /bin/sh
#
# get_proceeding_seq_for_peaks.sh
# Copyright (C) 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.
#


genome="../../../biodata/genomes/GRCh38/GRCh38_no_alt_plus_hs38d1_analysis_set/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.genome"
fa="../../../biodata/genomes/GRCh38/GRCh38_no_alt_plus_hs38d1_analysis_set/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"

lms=(L1 L2 TA TC)

lengths=(4 10)

for len in ${lengths[@]}
do
  for lm in ${lms[@]}
  do
    bedtools slop -s -l $len -r 0 -i "../analysis/peaks_$lm.bed" -g $genome > "../analysis/peaks_${lm}_slopped_${len}.bed"
    bedtools getfasta -s -nameOnly -fi $fa -bed "../analysis/peaks_${lm}_slopped_${len}.bed" > "../analysis/peaks_${lm}_slopped_${len}.fa"
  done
done
