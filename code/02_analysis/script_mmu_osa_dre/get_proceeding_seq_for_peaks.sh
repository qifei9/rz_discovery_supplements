#! /bin/sh
#
# get_proceeding_seq_for_peaks.sh
# Copyright (C) 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.
#


genome="../../../biodata/genomes/combined_GRCm39_GRCz11_IRGSP-1.0/combined_GRCm39_GRCz11_IRGSP-1.0.genome"
fa="../../../biodata/genomes/combined_GRCm39_GRCz11_IRGSP-1.0/combined_GRCm39_GRCz11_IRGSP-1.0.fa"

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
