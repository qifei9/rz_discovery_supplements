#! /bin/sh
#
# get_proceeding_seq_for_peaks.sh
# Copyright (C) 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.
#


genome="../../../biodata/genomes/combined_TAIR10.1_DhydRS2_ASM1142v1_ASM202278v1/combined_TAIR10.1_DhydRS2_ASM1142v1_ASM202278v1.genome"
fa="../../../biodata/genomes/combined_TAIR10.1_DhydRS2_ASM1142v1_ASM202278v1/combined_TAIR10.1_DhydRS2_ASM1142v1_ASM202278v1.fa"

lms=(L1 TA TC)

lengths=(4 10)

for len in ${lengths[@]}
do
  for lm in ${lms[@]}
  do
    bedtools slop -s -l $[len-1] -r 0 -i "../analysis/peaks_$lm.bed" -g $genome > "../analysis/peaks_${lm}_slopped_${len}.bed"
    bedtools getfasta -s -nameOnly -fi $fa -bed "../analysis/peaks_${lm}_slopped_${len}.bed" > "../analysis/peaks_${lm}_slopped_${len}.fa"
  done
done
