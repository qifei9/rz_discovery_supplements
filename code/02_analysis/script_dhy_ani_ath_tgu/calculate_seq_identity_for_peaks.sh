#! /bin/sh
#
# calculate_seq_identity_for_peaks.sh
# Copyright (C) 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.
#


script="../scripts/calculate_seq_identity.py"


lms=(L1 TA TC)

lengths=(4 10)

for len in ${lengths[@]}
do
  for lm in ${lms[@]}
  do
    python $script "../data/oligos_${lm}_${len}.fa" "../analysis/peaks_${lm}_slopped_${len}.fa" > "../analysis/peaks_${lm}_slopped_${len}_iden.tsv"
  done
done
