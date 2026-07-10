#! /bin/sh
#
# calculate_seq_identity_for_peaks.sh
# Copyright (C) 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.
#


script="../script/calculate_seq_identity.py"


genome="../../../biodata/genomes/combined_GRCm39_GRCz11_IRGSP-1.0/combined_GRCm39_GRCz11_IRGSP-1.0.genome"
fa="../../../biodata/genomes/combined_GRCm39_GRCz11_IRGSP-1.0/combined_GRCm39_GRCz11_IRGSP-1.0.fa"

lms=(L1 L2 TA TC)

lengths=(4 10)

for len in ${lengths[@]}
do
  for lm in ${lms[@]}
  do
    python $script "../data/oligos_${lm}_${len}.fa" "../analysis/peaks_${lm}_slopped_${len}.fa" > "../analysis/peaks_${lm}_slopped_${len}_iden.tsv"
  done
done
