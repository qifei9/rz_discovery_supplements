#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2023 qifei <qifei@IOG.localdomain>
#
# Distributed under terms of the MIT license.

"""

"""

import argparse
from Bio import SeqIO


def calculate_identity(ref_seq, seq):
    """
    Returns the number of identical characters between two sequences.
    Assumes the sequences are aligned.
    """

    sa, sb, sl = ref_seq.seq, seq.seq, len(ref_seq.seq)
    matches = [sa[i] == sb[i] for i in range(sl)]
    # seq_identity = (100 * sum(matches)) / sl
    seq_identity = sum(matches)

    return seq_identity


def main():
    """main"""

    parse = argparse.ArgumentParser()

    parse.add_argument('ref_seq', help='reference sequence')

    parse.add_argument('seq_set', help='sequences to calculate identity to ref seq')

    args = parse.parse_args()

    ref_seq = SeqIO.read(args.ref_seq, "fasta")

    with open(args.seq_set) as f:
        for seq in SeqIO.parse(f, 'fasta'):
            identity = calculate_identity(ref_seq, seq)
            print(ref_seq.id, seq.id, identity, sep='\t')


if __name__ == "__main__":
    main()

