#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: unique_fasta.py
# Create Date: 2018-03-15 15:19:46
#########################################

import sys
import argparse
from Bio import SeqIO
import operator

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="oligo fasta")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output fasta")
    parser.add_argument("-n", "--num", dest="num", type=int, default=0, help="num of seqs to output")
    args = parser.parse_args()

    seqs = {}
    for seq_record in SeqIO.parse(args.fasta, "fasta"):
        seq = str(seq_record.seq).upper()
        try:
            seqs[seq] += 1
        except KeyError:
            seqs[seq] = 1

    n_seq = 0
    with open(args.out, "w") as f:
        seqs_sorted = sorted(seqs, key=operator.itemgetter(1), reverse=True)
        if args.num > 0:
            seqs_sorted = seqs_sorted[0:args.num]
        for seq in seqs_sorted:
            f.write(">seq_{0:d}\n{1}\n".format(n_seq, seq))
            n_seq += 1

if __name__ == "__main__":
    sys.exit(main())

