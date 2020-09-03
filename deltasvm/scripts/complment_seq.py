#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: complment_seq.py
# Create Date: 2018-03-15 16:07:12
#########################################

import sys
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="input fasta")
    parser.add_argument("-s", "--seq", dest="seq", required=True, help="ref fasta")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output fasta")
    args = parser.parse_args()

    seqs = set()
    for seq_record in SeqIO.parse(args.seq, "fasta"):
        seqs.add(str(seq_record.seq).upper())

    with open(args.out, "w") as f:
        for seq_record in SeqIO.parse(args.fasta, "fasta"):
            seq = str(seq_record.seq).upper()
            if seq in seqs:
                continue
            f.write(">{0}\n{1}\n".format(seq_record.id, seq))

if __name__ == "__main__":
    sys.exit(main())

