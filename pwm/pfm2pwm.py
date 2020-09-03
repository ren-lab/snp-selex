#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: pfm2pwm.py
# Create Date: 2019-02-15 16:55:10
#########################################

import sys
import argparse
from Bio import motifs

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", dest="infile", required=True, help="input file")
    parser.add_argument("-n", "--name", dest="name", required=True, help="name of motif")
    args = parser.parse_args()

    name = args.name
    with open(args.infile) as handle:
        m = motifs.read(handle, "pfm")
        pwm = m.counts.normalize(pseudocounts={"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6})
        print(">{0}\t{1}".format(str(pwm.consensus), name))
        for i in range(len(pwm[0])):
            print("{0:f}\t{1:f}\t{2:f}\t{3:f}".format(pwm[0][i], pwm[1][i], pwm[2][i], pwm[3][i]))

if __name__ == "__main__":
    sys.exit(main())

