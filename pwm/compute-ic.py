import click
import csv
import numpy as np
from math import log
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import unambiguous_dna


def homer_parse(fstream):
    def build(name, freq):
        m = motifs.Motif(counts=freq)
        m.name = name
        return m
    
    nct  = "ACGT"
    name = ""
    mtx  = {a: [] for a in nct}

    for line in fstream:
        if line.startswith('>'):
            if name != '':
                yield build(name, mtx)
            name = line.rstrip().split()[1]
            mtx  = {a: [] for a in nct}
        else:
            score = [float(x) for x in line.rstrip().split()]
            for i, a in enumerate(nct):
                mtx[a].append(score[i])
    if name != '':
        yield build(name, mtx)


@click.command()
@click.option('--motif', '-m', type=click.File('r'), required=True)
@click.option('--out', '-o', type=click.File('w'), default='-')
def main(motif, out):
    writer = csv.DictWriter(out, delimiter='\t', fieldnames=['name', 'ic', 'mean_ic'])
    writer.writeheader()

    for m in homer_parse(motif):
        ic = _ic(m)
        writer.writerow(dict(name = _name(m), ic = ic, mean_ic = ic/m.length))


def _ic(motif):
    acc = 0

    pwm = motif.pwm
    # background = motif.background
    for i in range(pwm.length):
        for letter in "ACGT":
            p = pwm[letter, i]
            # b = 0.25 # Naive background
            # b = background[letter]
            acc += p * log(p, 2)
        acc += 2
    return acc

        
def _name(motif):
    name = motif.name.split('/')[0]
    return name.upper().replace('-', '_')


if __name__ == '__main__':
    main()
