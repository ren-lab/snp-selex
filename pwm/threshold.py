import click
import csv
from Bio import motifs

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

def progress(iter, freq=100):
    count = 0
    rotator = 0
    label = ("|", "/", "-", "\\")
    for i in iter:
        yield(i)
        count += 1
        if count % freq == 0:
            rotator = (rotator + 1) % 4
            click.echo("[%s] %6d\r" % (label[rotator], count), nl=False, err=True)

@click.command()
@click.option('--motif', '-m', type=click.File('r'), required=True)
@click.option('--out', '-o', type=click.File('w'), default='-')
def compute_threshold(motif, out):
    click.echo("Loading motifs matrix.", err=True)
    motifs = { _name(m) : m for m in progress(homer_parse(motif)) }

    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['tf', 'threshold'])

    click.echo("Progessing motifs.", err=True)
    for tf in motifs:
        motif = motifs[tf]
        motif.background = 0.4 #GC content of human genome
        threshold = _threshold(motif)
        writer.writerow([tf, threshold])

"""
Edited by Yunjiang
Determine threshold at fpr 0.01
"""
def _threshold(motif):
    pssm = motif.pssm
    background = {'A':0.3,'C':0.2,'G':0.2,'T':0.3}
    distribution = pssm.distribution(background=background, precision=10**4)
    threshold = distribution.threshold_fpr(0.01)

    return threshold

def _name(motif):
    name = motif.name.split('/')[0]
    return name.upper().replace('-', '_')

if __name__ == '__main__':
    compute_threshold()
