import click
import csv
import itertools
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqUtils import GC
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
            score = [float(x) for x in line.rstrip().split("\t")]
            for i, a in enumerate(nct):
                mtx[a].append(score[i])
    if name != '':
        yield build(name, mtx)

def fasta_iter(istream):
    name = None
    sequence = ''
    for line in istream:
        if line.startswith('>'):
            if name is not None:
                yield(name, sequence)
            name = line.rstrip().replace('>', '')
            sequence = ''
        else:
            sequence += line.rstrip()
    if name is not None:
        yield (name, sequence)

def mutation_iter(istream):
    for line in istream:
        if line.startswith('#'):
            continue

        oligo, pos, name, ref, alt = line.rstrip().split()
        if ('N' in ref) or ('N' in alt):
            continue
        yield (oligo, int(pos), name, ref, alt.split(','))

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
@click.option('--reference', '-r', type=click.File('r'), required=True)
@click.option('--mutation', '-M', type=click.File('r'), default='-')
@click.option('--out', '-o', type=click.File('w'), default='-')
def fit_motif(motif, reference, mutation, out):
    click.echo("Loading reference sequences.", err=True)
    oligos = { name: seq for (name, seq) in progress(fasta_iter(reference)) }

    click.echo("Loading motifs matrix.", err=True)
    motifs = { _name(m) : m for m in progress(homer_parse(motif)) }

    reader = mutation_iter(mutation)
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['oligo', 'rsid', 'ref', 'alt', 'tf', 'ref_score', 'alt_score',
                     'score', 'flank_gc', 'core_gc'])

    click.echo("Progessing mutations.", err=True)
    for ((oligo, pos, rsid, ref, alts), tf) in progress(itertools.product(reader, motifs)):
        sequence = oligos[oligo]
        motif = motifs[tf]
        refseq = sequence[:pos] + ref + sequence[(pos+1):]
        refat, refscore = _score(motif, refseq)
        refgc = GC(flank(refseq, refat, refat + len(motif)))
        refcore_gc = GC(sequence[refat:(refat+len(motif))])
        for alt in alts:
            altseq = sequence[:pos] + alt + sequence[(pos+1):]
            altat, altscore = _score(motif, altseq)
            if altscore > refscore:
                flank_gc = GC(flank(altseq, altat, altat + len(motif)))
                core_gc  = GC(sequence[altat:(altat+len(motif))])
            else:
                flank_gc = refgc
                core_gc = refcore_gc

            writer.writerow([oligo, rsid, ref, alt, tf, refscore, altscore,
                             refscore-altscore, flank_gc, core_gc])

def _score(motif, seq):
    seq = Seq(seq, unambiguous_dna)
    pssm = motif.pssm
    fw_scores = pssm.calculate(seq)
    rv_scores = pssm.calculate(seq.reverse_complement())

    fw_index = max(range(len(fw_scores)), key=fw_scores.__getitem__)
    rv_index = max(range(len(rv_scores)), key=rv_scores.__getitem__)
    if fw_scores[fw_index] > rv_scores[rv_index]:
        return fw_index, fw_scores[fw_index]
    else:
        index = len(seq) - len(motif) - rv_index
        return index, rv_scores[rv_index]

def flank(seq, start, stop, size=10):
    bgupper = max([0, start-size])
    endown  = min([len(seq), stop+size])
    return seq[bgupper:start] + seq[stop:endown]

def _name(motif):
    name = motif.name.split('/')[0]
    return name.upper().replace('-', '_')

if __name__ == '__main__':
    fit_motif()
