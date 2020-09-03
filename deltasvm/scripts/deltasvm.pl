#!/usr/bin/perl
#
# deltasvm.pl: A simple script that calculates deltaSVM scores
#
# Copyright (C) 2015 Dongwon Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Last revised 03/30/2015
use strict;

if ($#ARGV != 3) {
    print "\n";
    print " Usage: deltasvm.pl <ref_seq_file> <alt_seq_file> <svm_weight_file> <output_file>\n\n";
    print "   ref_seq_file - reference sequences to score in FASTA format\n";
    print "   alt_seq_file - altenative sequences to score in FASTA format\n";
    print "   (NOTE: sequences should be in the same order, also have the same name)\n";
    print "   svm_weight_file - svm weights from gkm-SVM (i.e. all possible k-mer scores)\n";
    print "   output_file - name of the output file\n\n";
    exit(0);
}

my $refseqf = $ARGV[0];
my $altseqf = $ARGV[1];
my $svmwf = $ARGV[2];
my $outf = $ARGV[3];
my $kmerlen = 0;

my %rc;
$rc{"A"} = 'T'; $rc{"T"} = 'A'; $rc{"C"} = 'G'; $rc{"G"} = 'C';
sub revcomp {
    my $seq = shift @_;
    my $rcseq = "";
    my $l = length($seq)-1;
    while ($l >= 0) {
        $rcseq = $rcseq . $rc{substr($seq, $l, 1)};
        $l-=1;
    }

    return $rcseq;
}

my $line;

#2. read svmwfile
print STDERR "reading svmwfile $svmwf..\n";
my %svmscore;
open IN, "<$svmwf";
my $bias;
chomp($bias=<IN>);
while (chomp($line=<IN>)) {
    if ($line=~/^#/) {next;}
    if ($line=~/^bias/) {next;}

    my @f = split /\s+/, $line;
    if ($#f == 1) {
        $svmscore{$f[0]} = $f[1];
        $svmscore{revcomp($f[0])} = $f[1];
    } elsif ($#f == 2) {
        $svmscore{$f[0]} = $f[2];
        $svmscore{$f[1]} = $f[2];
    } else {
        print STDERR "[ERROR] unkown format of svm weight file. quit.\n";
        exit(0);
    }
    if ($kmerlen == 0) {
        $kmerlen = length($f[0]);
    }
}
close IN;

#3. scan sequence file
my $seqid_ref;
my $seqid_alt;
my $seq_ref;
my $seq_alt;
open IN1, "<$refseqf";
open IN2, "<$altseqf";
open OUT, ">$outf";
print STDERR "calculating deltaSVM using $refseqf and $altseqf..\n";
while (chomp($seqid_ref=<IN1>)) {
    $seqid_ref = substr($seqid_ref, 1);
    chomp($seqid_alt=<IN2>); 
    $seqid_alt = substr($seqid_alt, 1);

    chomp($seq_ref = uc(<IN1>) );
    chomp($seq_alt = uc(<IN2>) );

    if ($seqid_ref ne $seqid_alt) {
        print STDERR "[ERROR] the sequence id in $refseqf ($seqid_ref) is different from the one in $altseqf ($seqid_alt). skip.\n";
    }

    my $score = 0;
    my $seqlen = (length($seq_ref) > length($seq_alt)) ? length($seq_ref) : length($seq_alt);
    foreach my $i (0..($seqlen-$kmerlen)) {
        my $kmer_ref = substr($seq_ref, $i, $kmerlen);
        my $kmer_alt = substr($seq_alt, $i, $kmerlen);
        $score += ($svmscore{$kmer_alt} - $svmscore{$kmer_ref});
        #print "$kmer_ref\t$kmer_alt\t$svmscore{$kmer_ref}\t$svmscore{$kmer_alt}\n";
    }
    printf OUT "$seqid_ref\t%.6f\n", $score;
}
close OUT;
close IN1;
close IN2;
