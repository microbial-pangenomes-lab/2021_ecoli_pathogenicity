#!/usr/bin/env python


import sys
import argparse
from Bio import SeqIO


def get_options():
    description = 'Fix mismatch between fasta and GFF'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('fasta',
                        help='Nucleotide fasta file')
    parser.add_argument('gff',
                        help='Annotation GFF file '
                             '(must have "##sequence-region" section)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    b = False
    contigs = []
    for l in open(options.gff):
        if not l.startswith('##sequence-region'):
            if b:
                break
            continue
        contig = l.rstrip().split()[1]
        contigs.append(contig)
        b = True

    seqs = [s for s in SeqIO.parse(options.fasta, 'fasta')]
    for name, s in zip(contigs, seqs):
        s.id = name
        s.description = ''

    SeqIO.write(seqs, sys.stdout, 'fasta')
