#!/usr/bin/env python


import os
import argparse
import pandas as pd


def get_options():
    description = 'Prepare all files for pyseer'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('unitigs',
                        help='Unitigs input file (to get a list of samples)')
    parser.add_argument('phenotypes',
                        help='Phenotypes and covariants file')
    parser.add_argument('similarities',
                        help='Similarity square matrix')
    parser.add_argument('distances',
                        help='Distances square matrix')
    parser.add_argument('lineages',
                        help='Lineages table')
    parser.add_argument('output',
                        help='Output directory')

    parser.add_argument('--invert',
                        action='store_true',
                        default=False,
                        help='Invert the "deces" and "choc" columns')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    samples = set(pd.read_csv(options.unitigs, sep='\t', index_col=0).index)

    p = pd.read_csv(options.phenotypes, sep='\t', index_col=0)
    if options.invert:
        p['deces'] = [1 if x == 0
                      else 0
                      for x in p['deces'].values]
        p['choc'] = [1 if x == 0
                     else 0
                     for x in p['choc'].values]

    s = pd.read_csv(options.similarities, sep='\t', index_col=0)
    d = pd.read_csv(options.distances, sep='\t', index_col=0)
    l = pd.read_csv(options.lineages, sep='\t', index_col=0, header=None)

    shared = samples.intersection(p.index).intersection(s.index).intersection(d.index).intersection(l.index)

    print(f'genomes: {len(samples)}, shared: {len(shared)}')

    p.loc[shared].to_csv(os.path.join(options.output, 'phenotypes.tsv'), sep='\t')
    s.loc[shared, shared].to_csv(os.path.join(options.output, 'similarity.tsv'), sep='\t')
    d.loc[shared, shared].to_csv(os.path.join(options.output, 'distances.tsv'), sep='\t')
    l.loc[shared].to_csv(os.path.join(options.output, 'lineages.tsv'), header=False, sep='\t')
