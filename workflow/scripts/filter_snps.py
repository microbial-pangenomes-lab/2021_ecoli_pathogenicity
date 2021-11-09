#!/usr/bin/env python


import sys
import argparse
import itertools
import numpy as np
import pandas as pd


def get_options():
    description = 'Filter and annotate SNP table'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('pangenome',
                        help='Panaroo csv output (ending in "_roary.csv")')
    parser.add_argument('snps',
                        help='SNPs TSV file from snippy')
    parser.add_argument('focal',
                        help='Focal strain')

    parser.add_argument('--reference',
                        default='K-12',
                        help='Reference strain (default: %(default)s)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    g = pd.read_csv(options.pangenome, index_col=0, sep=',')

    g = g[g.columns[13:]]

    gn = {}
    gg = {}
    for idx, (gg1, gg2) in g[[options.focal, options.reference]].iterrows():
        if str(gg1) == 'nan':
            continue
        for g1 in gg1.split(';'):
            gn[g1] = idx
        if str(gg2) == 'nan':
            continue
        for g1, g2 in itertools.product(gg1.split(';'),
                                        gg2.split(';')):
            if str(g2) != 'nan':
                gg[g1] = g2

    m = pd.read_csv(options.snps, sep='\t')

    m = m[m['EFFECT'].str.contains('missense')].copy()

    m['K-12'] = [gg[x] if x in gg
                 else np.nan
                 for x in m['LOCUS_TAG'].values]

    m['ALT-GENE'] = [gn[x] if x in gn
                     else np.nan
                     for x in m['LOCUS_TAG'].values]

    m = m[['CHROM', 'POS', 'REF', 'ALT', 'STRAND',
           'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'K-12',
           'ALT-GENE', 'PRODUCT',]]

    m.to_csv(sys.stdout, sep='\t', index=False)
