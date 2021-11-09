#!/usr/bin/env python


import os
import sys
import argparse
import pandas as pd


def get_options():
    description = 'Prepare phenotype file from simulation'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('phenotypes',
                        help='Phenotypes file from BacGWASim')

    parser.add_argument('--sample',
                        default=None,
                        type=float,
                        help='Proportion of samples to randomly '
                             'select (default: all)')
    parser.add_argument('--seed',
                        default=42,
                        type=int,
                        help='Seed for random number generation'
                             '(default: %(default)d)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    p = pd.read_csv(options.phenotypes, sep=' ', header=None)
    p = p[[0, 2]]
    p.columns = ['sample', 'phenotype']
    p['phenotype'] = [0 if x == 1 else 1
                      for x in p['phenotype'].values]

    if options.sample is not None:
        p = p.sample(frac=options.sample, random_state=options.seed)

    p.to_csv(sys.stdout, sep='\t', index=False)
