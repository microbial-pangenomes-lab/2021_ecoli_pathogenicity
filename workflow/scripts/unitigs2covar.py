#!/usr/bin/env python


import sys
import math
import random
import argparse
import numpy as np
import pandas as pd


def get_options():
    description = 'Compute variance/covariance matrix from all unitigs'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('phenotypes',
                        help='Phenotypes file (tab delimited, to define which strains are there)')
    parser.add_argument('unitigs',
                        help='Unitigs txt file')

    parser.add_argument('--sample',
                        type=float,
                        default=1,
                        help='What fraction of unitigs to sample (default: %(default).2f)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.phenotypes, sep='\t', index_col=0)
    idx = set(m.index)
    all_vars = []
    i = 0
    j = 0
    for line_in in open(options.unitigs):
        var_name, strains = (line_in.split()[0],
                             line_in.rstrip().split('|')[1].lstrip().split())
        strains = {x.split(':')[0] for x in strains}
        correct = math.sqrt(len(var_name))
        d = {x: 1 / correct
             for x in strains
             if x in idx}
        for x in idx:
            if x not in d:
                d[x] = 0
        af = sum(d.values()) / len(idx)
        if af == 1:
            continue
        j += 1
        if options.sample < 1 and random.random() > options.sample:
            continue
        i += 1
        all_vars.append(pd.DataFrame([d[x] for x in idx if x in d],
                                  columns=[i],
                                  index=idx).T)
        if not i % 1000:
            sys.stderr.write(f'read {i} unitigs ({j} eligible)\n')
    m = pd.concat(all_vars)
    del all_vars
        
    m.cov().to_csv(sys.stdout, sep='\t')
