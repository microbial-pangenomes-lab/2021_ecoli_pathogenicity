#!/usr/bin/env python


import sys
import argparse
import numpy as np
import pandas as pd


def get_options():
    description = 'Compute variance/covariance matrix from all unitigs'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('unitigs',
                        help='Unitigs Rtab file (space delimited)')
    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.unitigs, sep=' ')

    #k = np.matmul(m.T.values, m.values)
    #k = pd.DataFrame(k, index=m.columns, columns=m.columns)
    #k.to_csv(sys.stdout, sep='\t')
    
    m.cov().to_csv(sys.stdout, sep='\t')
