#!/usr/bin/env python


import os
import sys
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm


def get_options():
    description = 'Univariate analysis'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('phenotypes',
                        help='Phenotypes and covariants file')
    parser.add_argument('output',
                        help='Output directory for multivariate analysis input')
    
    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()
   
    m = pd.read_csv(options.phenotypes, sep='\t', index_col=0)
    # min - max scaling of age
    m['age'] = (m['age'] - m['age'].min()) / (m['age'].max() - m['age'].min())

    out = []
    targets = set(('deces', 'choc', 'passage_en_rea'))
    for target in targets:
        n = m.drop(columns=targets.difference((target, ))).copy()
        endog = n[target]
        keep = set([target, 'septicoli'])
        for c in n.columns.difference([target, 'septicoli']):
            exog = n[[c, 'septicoli']]
            exog = sm.add_constant(exog)
            mod = sm.Logit(endog, exog)
            res = mod.fit(disp=0, max_iter=100)
            t = res.summary(alpha=0.05)
            try:
                pval = res.pvalues[1]
                coef = np.power(2, res.params[1])
                # extract 95% CI
                lower_bound = np.power(2, float(t.tables[1].data[2][-2]))
                upper_bound = np.power(2, float(t.tables[1].data[2][-1]))
            except IndexError:
                pval = np.nan
                coef = np.nan
                lower_bound = np.nan
                upper_bound = np.nan
            if pval <= 0.1:
                keep.add(c)
            # save results: coefficient, 95% CI, p-value
            out.append(('full', target, c, coef, lower_bound, upper_bound, pval))
        # save resulting table
        n[sorted(keep)].rename(columns={target:
            'target'}).to_csv(os.path.join(options.output, f'{target}.tsv'),
                               sep='\t', index=False)
        for study, flag in zip(['septicoli', 'colibafi'],
                               [1, 0]):
            n = m[m['septicoli'] == flag].drop(
                    columns=targets.difference((target, )).union(('septicoli', ))).copy()
            endog = n[target]
            keep = set([target])
            for c in n.columns.difference([target]):
                exog = n[c]
                exog = sm.add_constant(exog)
                mod = sm.Logit(endog, exog)
                res = mod.fit(disp=0, max_iter=100)
                t = res.summary(alpha=0.05)
                try:
                    pval = res.pvalues[1]
                    coef = np.power(2, res.params[1])
                    # extract 95% CI
                    lower_bound = np.power(2, float(t.tables[1].data[2][-2]))
                    upper_bound = np.power(2, float(t.tables[1].data[2][-1]))
                except IndexError:
                    pval = np.nan
                    coef = np.nan
                    lower_bound = np.nan
                    upper_bound = np.nan
                if pval <= 0.1:
                    keep.add(c)
                # save results: coefficient, 95% CI, p-value
                out.append((study, target, c, coef, lower_bound, upper_bound, pval))
            # save resulting table
            n[sorted(keep)].rename(columns={target:
                'target'}).to_csv(os.path.join(options.output,
                                                f'{target}_{study}.tsv'),
                               sep='\t', index=False)

    r = pd.DataFrame(out, columns=['dataset', 'target', 'variable',
                                   'odds-ratio', 'odds-ratio-lower',
                                   'odds-ratio-higher', 'pvalue'])
    r.to_csv(sys.stdout, sep='\t', index=False)
