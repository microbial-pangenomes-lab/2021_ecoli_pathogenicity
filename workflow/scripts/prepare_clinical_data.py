#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

if __name__ == "__main__":
    table, out1, out2, out3 = sys.argv[1:5]
    
    t = pd.read_csv(table, sep=';')
    t = t.replace(to_replace={'Oui': 1, 'Non':0})

    t['appropriate_antibiotic'] = 1
    t.loc[t[(t['from'] == 'septi') & np.isnan(t['delai_atb'])].index, 'appropriate_antibiotic'] = 0
    t.loc[t[(t['from'] == 'septi') & (t['delai_atb'] >= 2)].index, 'appropriate_antibiotic'] = 0
    t['septicoli'] = [1 if x == 'septi'
                      else 0 if x == 'coli'
                      else np.nan
                      for x in t['from'].values]
    t['female'] = [1 if x == 'Femme'
                   else 0 if x == 'Homme'
                   else np.nan
                   for x in t['sexe'].values]
    t['community'] = [1 if x == 'Communautaire'
                      else 0 if x == 'Liee aux soins'
                      else np.nan
                      for x in t['type_inf2'].values]
    t['plurimicrobial'] = [1 if x == 'Pluri'
                           else 0 if x == 'Mono'
                           else np.nan
                           for x in t['plurimicrobien'].values]
    t['source'] = [0 if x == 'Au domicile'
                   else 1 if x == 'Institution / soins de suite'
                   else 2 if x == 'Hosp. service aigue'
                   else np.nan
                   for x in t['lieu_pat2'].values]

    # drop the columns that were re-encoded
    t = t.drop(columns=['from', 'sexe', 'type_inf2',
                        'plurimicrobien', 'lieu_pat2',
                        'delai_atb'])
    # drop some columns that have too many missing values
    t = t.drop(columns=['transplantation', 'neutropenie', 'antibio',
                        'grossesse', 'bmi', 'localisation_secondaire'])
    # also remove the recoded "source" column
    t = t.drop(columns=['source'])
    
    t = t.set_index('code_rangement').reset_index()
    
    y = t.copy()
    for c, v in y.nunique().iteritems():
        if v == 2:
            y[c] = y[c].astype('Int64')

    w = t.dropna().copy()
    for c, v in w.nunique().iteritems():
        if v == 2:
            w[c] = w[c].astype(int)

    t = t[['code_rangement', 'deces', 'choc', 'immuno',
           'plurimicrobial', 'septicoli', 'appropriate_antibiotic'] +
          [x for x in t.columns
           if x.startswith('pe_')]]

    z = t.dropna().copy()
    for c, v in z.nunique().iteritems():
        if v == 2:
            z[c] = z[c].astype(int)

    z = z[['code_rangement', 'deces', 'choc'] + sorted(z.columns.difference(['code_rangement', 'deces', 'choc']))]
    w = w[['code_rangement', 'deces', 'choc'] + sorted(w.columns.difference(['code_rangement', 'deces', 'choc']))]
    
    w.to_csv(out1, sep='\t', index=False)
    z.to_csv(out2, sep='\t', index=False)
    y.to_csv(out3, sep='\t', index=False)
