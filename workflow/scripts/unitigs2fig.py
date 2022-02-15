#!/usr/bin/env python
# Copyright (C) <2015> EMBL-European Bioinformatics Institute

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Neither the institution name nor the name roary_plots
# can be used to endorse or promote products derived from
# this software without prior written permission.
# For written permission, please contact <marco@ebi.ac.uk>.

# Products derived from this software may not be called roary_plots
# nor may roary_plots appear in their names without prior written
# permission of the developers. You should have received a copy
# of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.

__author__ = "Marco Galardini"
__version__ = '0.1.0'

def get_options():
    import argparse

    # create the top-level parser
    description = "Create plots from mapped unitigs"
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument('tree', action='store',
                        help='Newick Tree file')
    parser.add_argument('lineage', action='store',
                        help='Lineages file')
    parser.add_argument('spreadsheet', action='store',
                        help='Mapped unitigs file')
    parser.add_argument('base', action='store',
                        help='output base path')

    parser.add_argument('--labels', action='store_true',
                        default=False,
                        help='Add node labels to the tree (up to 10 chars)')
    
    parser.add_argument('--version', action='version',
                         version='%(prog)s '+__version__)

    return parser.parse_args()

if __name__ == "__main__":
    options = get_options()

    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import ListedColormap

    sns.set_style('white')

    import os
    import pandas as pd
    import numpy as np
    from Bio import Phylo

    t = Phylo.read(options.tree, 'newick')

    # Max distance to create better plots
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    phc = {'A': '#017ab4',
           'B1': '#009048',
           'B2': '#d22313',
           'C': '#a6d9ad',
           'D': '#cacb3b',
           'E': '#d34b90',
           'F': '#fccfa2',
           'G': '#f29f4e',
           'E.clade': '#6e7880',}

    l = pd.read_csv(options.lineage, sep='\t', header=None)
    l[2] = [x if 'clade' not in x else 'E.clade' for x in l[1].values]
    ld = {x: y for x, y in l[[0, 2]].values}
    ul = {x: i for i, x in enumerate(set(l[2].values))}
    uc = [phc[x] for x in sorted(ul, key=lambda x: ul[x])]
    
    umap = ListedColormap(uc)
    uvector = np.array([[ul[ld[x.name]]] for x in t.get_terminals()])

    m = pd.read_csv(options.spreadsheet, sep='\t')
    m['pres'] = 1
    m = m.pivot_table(index='unitig', columns='strain', values='pres')
    m = m.fillna(0)

    # add missing samples with no unitigs mapped
    for x in t.get_terminals():
        if x.name not in m.columns:
            m[x.name] = 0

    roary = m

    # Sort the matrix by the sum of strains presence
    idx = roary.sum(axis=1).sort_values(ascending=False).index
    roary_sorted = roary.loc[idx]

    # Sort the matrix according to tip labels in the tree
    roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]

    # Plot presence/absence matrix against the tree
    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))

        ax1=plt.subplot2grid((1,40), (0, 11), colspan=29)
        a=ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues,
                   vmin=0, vmax=1,
                   aspect='auto',
                   interpolation='none',
                   rasterized=True )
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.axis('off')

        ax = fig.add_subplot(1,3,1)
        ax=plt.subplot2grid((1,40), (0, 10), colspan=1)
        fig.subplots_adjust(wspace=0, hspace=0)
        b=ax.matshow(uvector, cmap=umap,
                   #vmin=0, vmax=1,
                   aspect='auto',
                   interpolation='none',
                   rasterized=True )
        ax.set_yticks([])
        ax.set_xticks([])
        ax.axis('off')
        
        ax2 = fig.add_subplot(1,3,2)
        # matplotlib v1/2 workaround
        try:
            ax2=plt.subplot2grid((1,40), (0, 0), colspan=10, facecolor='white')
        except AttributeError:
            ax2=plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('Unitigs presence/absence matrix\n(%d unitigs)'%roary.shape[0])

        if options.labels:
            fsize = 12 - 0.1*roary.shape[1]
            if fsize < 7:
                fsize = 7
            with plt.rc_context({'font.size': fsize}):
                Phylo.draw(t, axes=ax2, 
                           show_confidence=False,
                           label_func=lambda x: str(x)[:10],
                           xticks=([],), yticks=([],),
                           ylabel=('',), xlabel=('',),
                           xlim=(-mdist*0.1,mdist+mdist*0.45-mdist*roary.shape[1]*0.001),
                           axis=('off',),
                           title=('Tree\n(%d strains)'%roary.shape[1],), 
                           do_show=False,
                          )
        else:
            Phylo.draw(t, axes=ax2, 
                       show_confidence=False,
                       label_func=lambda x: None,
                       xticks=([],), yticks=([],),
                       ylabel=('',), xlabel=('',),
                       xlim=(-mdist*0.1,mdist+mdist*0.1),
                       axis=('off',),
                       title=('Tree\n(%d strains)'%roary.shape[1],),
                       do_show=False,
                      )
        plt.savefig(f'{options.base}', dpi=300)
        plt.clf()
