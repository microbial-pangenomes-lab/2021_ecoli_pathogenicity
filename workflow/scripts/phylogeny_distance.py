#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Extract a distance matrix from a phylogeny'''


def get_options():
    import argparse

    description = 'Extract a distance matrix from a phylogeny'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('phylogeny',
                        help='Tree file')

    parser.add_argument('--format',
                        default="newick",
                        help="Format of tree file [Default: newick]")
    method_group = parser.add_mutually_exclusive_group()
    method_group.add_argument('--calc-C',
                              action='store_true',
                              help='Produce var-covar matrix C (as from PDDIST). '
                                   'Always uses branch lengths.')
    method_group.add_argument('--topology',
                              action='store_true',
                              default=False,
                              help='Ignore branch lengths, and only use topological '
                                   'distances')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import sys
    import re
    import pandas as pd
    import dendropy

    sys.stderr.write('Loading tree\n')
    tree = dendropy.Tree.get(
        path=options.phylogeny,
        schema=options.format,
        preserve_underscores=True)

    d = {}
    sys.stderr.write('Generating distance matrix\n')
    pdm = tree.phylogenetic_distance_matrix()
    sys.stderr.write('Computing distances\n')
    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        d[taxon1.label] = d.get(taxon1.label, {})
        for taxon2 in tree.taxon_namespace:
            if taxon2.label not in d[taxon1.label].keys():
                if options.calc_C:
                    mrca = pdm.mrca(taxon1, taxon2)
                    d[taxon1.label][taxon2.label] = mrca.distance_from_root()
                elif options.topology:
                    d[taxon1.label][taxon2.label] = pdm.path_edge_count(taxon1, taxon2)
                else:
                    d[taxon1.label][taxon2.label] = pdm.patristic_distance(taxon1, taxon2)

    sys.stderr.write('Assembling final dataframe\n')
    m = pd.DataFrame(d)
    sys.stderr.write('Saving dataframe\n')
    m.to_csv(sys.stdout,
             sep='\t')

