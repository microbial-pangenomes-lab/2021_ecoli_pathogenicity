#!/usr/bin/env python


import dendropy
import argparse
import pandas as pd
from Bio import Phylo
import seaborn as sns
from matplotlib import colors


itol = '''DATASET_COLORSTRIP
#In colored strip datasets, each ID is associated to a color box/strip and can have an optional label. Color can be specified in hexadecimal, RGB or RGBA notation. When using RGB or RGBA notation, you cannot use COMMA as the dataset separator

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
#SEPARATOR COMMA
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL label1

#dataset color (can be changed later)
COLOR #ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#If COLOR_BRANCHES is set to 1, branches of the tree will be colored according to the colors of the strips above the leaves.
#When all children of a node have the same color, it will be colored the same, ie. the color will propagate inwards towards the root.
COLOR_BRANCHES 0


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE Dataset_legend
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_SHAPES 1 1 2 2
#LEGEND_COLORS #ff0000 #00ff00 rgba(0,255,0,0.5) #0000ff
#LEGEND_LABELS value1 value2 value3 value4
#LEGEND_SHAPE_SCALES 1 1 0.5 1

#width of the colored strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the color strip 
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 0

#display or hide the dataset label above the colored strip
#SHOW_LABELS 1

#text label size factor
#SIZE_FACTOR 1

#text label rotation
#LABEL_ROTATION 0

#text label shift in pixels (positive or negative)
#LABEL_SHIFT 0
'''


data = '''#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages

#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA'''


def get_options():
    description = 'Prepare an iTOL dataset'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('tree',
                        help='Tree in newick format')
    parser.add_argument('phenotypes',
                        help='Phenotypes file')
    parser.add_argument('phylogroups',
                        help='Phylogroups file')
    parser.add_argument('lineages',
                        help='Lineages file')
    
    parser.add_argument('--remove-clades',
                        action='store_true',
                        default=False,
                        help='Remove samples that are not E. coli sensu strictu')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    # tree
    tree = Phylo.read(options.tree, 'newick')
    # remove two samples:
    for c in tree.get_terminals():
        if c.name == 'R1B13F1' or c.name == 'R1B1E2':
            tree.prune(c)
    Phylo.write(tree, 'tree.nwk', 'newick')

    # phenotypes
    p = pd.read_csv(options.phenotypes, sep='\t', index_col=0)
    psepti = p['septicoli'].to_dict()
    pdeces = p['deces'].to_dict()
    pchoc = p['choc'].to_dict()
    prea = p['passage_en_rea'].to_dict()
    puri = p['pe_urinaire'].to_dict()
    pdige = p['pe_digestive'].to_dict()

    study_color = {0: '#7fc97f',
                   1: '#beaed4'}
    red_color = {0: '#fff4f2',
                 1: '#840000'}
    orange_color = {0: '#fff4f2',
                    1: '#f97306'}
    purple_color = {0: '#fff4f2',
                    1: '#9e43a2'}
    pee_color = {0: '#ffffcb',
                 1: '#bf9005'}
    brown_color = {0: '#ffffcb',
                   1: '#653700'}

    lg = {x.rstrip().split(',')[0].split('/')[-1].split('.')[0]: x.rstrip().split(',')[1]
          for x in open(options.lineages)
          if x.split(',')[0] != 'Taxon'}
    lgc = {x: y
           for x,y in zip(set(lg.values()),
                          sns.color_palette('hsv',
                                            len(set(lg.values()))))}
    
    clades = {x.rstrip().split(',')[0]: x.rstrip().split(',')[1]
          for x in open(options.phylogroups)
          if x.split(',')[0] != 'Strain'}
    ph = {x.rstrip().split(',')[0]: (x.rstrip().split(',')[1]
          if 'clade' not in  x.rstrip().split(',')[1]
          else 'E.clade')
          for x in open(options.phylogroups)
          if x.split(',')[0] != 'Strain'}
    phc = {'A': '#7995cb',
           'B1': '#5fad3e',
           'B2': '#b33740',
           'C': '#8a58a4',
           'D': '#f8f457',
           'E': '#935c29',
           'F': '#f8b515',
           'G': '#1d26db',
           'E.clade': '#65ae9a',}
    phc = {'A': '#017ab4',
           'B1': '#009048',
           'B2': '#d22313',
           'C': '#a6d9ad',
           'D': '#cacb3b',
           'E': '#d34b90',
           'F': '#fccfa2',
           'G': '#f29f4e',
           'E.clade': '#6e7880',}

    if options.remove_clades:
        for c in tree.get_terminals():
            if ph.get(c.name) == 'E.clade' and clades.get(c.name) != 'cladeI':
                tree.prune(c)
        Phylo.write(tree, 'tree.nwk', 'newick')
        node = 'R1B5C5'
        tree_d = dendropy.Tree.get(file=open('tree.nwk'), schema='newick')
        outgroup_node = tree_d.find_node_with_taxon_label(node)
        tree_d.to_outgroup_position(outgroup_node, update_bipartitions=True)
        tree_d.ladderize()
        tree_d.write(path='rerooted.nwk', schema='newick')
    else:
        node = 'H1-003-0086-C-F'
        tree_d = dendropy.Tree.get(file=open('tree.nwk'), schema='newick')
        outgroup_node = tree_d.find_node_with_taxon_label(node)
        tree_d.to_outgroup_position(outgroup_node, update_bipartitions=True)
        tree_d.ladderize()
        tree_d.write(path='rerooted.nwk', schema='newick')

    # study
    f = open('itol_study.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = psepti[c.name]
        color = study_color.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # phylogroups
    phd = sorted({ph[c.name] for c in tree.get_terminals()})
    phd_colors = [phc[x] for x in phd]
    phd_shapes = ['1' for x in phd]

    phd = ' '.join(phd)
    phd_colors = ' '.join(phd_colors)
    phd_shapes = ' '.join(phd_shapes)

    legend=f'''
    LEGEND_TITLE Species
    LEGEND_COLORS {phd_colors}
    LEGEND_LABELS {phd}
    LEGEND_SHAPES {phd_shapes}
    '''

    f = open('itol_phylogroup.txt', 'w')
    f.write(itol + '\n')
    f.write(legend + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = ph[c.name]
        color = phc.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # lineages
    f = open('itol_lineages.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = lg[c.name]
        color = lgc.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # deces
    f = open('itol_deces.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = pdeces[c.name]
        color = red_color.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # choc
    f = open('itol_choc.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = pchoc[c.name]
        color = orange_color.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # ICU
    f = open('itol_icu.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = prea[c.name]
        color = purple_color.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # urinaire
    f = open('itol_urinaire.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = puri[c.name]
        color = pee_color.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
    
    # digestive
    f = open('itol_digestive.txt', 'w')
    f.write(itol + '\n')
    f.write(data + '\n')
    for c in tree.get_terminals():
        name = c.name
        phylogroup = pdige[c.name]
        color = brown_color.get(phylogroup)
        if color is None:
            color = '#ffffff'
        else:
            color = colors.to_hex(color)
        f.write(f'{name} {color} {phylogroup}\n')
    f.close()
