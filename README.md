# Snakemake workflow: E. coli GWAS 2021

Snakemake workflow to reproduce the GWAS analysis on genetic determinants of
_E. coli_ bloodstream infections. Given the sensitive information used as 
covariates in this analysis, users wishing to replicate this study should
generate the data contained in the missing `data/GWAS_20201217_MEO.csv`
file. Summary statistics about each variable can be found as supplementary
material in the preprint (to be available soon).

## Input genomes

The input assemblies can be found under the following bioproject accessions:
[PRJEB39260](https://www.ebi.ac.uk/ena/browser/view/PRJEB39260) and
[PRJEB35745](https://www.ebi.ac.uk/ena/browser/view/PRJEB39260). The fasta
files should be placed inside the `data/fastas` directory, and annotated GFFs
should be placed inside the `data/gffs` directory.

## Usage

To run the three main GWAS analysis (naive, using covariates and burden test)
and the power analysis:

    snakemake -p annotate_summary map_summary_nc pyseer_rare pyseer_simulation --cores 36 --use-conda

Snakemake will install the appropriate packages for each step as conda environments.
 symbolic link to a directory containing the eggnog-mapper database should be
placed in `data/eggnog-mapper`, as well as a symbolic link to the unzipped fasta
file from uniref50 (`data/uniref50.fasta`). 

## Author

Marco Galardini (marco.galardini@twincore.de)

