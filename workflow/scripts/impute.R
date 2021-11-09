#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
 
parser$add_argument("table",
                    help="Input table")
parser$add_argument("output",
                    help="Output table")
parser$add_argument("--iterations",
                    type="integer",
                    default=15, 
                    help="Rounds of imputation (default %(default)d)")

args <- parser$parse_args()

suppressPackageStartupMessages(library("mice"))

m = read.delim(args$table)
imp = mice(m, m=args$iterations, print=F)
dat = complete(imp)
write.table(dat, args$output, sep='\t', row.names=F, quote=F)
