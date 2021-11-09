#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
 
parser$add_argument("table",
                    help="Input table")
parser$add_argument("output",
                    help="Output table")

args <- parser$parse_args()

suppressPackageStartupMessages(library("MASS"))

t = read.table(args$table, header=T, sep='\t')
fit <- glm('target ~ .' , data=t,  'binomial')
step <- stepAIC( fit  , direction="backward")
ci <- confint(step)
sum <- summary(step)
coef <- sum$coefficients
coef <- data.frame(coef)
ci <- data.frame(ci)
coef['name'] = row.names(coef)
ci['name'] = row.names(ci)
df <- merge(coef, ci)
write.table(df, args$output, sep='\t', row.names=F, quote=F)
