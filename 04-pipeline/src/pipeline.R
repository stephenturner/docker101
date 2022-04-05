#!/usr/bin/env Rscript

# Read input /tmp from pipeline.sh
x <- read.table("/tmp/vcf.tsv")

# Select only certain columns and rename them using base R
x <- x[,c(1,2,4,5,10)]
names(x) <- c("chr", "pos", "ref", "alt", "gt")

# get genotypes by chromosome and write to container's /tmp/out.csv
table(x$chr, x$gt) |> 
    addmargins() |> 
    write.csv("/tmp/out.csv")