#!/bin/bash
set -e

# Get the data from the gzipped vcf and write to tsv
echo "Getting VCF data out of ${SAMPLE}.vcf.gz..."
zcat $SAMPLE.vcf.gz | grep -v "#" > /tmp/vcf.tsv

# Run the R script to count genotypes by chromosome
echo "Counting genotypes on chromosomes..."
Rscript /src/pipeline.R

# Copy the R script's output from /tmp to /data in the container.
# if some directory on the host is mapped to /data on the container, 
# The output will be written here.
echo "Writing output to /data/${SAMPLE}.out.csv..."
cp /tmp/out.csv /data/${SAMPLE}.out.csv