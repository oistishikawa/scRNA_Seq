#!/bin/bash
#SBATCH -p compute
#SBATCH -t 96:00:00
#SBATCH --mem=32G
#SBATCH -c 8

# if you want to analyze published data, you can find GEO accession number in the paper.
# Then you can serch it on NCBI database and find number like SRR*******
# You can download .sra file using prefetch

module load sra-tools/2.8.2-1

prefetch SRR8207938
prefetch SRR8207939
