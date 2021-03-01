#!/bin/bash
#SBATCH -p compute
#SBATCH -t 10:00:00
#SBATCH --mem=32G
#SBATCH -c 8

export PATH=/home/m/masato-hirota/tools/cellranger-4.0.0:$PATH

# please access cellranger website https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate to get further information
# You have to create Aggregation.csv file to integrate data
# id is just output folder name you can change this freely

cellranger aggr --id=Gulfiya \
                --csv=/bucket/IshikawaU/Masato/Gulfiya_RNA/Aggregation.csv \
                --normalize=mapped

