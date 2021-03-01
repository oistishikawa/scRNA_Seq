#!/bin/bash
#SBATCH -p compute
#SBATCH -t 96:00:00
#SBATCH --mem=32G
#SBATCH -c 8

export PATH=/home/m/masato-hirota/tools/cellranger-4.0.0:$PATH

# id is output folder name
# transcriptome is the path to reference genome
# fastqs is the path to folder contain your fastq files
# sample is the name of your sample e.g. if you have 〇〇_S1_L001_R1_001.fastq.gz, your sample name is 〇〇.

cellranger count --id=D4.5_P14_Arm_1 \
                 --transcriptome=/bucket/IshikawaU/Masato/analysis/data/scRNAseq
                 --fastqs=/flash/IshikawaU/Masato/scRNAseq_Exhaustion/data/fastq
                 --sample=scRNA_D4.5_P14_Arm_1 \
                 --localcores=8 \
                 --localmem=32
