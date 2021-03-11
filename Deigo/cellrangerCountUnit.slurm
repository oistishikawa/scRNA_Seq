#!/bin/bash
#SBATCH -p compute
#SBATCH -t 96:00:00
#SBATCH --mem=32G
#SBATCH -c 8

module load cellranger/6.0.0
# id is output folder name
# transcriptome is the path to reference genome
# fastqs is the path to folder contain your fastq files
# sample is the name of your sample e.g. if you have 〇〇_S1_L001_R1_001.fastq.gz, your sample name is 〇〇.

#enter your output folder here
id=(your_output_folder)

#enter your scRNA data dir here
fastqs=(your_fastq)

#enter the prefix of your sample file
sample=(your_sample_prefix)

cellranger count --id=$id \
                 --transcriptome=/bucket/IshikawaU/refdata-gex-mm10-2020-A \
                 --fastqs=$fastqs \
                 --sample=$sample \
                 --localcores=8 \
                 --localmem=32
