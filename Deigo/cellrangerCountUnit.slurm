#!/bin/bash
#SBATCH -p compute
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --array=0-1
#SBATCH --output=cellranger_%j.out

# load the module
module load cellranger/6.0.0

ilist=(ID*)
input=${ilist[$SLURM_ARRAY_TASK_ID]}

# id is output folder name, here we name our output folder
output="output"
id="$input$output"

# fastqs is the path to folder contain your fastq files
fastqs=$input

# sample is the name of your sample e.g. if you have 〇〇_S1_L001_R1_001.fastq.gz, your sample name is 〇〇.
sample=$input

# transcriptome is the path to reference genome

cellranger count --id=$id \
                 --transcriptome=/bucket/IshikawaU/cellranger/6.0.0/refdata-gex-mm10-2020-A \
                 --fastqs=$fastqs \
                 --sample=$sample \
                 --localcores=8 \
                 --localmem=32
