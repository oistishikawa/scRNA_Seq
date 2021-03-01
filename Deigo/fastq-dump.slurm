#!/bin/bash
#SBATCH -p compute
#SBATCH -t 96:00:00
#SBATCH --mem=32G
#SBATCH -c 8

# After you download .sra file you can get fastq file using the code below

module load sra-tools/2.8.2-1

fastq-dump --split-files SRR8207939.sra

# Usually, we can get three files like test_sample1_S1_L001_I1_001.fastq, test_sample1_S1_L001_R1_001.fastq, test_sample1_S1_L001_R2_001.fastq
# However, sometimes fastq-dump generate only one file. then you can download oliginal bam file from ncbi database e.g. https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8207939
# Then you can convert bamfile to fastq using bamtofastq
