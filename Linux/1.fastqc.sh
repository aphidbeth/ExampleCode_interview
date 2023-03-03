#!/bin/bash

#SBATCH --ntasks=5
#SBATCH --job-name=aphrad19_fastqc_raw
#SBATCH --mail-user=b.moore.17@abdn.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logfiles/aphrad19_fastqc_raw.out

# run fastqc on raw files 

module load fastqc/0.11.9

for file in *fastq.gz

do 
fastqc  ./raw/$file -o ./QC  --nogroup -t 5

done
