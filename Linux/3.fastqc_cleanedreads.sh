#!/bin/bash

#SBATCH --ntasks=5
#SBATCH --job-name=aphrad19_fastqc_cleaned
#SBATCH --mail-user=b.moore.17@abdn.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logfiles/fastqc_cleanedreads.out


# set working directory to file location
cd ~/sharedscratch/AphRad19/

# run fastqc on raw files 

module load fastqc/0.11.9

for sample in $(ls ./cleaned)


do

fastqc  ./cleaned/$sample -o ./QC/demultiplexed/ --nogroup -t 5

done
