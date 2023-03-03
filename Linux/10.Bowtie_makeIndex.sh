#!/bin/bash

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -o ./slurmOutput/BowtieIndexSA3.out # STDout

# Build an index with bowtie 

# Load modules 
module load bowtie2/2.4.2

homedir=/uoa/home/r03bm17/sharedscratch/AphRad19

bowtie2-build $homedir/genomes/sa3_genome/S.avenae-primary-V1.0.fna $homedir/genomes/bowtiealign/sa3 --threads 10