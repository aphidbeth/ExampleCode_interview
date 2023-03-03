#!/bin/bash

#$SBATCH --output= ./logfiles/MultiQC.out # STDOUT  
#$SBATCH --ntasks-per-node=2
#$SBATCH --job-name=MultiQC 
#$SBATCH --mail-user=b.moore.17@abdn.ac.uk

# Using MultiQC to collate all of the fastqc files from the demultiplexed reads  

# set working directory to file location
cd ~/sharedscratch/AphRad19/

# Load multiqc
module load multiqc/1.12  

# Run multiQC 
multiqc QC/demultiplexed/


