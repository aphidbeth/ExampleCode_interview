#!/bin/bash

#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem 100G
#SBATCH --job-name=aphrad19_demultiplex1
#SBATCH --mail-user=b.moore.17@abdn.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logfiles/demultiplex_plate1.out

echo "Time: $(date)"

# Process radtags script pool 1 

# set working directory to file location 
cd ~/sharedscratch/AphRad19/

# load stacks 
module load stacks/2.53

# Run process radtags

process_radtags -p ./raw/plate1 \
		-b ./raw/plate1/plate1_barcodes.txt \
		--inline_null\
		-i gzfastq\
		--renz_1 pstI --renz_2 mseI \
		-c -r -q \
		-o ./cleaned/

 
