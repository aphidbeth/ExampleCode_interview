#!/bin/bash

#SBATCH -o ./slurmOutput/slurm.BOWTIEALIGN_SA3%j.out # STDout
#SBATCH --job-name=UstacksArray
#SBATCH --mail-user=b.moore.17@abdn.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=10
#SBATCH --mem 100G
#SBATCH --time=2-00:00:00

homedir=/uoa/home/r03bm17/sharedscratch/AphRad19

while getopts ":t:p:b:" opt; do
  case ${opt} in
    t )
      filename=${OPTARG}
      ;;
  esac
done


echo batch table file $filename
if [ -z ${filename+x} ];
then
   echo "batch file is unset, use -t"
   exit -1
fi

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $filename)


IFS=', ' read -r -a array <<< "$LINE"

sample="${array[0]}"
file="${array[1]}"
echo $file
echo $sample

	 
### RUN BOWTIE ALIGN ###	 
	 
 
# Load modules 
module load bowtie2/2.4.2
module load samtools/1.14

#setwd
pwd
cd $homedir

echo "running bowtie align on sample=$sample"

# check the file assignment
echo $file
echo $sample

bowtie2 -x $homedir/genomes/bowtiealign/sa3 -q -U $homedir/cleaned/$file | samtools sort -o $homedir/refmap/bowtie_align_out_sa3/$sample.bam -@ 10 -

echo "$sample finished bowtie"
