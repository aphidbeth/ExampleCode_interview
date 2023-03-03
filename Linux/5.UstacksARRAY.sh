#!/bin/bash
#SBATCH -o ./slurmOutput/slurm.%x_%j.out # STDOUT
#SBATCH --job-name=UstacksArray
#SBATCH --mail-user=b.moore.17@abdn.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=10
#SBATCH --mem 100G
#SBATCH --time=2-00:00:00

output_dir=/uoa/home/r03bm17/sharedscratch/AphRad19/denovo/paramTesting

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

catalog="${array[0]}"
M="${array[1]}"
m="${array[2]}"
n="${array[3]}"
popmap="${array[4]}"
Samplename="${array[5]}"   
SampleID="${array[6]}" 	
	 
catalogDir=$output_dir/$catalog
     

if [ ! -d "$output_dir/$catalog" ]; then
        mkdir $output_dir/$catalog
		echo "Making $catalogDir directory"
     fi
	 
echo "running ustacks for catalog=$catalog, m=$m, M=$M, n=$n ID=$SampleID, SampleName"
	 
	 
	 
### RUN STACKS ###	 
	 
	 
# Loading stacks

module load stacks/2.53

# Setting working directory

cd /uoa/home/r03bm17/sharedscratch/AphRad19/

# Run ustacks #------------------------------------------------------------------------

catalog="${array[0]}"
M="${array[1]}"
m="${array[2]}"
n="${array[3]}"
popmap="${array[4]}"
Samplename="${array[5]}"   
SampleID="${array[6]}"
 
ustacks -f  cleaned/$Samplename -o $catalogDir -i 1 --name $Samplename \
-M $M -m $m -N 0 


echo $SampleID
#-------------------------------------------------------------------------------------



### RENAMING FILES ###

echo "Starting renaming"

for file in $(ls $catalogDir/${Samplename}*)
do
echo $file
gunzip $file
echo "gunzipped"
name=${file:0:-3}
echo "samp name extracted"
head -n 1 $name > ${name}.header
echo "header"
tail -n 1 $name > ${name}.footer
echo "footer"
cat $name | head -n -1 | tail -n +2 > ${name}.midsection
echo "midsection"
cut -f 2- ${name}.midsection> ${name}.midsectionmod
echo "midsection modified"
var=$(cat $name | wc -l | awk '{print $1-2}' | xargs)
echo "nrows counted"
yes $SampleID | head -n $var > ${name}.newIDcol
echo "new ID col created"
paste ${name}.newIDcol ${name}.midsectionmod >${name}.midsectionmod2
echo "new ID col added to midsection"
cat ${name}.header ${name}.midsectionmod2 ${name}.footer > ${name}
echo "header and footer added"
gzip ${name}
echo "gzipped"
echo $file "done"
rm ${name}.header 
rm ${name}.midsectionmod2 
rm ${name}.footer
rm ${name}.newIDcol
rm ${name}.midsection
rm ${name}.midsectionmod
echo "temporary files removed"
done

echo "renaming done" 
