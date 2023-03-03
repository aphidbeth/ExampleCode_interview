#!/bin/bash

#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem 100G
#SBATCH -o ./slurmOutput/GstacksRefMapSA3.out # STDout
#SBATCH --mail-user=b.moore.17@abdn.ac.uk
#SBATCH --mail-type=END,FAIL

# Setting working directory

homedir=/uoa/home/r03bm17/sharedscratch/AphRad19
cd $homedir

# Load modules
module load stacks/2.53
module load plink/1.09b6.26

################################################
# Steps for reference aligned Loci in stacks
################################################

# Step 1) Run script 10.Bowtie_makeIndex.sh (build the index for bowtie) DONE
#---------


# Step 2) Run script 11.1_BowtieArray.sh (align the cleaned reads against our index) -> BAM files in folder refmap/bowtie_align_out_sa3 DONE
#--------


# Step 3) Import into stacks and run gstacks to create a catalog of loci:
#--------
gstacks -I $homedir/refmap/bowtie_align_out_sa3/ -M $homedir/popmaps/popmap1.txt --unpaired -O $homedir/refmap/stacks_out_sa3/ -t 10


# Step 4) Run populations to get plink output files 
#---------
populations -P $homedir/refmap/stacks_out_sa3/ -M $homedir/popmaps/popmap1.txt -O $homedir/refmap/popRes_unfiltered/  -R 0.8 --plink -t 10> ./slurmOutput/log.pop


# Step 5) Use plink to remove missing data. mind = 0.5 missingness per ID,  geno = 0.2 missingness per loci.
#--------
echo "Use plink for filtering"

plink --file $homedir/refmap/PopsPlinkFiltered_sa3/populations.plink  --allow-extra-chr --mind 0.5 --geno 0.20 \
  --out $homedir/refmap/PopsPlinkFiltered_sa3/cleaned --recode &> $homedir/slurmOutput/PopsPlinkFilterp2_sa3.oe

echo "Create list of loci that pass filters"
cut -f 2 $homedir/refmap/PopsPlinkFiltered_sa3/cleaned.map | sed 's/_/\t/' > ./Strict_whitelist.txt

echo "Make a list of the IDs to remove"
cut -f 2 $homedir/refmap/PopsPlinkFiltered_sa3/cleaned.irem > $homedir/refmap/PopsPlinkFiltered_sa3/ids_to_remove.txt

echo "remove IDs from original popmap and save as new popmap"
grep -F -v -f  $homedir/refmap/PopsPlinkFiltered_sa3/ids_to_remove.txt $POPMAP > $homedir/popmaps/plinkFiltered_popmapsa_sa3.txt


# Step 6) Rerun populations- at this stage use whitelist & new popmap - to output files for downstream analysis
#--------
populations -P $src/refmap/stacks_out_sa3/ -M $src/popmaps/plinkFiltered_popmapsa_sa3 -O $src/refmap/populations_output/popRes_filtered/ ./Strict_whitelist.txt \
				--write-single-snp --plink --structure --phylip --vcf -t 10 > ./slurmOutput/log.pop_rerun


# Step 7) downstream analysis => Rscripts
#--------





