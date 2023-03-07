Explanation of files in the Linux folder

These are the bash scripts I wrote to process my raw reads from RAD-sequencing into commonly used output files for downstream analysis. 
I mainly use STACKS for this, but also call on fastqc, multiqc, bowtie, samtools and plink.

I outline briefly what each script does below: 
 
 - 1.fastqc.sh : Runs fastqc on the two sequencing plates
 - 2.1_clean_plate1.sh : uses stacks to remove reads with low quality and seperates reads into fastq file per sample
 - 2.2_clean_plate2.sh : same as 2.1_clean_plate1.sh but for the second plate
 - 3.fastqc_cleanedreads.sh : Run fastq on the cleaned fastq files. Allows identification of any abnormal samples that might have slipped through initial qc.
 - 4.Multiqc.sh : compile all the fastqc reports into an more easily readable html summary 
 - 5.Ustacks.ARRAY.sh : Runs the Ustacks module of STACKS on the raw reads which compiles them into putative loci. 
						The batch array allows the process to be run for a variety of parameter combinations (listed in UstacksParamTable.csv). 
							Please note I did not write lines 12-31 (these were provided by a collegue) which I adapted for my scripts.
Scripts 6-9 are missing from this folder as they are the remainder of the "denovo pipeline" which is still being worked on. 
- 10.Bowtie_makeIndex.sh : makes an index for the Sitobion avenae genome using default bowtie-2 settings.
- 11.1_BowtieArray.sh : Aligns reads to genome using bowtie-2 and sorts the reads with samtools and outputs BAM files. 
						Uses an array to speed up the process for each sample. The renaming/formatting of the files
						afterwards is necessary so the files can be read by STACKS. 
- 12.Refmap_Gstacks.sh : Runs the remainder of the STACKS pipeline on the BAM files outputted from 11.1_BowtieArray.sh. Includes plink filtering for missingness. 


