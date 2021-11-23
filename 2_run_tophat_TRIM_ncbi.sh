#!/bin/bash
#SBATCH --job-name=Mapping
#SBATCH --mem-per-cpu=4G
#SBATCH --time=40:00:00
#SBATCH --output out/maptrim-%A-%a.out
#SBATCH --error err/maptrim-%A-%a.err
#SBATCH --array=1-12
#SBATCH --qos=conesa-b
#SBATCH --nodes=1
#SBATCH --ntasks=4

#Cecile Pereira
#last modification 14 october 2016

#output directory
outdir=../TOPHAT_150_Trim_ncbi_27oct/
ref=../1_genome/GCF_000005845.2_ASM584v2_genomic.fna
datadir=../TRIM/

#extract the line $SLURM_ARRAY_TASK_ID from the file containing the names of the samples (one name by line)
input='list_samples_names.txt'

date;hostname;pwd

mkdir ${outdir}

#module load tophat

#distance between R1 and R2?
#number of mismatches: default 2, test 5
while IFS= read -r i;
	do
	mkdir ${outdir}/${i};
	./../tophat-2.1.1/tophat2 -r 150 -p 4 -o ${outdir}/${i}/ ${ref} ${datadir}${i}_1P ${datadir}${i}_2P; 

done < "$input"
