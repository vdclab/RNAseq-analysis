#!/bin/bash
#SBATCH --job-name=Trimm
#SBATCH --mem-per-cpu=4G
#SBATCH --time=40:00:00
#SBATCH --output out/trim-%A-%a.out
#SBATCH --error err/trim-%A-%a.err
#SBATCH --array=1-12
#SBATCH --qos=conesa-b
#SBATCH --nodes=1
#SBATCH --ntasks=4

#Cecile Pereira
#last modification 11 october 2016

#output directory
#outdir=/ufrc/conesa/cecilepereira/Leticia_tgt/TRIM/
#datadir=/ufrc/conesa/cecilepereira/Leticia_tgt/FASTQfiles/

outdir=../TRIM/
datadir=../FASTQfiles/

#extract the line $SLURM_ARRAY_TASK_ID from the file containing the names of the samples (one name by line)
#samplename=`sed -n "${SLURM_ARRAY_TASK_ID} p" list_samples_names.txt`
input='list_samples_names.txt'

date;hostname;pwd

mkdir ${outdir}

#module load tophat
#module load trimmomatic
while IFS= read -r i
	do java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 ${datadir}/${i}R1_001.fastq ${datadir}/${i}R2_001.fastq -baseout ${outdir}/${i} HEADCROP:9
done < "$input"
