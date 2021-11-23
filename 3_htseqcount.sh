#!/bin/bash
#SBATCH --job-name=htseqc
#SBATCH --mail-user=cecilepereira@ufl.edu
#SBATCH --output out/htseqcR-%A-%a.out
#SBATCH --error err/htseqcR-%A-%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6Gb
#SBATCH -t 1-24:00:00
#SBATCH --array=1-12
#SBATCH --dependency=afterok:2310060

#Cecile Pereira

#set directories
DIRBAM=../TOPHAT_150_Trim_ncbi_26oct/
DIRBAMT=../HTSEQC_150_Trim_ncbi_26oct/
DIRREF=../1_genome/

mkdir ${DIRBAMT}

#extract the line $SLURM_ARRAY_TASK_ID from the file containing the names of the samples (one name by line)
#samplename=`sed -n "${SLURM_ARRAY_TASK_ID} p" list_samples_names.txt`
input='list_samples_names.txt'

#module load htseq
#module load samtools

#sort BAM file 
while IFS= read -r i;
	do
	samtools sort -n ${DIRBAM}/${i}/accepted_hits.bam -o ${DIRBAMT}/ah_${i}_sortedbyname > ${DIRBAMT}/ah_${i}_sortedbyname

	#count reads 
	#s: stranded: no
	#r sorted by gene name
	htseq-count -t gene -r name -i Name -f bam -s no ${DIRBAMT}/ah_${i}_sortedbyname ${DIRREF}/GCF_000005845.2_ASM584v2_genomic.gff > ${DIRBAMT}/ah_${i}_counts.txt
done < "$input"