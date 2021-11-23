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
DIRINP=../counts/
FILEOUT=../results/table_reads_counts.csv

DIROUT=../results/

mkdir $DIROUT

perl Extract_tablecounts_from_counts.pl	$DIRINP > $FILEOUT