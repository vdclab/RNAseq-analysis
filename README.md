# RNAseq-analysis
The repository has all the scripts used for the DE Analysis.

## list_samples_names.txt
List of the samples names

## 0_fastqc.sh
Run FasQC app to check the quality of the data.

## 1_run_trimmomatic.sh
Run trimmoatic to remove the adaptors of the data.

## 2_run_tophat_TRIM_ncbi.sh
Mapping the reads to the genome (before this step you have should created the genome index).

## 3_htseqcount.sh
Perform the read-counts.

## 4_tablecounts.sh and Extract_tablecounts_from_counts.pl
It uses the file "Extract_tablecounts_from_counts.pl" to extract the counts for each gene.

## DE_analysis.R
R Script with all the code to process the data and perform the DE Analysis.

Author: Cècile Pereira - cecile.pereira.bibs(at)gmail.com
