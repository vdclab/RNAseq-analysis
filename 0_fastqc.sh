#!/bin/bash

liste=../FASTQfiles/*.fastq

for i in $liste;
	do ../FastQC_app/fastqc $i -o ~psalguero/Proyectos/miRNA_Cecile/FastQC
done
