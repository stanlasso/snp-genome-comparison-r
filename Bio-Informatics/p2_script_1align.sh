#!/bin/bash
cd data
ls *.fastq | sed -e "s/_Human_umapped_R[12]\.fastq//g" | sort | uniq > ../dataseq.txt
cd ..
mkdir -p align
for seq in `cat dataseq.txt` 
do 
bwa mem -M -R "@RG\tID:${seq}\tSM:${seq}\tLB:${seq}\tPL:ILLUMIN" -o align/${seq}.sam ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta data/${seq}_Human_umapped_R1.fastq data/${seq}_Human_umapped_R2.fastq
done
