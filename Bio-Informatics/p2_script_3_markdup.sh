#!/bin/bash
cd Global/align_sort
# Lister les fichiers *.bam, enlever l'extension .bam, trier et enlever les doublons, et sauvegarder dans un fichier
ls *.bam | sed -e "s/\.bam//g" | sort | uniq > ../../data_sort_sam.txt

cd ../../
mkdir -p Global/align_sort_markdup
mkdir -p Mapped/align_sort_mapped_markdup
for seq in `cat data_sort_sam.txt` 
do 
# markduplicates des fichiers 
#picard MarkDuplicates I=Global/align_sort/${seq}.bam O=Global/align_sort_markdup/${seq}.dup.bam METRICS_FILE=Global/align_sort_markdup/metrics.${seq}.txt
# markduplicates des fichiers de lecture mappées
picard MarkDuplicates I=Mapped/align_sort_mapped/${seq}_mapped.bam O=Mapped/align_sort_mapped_markdup/${seq}_mapped.dup.bam METRICS_FILE=Mapped/align_sort_mapped_markdup/metrics.${seq}_mapped.txt
done

