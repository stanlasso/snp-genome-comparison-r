#!/bin/bash
cd align
# Lister les fichiers *.sam, enlever l'extension .sam, trier et enlever les doublons, et sauvegarder dans un fichier
ls *.sam | sed -e "s/\.sam//g" | sort | uniq > ../data_sam.txt

cd ..
mkdir -p Global/align_sort
mkdir -p Mapped/align_sort_mapped
for seq in `cat data_sam.txt` 
do
# trie des fichiers d'alignement  & stockage dans Global/align_sort/
picard SortSam I=align/${seq}.sam O=Global/align_sort/${seq}_sort.bam SORT_ORDER=coordinate
# extraire uniquement les lectures mappées (alignées sur la référence) & stockage dans Mapped/align_sort_mapped/
samtools view -b -F 4 Global/align_sort/${seq}_sort.bam > Mapped/align_sort_mapped/${seq}_sort_mapped.bam
done
