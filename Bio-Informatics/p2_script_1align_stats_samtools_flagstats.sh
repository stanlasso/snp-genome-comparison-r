#!/bin/bash

# Dossier contenant les fichiers sam
sam_dir= pwd
# Fichier de sortie
output_file="/home/sherif/Bureau/yacouba_p2/Align_stats/statistiques_flagstat.tsv"

# En-tête du fichier tableau
#echo -e "Fichier\tTotal_reads\tMapped_reads\tMapped_percent" > "$output_file"

# Boucle sur chaque fichier sam
for sam in `ls *.sam`; do
    filename=$(basename "$sam")
    stats=$(samtools flagstat "$sam")
    
    total=$(echo "$stats" | head -n 1 | cut -d ' ' -f 1)
    mapped=$(echo "$stats" | grep "mapped (" | cut -d ' ' -f 1)
    percent=$(echo "$stats" | grep "mapped (" | sed -E 's/.*\(([^)]+)\).*/\1/')
    
    echo -e "${filename}\t${total}\t${mapped}\t${percent}" >> "$output_file"
done



####################################### samtools flagstats Details ###################
sam_dir= pwd
# Fichier de sortie
output_file="/home/sherif/Bureau/yacouba_p2/Align_stats/statistiques_flagstat_details.txt"

# Pour chaque fichier .sam dans le dossier
for sam in `ls *.sam`; do
    echo "Traitement de : $sam" | tee -a "$output_file"
    
    # Statistiques simples avec samtools flagstat
    samtools flagstat "$sam" >> "$output_file"
    
    echo -e "\n==============================\n" >> "$output_file"
done
