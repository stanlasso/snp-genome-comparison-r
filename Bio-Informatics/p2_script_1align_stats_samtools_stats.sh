#!/bin/bash

sam_dir= pwd
output_file="/home/sherif/Bureau/yacouba_p2/Align_stats/statistiques_stats.tsv"

# En-tête
#echo -e "/home/sherif/Bureau/yacouba_p2/Align_stats" > "$output_file"

for sam in `ls *.sam`; do
    filename=$(basename "$sam")
    stats=$(samtools stats "$sam")

    total=$(echo "$stats" | grep ^SN | grep "raw total sequences" | cut -f 3)
    mapped=$(echo "$stats" | grep ^SN | grep "reads mapped:" | cut -f 3)
    qc_failed=$(echo "$stats" | grep ^SN | grep "reads QC failed" | cut -f 3)
    duplicates=$(echo "$stats" | grep ^SN | grep "reads duplicated" | cut -f 3)
    avg_length=$(echo "$stats" | grep ^SN | grep "average length" | cut -f 3)

    echo -e "${filename}\t${total}\t${mapped}\t${qc_failed}\t${duplicates}\t${avg_length}" >> "$output_file"
done

####################################### samtools Stats Details ###################
sam_dir= pwd
# Fichier de sortie
output_file="/home/sherif/Bureau/yacouba_p2/Align_stats/statistiques_stats_details.txt"

# Pour chaque fichier .sam dans le dossier
for sam in `ls *.sam`; do
    echo "Traitement de : $sam" | tee -a "$output_file"
    
    # Statistiques simples avec samtools flagstat
    samtools stats "$sam" | grep ^SN | cut -f 2- >> "$output_file"
    
    echo -e "\n==============================\n" >> "$output_file"
done
