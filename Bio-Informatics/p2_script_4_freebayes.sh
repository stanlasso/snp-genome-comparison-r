#!/bin/bash

#################################################################################################
##################################### Global variant calling ####################################

mkdir -p Global/variant_calling/all_variant_calling

cd Global/align_sort_markdup
# Lister les fichiers *.dup.bam, enlever l'extension .dup.bam, trier et enlever les doublons, et sauvegarder dans un fichier
ls *.dup.bam > ../../data_markdup_global.txt

# indextation des fichiers .dup.bam

for seq in `cat ../../data_markdup_global.txt` 
do
samtools index ${seq}
done


# freebayes
freebayes -f ../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --ploidy 2 --bam-list ../../data_markdup_global.txt  > ../variant_calling/all_variant_calling/p2_all_pf_sample.vcf

# vcffile one one
for seq in `cat ../../data_markdup_global.txt` 
do 
freebayes -f ../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --ploidy 2 ${seq} > ../variant_calling/"${seq%.dup.bam}_dup.vcf"
done
cd ../../

################################################################################################
##################################### Mapped variant calling ###################################

mkdir -p Mapped/variant_calling_mapped/all_variant_calling
cd Mapped/align_sort_mapped_markdup
# Lister les fichiers *.dup.bam, enlever l'extension .dup.bam, trier et enlever les doublons, et sauvegarder dans un fichier
ls *.dup.bam > ../../data_markdup_mapped.txt


# indextation des fichiers .dup.bam

for seq in `cat ../../data_markdup_mapped.txt` 
do
samtools index ${seq}
done
# freebayes
freebayes -f ../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --ploidy 2 --bam-list ../../data_markdup_mapped.txt  > ../variant_calling_mapped/all_variant_calling/p2_all_pf_sample_mapped.vcf

# vcffile one one
for seq in `cat ../../data_markdup_mapped.txt` 
do 
freebayes -f ../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --ploidy 2 ${seq} > ../variant_calling_mapped/"${seq%.dup.bam}_dup.vcf"
done

cd ../../
