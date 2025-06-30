#!/bin/bash

# Script à faire copier-coller puis exécuter

######################## Prétraitement fichier gff #################

cd ~/Bureau/yacouba_p2/
grep -v "#" ref/PlasmoDB-62_Pfalciparum3D7.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > ref/PlasmoDB-62_Pfalciparum3D7.gff.gz
tabix -p gff ref/PlasmoDB-62_Pfalciparum3D7.gff.gz 

######################## Prétraitement fichier gff si utilisé #################
#grep -v "#" ref/PlasmoDB-62_Pfalciparum3D7.bed | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > ref/PlasmoDB-62_Pfalciparum3D7.bed.gz
#tabix -p bed myData.bed.gz

#########################################################################
mkdir -p Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/outvep/ 
mkdir -p Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/echantillon/outvep/ 
sudo chmod -R 777 ~/Bureau/yacouba_p2

# exécuter la ligne ci-dessous seule d'abord
sudo docker run -t -i -v $HOME/Bureau/yacouba_p2:/data ensemblorg/ensembl-vep

cd Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/

###################################################################################
############################### Annotation VEP #####################################

#*************************** groupe
for seq in `ls *.vcf`
do
vep -i "$seq" -o outvep/"${seq%.vcf}"_out62.vcf --vcf --force_overwrite --fasta ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum -e --gff ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz
done

#*************************** échantillon
cd Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/echantillon/
for seq in `ls *.vcf`
do
vep -i "$seq" -o outvep/"${seq%.vcf}"_out62.vcf --vcf --force_overwrite --fasta ../../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum -e --gff ../../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz
done
