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

mkdir -p Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/outvep22_tab/
mkdir -p Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/echantillon/outvep22_tab/
sudo chmod -R 777 ~/Bureau/yacouba_p2

# exécuter la ligne ci-dessous seule d'abord
sudo docker run -t -i -v $HOME/Bureau/yacouba_p2:/data ensemblorg/ensembl-vep

cd Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/

###################################################################################
############################### Annotation VEP #####################################

#****************************************** Sortie tabulaire ****************************************************

#*************************** groupe

vep -i bcf_snp_AA.vcf -o outvep22_tab/bcf_snp_AA_out62.txt --tab --force_overwrite --fasta ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum --gff ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz -e

vep -i bcf_snp_AS.vcf -o outvep22_tab/bcf_snp_AS_out62.txt --tab --force_overwrite --fasta ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum --gff ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz -e

vep -i bcf_snp_SS.vcf -o outvep22_tab/bcf_snp_SS_out62.txt --tab --force_overwrite --fasta ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum --gff ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz -e


#*************************** échantillon

sudo docker run -t -i -v $HOME/Bureau/yacouba_p2:/data ensemblorg/ensembl-vep
cd Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/echantillon/
mkdir outvep22_tab

for seq in `ls *.vcf`
do
vep -i "$seq" -o outvep22_tab/"${seq%.vcf}"_out62.vcf --tab --force_overwrite --fasta ../../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum --gff ../../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz -e
done
