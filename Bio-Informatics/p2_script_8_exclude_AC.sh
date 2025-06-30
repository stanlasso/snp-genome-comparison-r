cd /Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok
bcftools view -s "^ANK91S12" bcf_snp_merge30filter55vcf.vcf > bcf_snp_merge30filter55vcf_moins_AC.vcf

# Annotation
sudo chmod -R 777 ~/Bureau/yacouba_p2
sudo docker run -t -i -v $HOME/Bureau/yacouba_p2:/data ensemblorg/ensembl-vep

cd Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/

vep -i "bcf_snp_merge30filter55vcf_moins_AC.vcf" -o outvep/"bcf_snp_merge30filter55vcf_moins_AC_moins_AC"_out62.vcf --vcf --force_overwrite --fasta ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7_Genome.fasta --species plasmodium_falciparum -e --gff ../../../../../../ref/PlasmoDB-62_Pfalciparum3D7.gff.gz
