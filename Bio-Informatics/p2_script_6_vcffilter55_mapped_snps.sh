#!/bin/bash
cd ~/Bureau/yacouba_p2 
mkdir -p Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/echantillon/
cd Mapped/variant_calling_mapped/vcffilter55
cp -f *.vcf Gene/
cd Gene/

for seq in `ls *.vcf`
do
bgzip -c ${seq} > ${seq}.gz
tabix -p vcf ${seq}.gz 
done

############################################################################
############################### Fusion #####################################

#*************************** All vcf

bcftools merge Ank086S142_sort_mapped_dup_FILTERED55.vcf.gz Ank087S143_sort_mapped_dup_FILTERED55.vcf.gz ANK96S76_sort_mapped_dup_FILTERED55.vcf.gz ANK29S107_sort_mapped_dup_FILTERED55.vcf.gz ank123_sort_mapped_dup_FILTERED55.vcf.gz ANK67S102_sort_mapped_dup_FILTERED55.vcf.gz ANK33S100_sort_mapped_dup_FILTERED55.vcf.gz ANK28S106_sort_mapped_dup_FILTERED55.vcf.gz ANK55S101_sort_mapped_dup_FILTERED55.vcf.gz ANK54S108_sort_mapped_dup_FILTERED55.vcf.gz TC45CredepositS190_sort_mapped_dup_FILTERED55.vcf.gz TC41AredepositS193_sort_mapped_dup_FILTERED55.vcf.gz TC41BredepositS194_sort_mapped_dup_FILTERED55.vcf.gz TC45AredepositS188_sort_mapped_dup_FILTERED55.vcf.gz TC45BredepositS189_sort_mapped_dup_FILTERED55.vcf.gz ANK113S77_sort_mapped_dup_FILTERED55.vcf.gz YOP11S103_sort_mapped_dup_FILTERED55.vcf.gz ANK125S5_sort_mapped_dup_FILTERED55.vcf.gz YOP61S8_sort_mapped_dup_FILTERED55.vcf.gz YOP06S99_sort_mapped_dup_FILTERED55.vcf.gz YOP57S78_sort_mapped_dup_FILTERED55.vcf.gz YOP63S9_sort_mapped_dup_FILTERED55.vcf.gz YOP07S98_sort_mapped_dup_FILTERED55.vcf.gz YOP02S98_sort_mapped_dup_FILTERED55.vcf.gz YOP53S6_sort_mapped_dup_FILTERED55.vcf.gz YOP32S99_sort_mapped_dup_FILTERED55.vcf.gz YOP55S7_sort_mapped_dup_FILTERED55.vcf.gz TC31BredepositS192_sort_mapped_dup_FILTERED55.vcf.gz TC31AredepositS191_sort_mapped_dup_FILTERED55.vcf.gz ANK91S12_sort_mapped_dup_FILTERED55.vcf.gz > fusion/merge30filter55vcf.vcf


#************************ symptomatique vcf (as ac ss)

bcftools merge TC45CredepositS190_sort_mapped_dup_FILTERED55.vcf.gz TC41AredepositS193_sort_mapped_dup_FILTERED55.vcf.gz TC41BredepositS194_sort_mapped_dup_FILTERED55.vcf.gz TC45AredepositS188_sort_mapped_dup_FILTERED55.vcf.gz TC45BredepositS189_sort_mapped_dup_FILTERED55.vcf.gz ANK113S77_sort_mapped_dup_FILTERED55.vcf.gz YOP11S103_sort_mapped_dup_FILTERED55.vcf.gz ANK125S5_sort_mapped_dup_FILTERED55.vcf.gz YOP61S8_sort_mapped_dup_FILTERED55.vcf.gz YOP06S99_sort_mapped_dup_FILTERED55.vcf.gz YOP57S78_sort_mapped_dup_FILTERED55.vcf.gz YOP63S9_sort_mapped_dup_FILTERED55.vcf.gz YOP07S98_sort_mapped_dup_FILTERED55.vcf.gz YOP02S98_sort_mapped_dup_FILTERED55.vcf.gz YOP53S6_sort_mapped_dup_FILTERED55.vcf.gz YOP32S99_sort_mapped_dup_FILTERED55.vcf.gz YOP55S7_sort_mapped_dup_FILTERED55.vcf.gz TC31BredepositS192_sort_mapped_dup_FILTERED55.vcf.gz TC31AredepositS191_sort_mapped_dup_FILTERED55.vcf.gz > fusion/symptomatique30filter55vcf.vcf

# *******************************************  AA
bcftools merge Ank086S142_sort_mapped_dup_FILTERED55.vcf.gz Ank087S143_sort_mapped_dup_FILTERED55.vcf.gz ANK96S76_sort_mapped_dup_FILTERED55.vcf.gz ANK29S107_sort_mapped_dup_FILTERED55.vcf.gz ank123_sort_mapped_dup_FILTERED55.vcf.gz ANK67S102_sort_mapped_dup_FILTERED55.vcf.gz ANK33S100_sort_mapped_dup_FILTERED55.vcf.gz ANK28S106_sort_mapped_dup_FILTERED55.vcf.gz ANK55S101_sort_mapped_dup_FILTERED55.vcf.gz ANK54S108_sort_mapped_dup_FILTERED55.vcf.gz > fusion/AA.vcf 

#*******************************************  SS
bcftools merge YOP06S99_sort_mapped_dup_FILTERED55.vcf.gz YOP57S78_sort_mapped_dup_FILTERED55.vcf.gz YOP63S9_sort_mapped_dup_FILTERED55.vcf.gz YOP07S98_sort_mapped_dup_FILTERED55.vcf.gz YOP02S98_sort_mapped_dup_FILTERED55.vcf.gz YOP53S6_sort_mapped_dup_FILTERED55.vcf.gz YOP32S99_sort_mapped_dup_FILTERED55.vcf.gz YOP55S7_sort_mapped_dup_FILTERED55.vcf.gz TC31BredepositS192_sort_mapped_dup_FILTERED55.vcf.gz TC31AredepositS191_sort_mapped_dup_FILTERED55.vcf.gz > fusion/SS.vcf 

#*******************************************  AS
bcftools merge TC45CredepositS190_sort_mapped_dup_FILTERED55.vcf.gz TC41AredepositS193_sort_mapped_dup_FILTERED55.vcf.gz TC41BredepositS194_sort_mapped_dup_FILTERED55.vcf.gz TC45AredepositS188_sort_mapped_dup_FILTERED55.vcf.gz TC45BredepositS189_sort_mapped_dup_FILTERED55.vcf.gz ANK113S77_sort_mapped_dup_FILTERED55.vcf.gz YOP11S103_sort_mapped_dup_FILTERED55.vcf.gz ANK125S5_sort_mapped_dup_FILTERED55.vcf.gz YOP61S8_sort_mapped_dup_FILTERED55.vcf.gz > fusion/AS.vcf

############################################################################
############################### Filtre snps ##################################

#*******************************************  Ecantillon extraction SNPs
for seq in `ls Ank086S142_sort_mapped_dup_FILTERED55.vcf.gz Ank087S143_sort_mapped_dup_FILTERED55.vcf.gz ANK96S76_sort_mapped_dup_FILTERED55.vcf.gz ANK29S107_sort_mapped_dup_FILTERED55.vcf.gz ank123_sort_mapped_dup_FILTERED55.vcf.gz ANK67S102_sort_mapped_dup_FILTERED55.vcf.gz ANK33S100_sort_mapped_dup_FILTERED55.vcf.gz ANK28S106_sort_mapped_dup_FILTERED55.vcf.gz ANK55S101_sort_mapped_dup_FILTERED55.vcf.gz ANK54S108_sort_mapped_dup_FILTERED55.vcf.gz TC45CredepositS190_sort_mapped_dup_FILTERED55.vcf.gz TC41AredepositS193_sort_mapped_dup_FILTERED55.vcf.gz TC41BredepositS194_sort_mapped_dup_FILTERED55.vcf.gz TC45AredepositS188_sort_mapped_dup_FILTERED55.vcf.gz TC45BredepositS189_sort_mapped_dup_FILTERED55.vcf.gz ANK113S77_sort_mapped_dup_FILTERED55.vcf.gz YOP11S103_sort_mapped_dup_FILTERED55.vcf.gz ANK125S5_sort_mapped_dup_FILTERED55.vcf.gz YOP61S8_sort_mapped_dup_FILTERED55.vcf.gz YOP06S99_sort_mapped_dup_FILTERED55.vcf.gz YOP57S78_sort_mapped_dup_FILTERED55.vcf.gz YOP63S9_sort_mapped_dup_FILTERED55.vcf.gz YOP07S98_sort_mapped_dup_FILTERED55.vcf.gz YOP02S98_sort_mapped_dup_FILTERED55.vcf.gz YOP53S6_sort_mapped_dup_FILTERED55.vcf.gz YOP32S99_sort_mapped_dup_FILTERED55.vcf.gz YOP55S7_sort_mapped_dup_FILTERED55.vcf.gz TC31BredepositS192_sort_mapped_dup_FILTERED55.vcf.gz TC31AredepositS191_sort_mapped_dup_FILTERED55.vcf.gz ANK91S12_sort_mapped_dup_FILTERED55.vcf.gz`
do
bcftools view -v snps ${seq} > fusion/bcftools_snp_ok/echantillon/bcf_snp_${seq%.gz}
done

#*******************************************  Groupe extraction SNPs

cd fusion

for seq in `ls *.vcf`
do
bgzip -c ${seq} > ${seq}.gz
tabix -p vcf ${seq}.gz 
bcftools view -v snps ${seq}.gz > bcftools_snp_ok/bcf_snp_${seq%.gz}
done

############################################ Comparaison vcf genotype vs AA	#################################

cd bcftools_snp_ok

bgzip -c bcf_snp_AA.vcf  > bcf_snp_AA.vcf.gz
bgzip -c bcf_snp_AC.vcf  > bcf_snp_AC.vcf.gz
bgzip -c bcf_snp_AS.vcf  > bcf_snp_AS.vcf.gz
bgzip -c bcf_snp_SS.vcf  > bcf_snp_SS.vcf.gz
bgzip -c bcf_snp_merge30filter55vcf.vcf > bcf_snp_merge30filter55vcf.vcf.gz
bgzip -c bcf_snp_symptomatique30filter55vcf.vcf > bcf_snp_symptomatique30filter55vcf.vcf.gz

tabix -p vcf bcf_snp_AA.vcf.gz 
tabix -p vcf bcf_snp_AC.vcf.gz 
tabix -p vcf bcf_snp_AS.vcf.gz 
tabix -p vcf bcf_snp_SS.vcf.gz 
tabix -p vcf bcf_snp_merge30filter55vcf.vcf.gz
tabix -p vcf bcf_snp_symptomatique30filter55vcf.vcf.gz

bcftools merge bcf_snp_AA.vcf.gz bcf_snp_AS.vcf.gz > bcf_snp_AA_vs_AS.vcf
bcftools merge bcf_snp_AA.vcf.gz bcf_snp_AC.vcf.gz > bcf_snp_AA_vs_AC.vcf
bcftools merge bcf_snp_AA.vcf.gz bcf_snp_SS.vcf.gz > bcf_snp_AA_vs_SS.vcf
bcftools merge bcf_snp_AS.vcf.gz bcf_snp_SS.vcf.gz > bcf_snp_AS_vs_SS.vcf
