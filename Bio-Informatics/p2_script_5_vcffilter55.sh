#!/bin/bash

#################################################################################################
##################################### Global filter variant calling ####################################

mkdir -p Global/variant_calling/vcffilter55/all_variant_calling

cd Global/variant_calling

ls *.vcf > ../../data_global_vcf.txt

vcffilter -f "DP > 5 & QUAL > 5" all_variant_calling/p2_all_pf_sample.vcf > vcffilter55/all_variant_calling/p2_all_pf_sample_FILTERED55.vcf

for seq in `cat ../../data_global_vcf.txt`
do
vcffilter -f "DP > 5 & QUAL > 5" "${seq}" > vcffilter55/"${seq%.vcf}"_FILTERED55.vcf 

done

cd ../../



#################################################################################################
##################################### Mapped filter variant calling ####################################

mkdir -p Mapped/variant_calling_mapped/vcffilter55/all_variant_calling

cd Mapped/variant_calling_mapped/

ls *.vcf > ../../data_mapped_vcf.txt

vcffilter -f "DP > 5 & QUAL > 5" all_variant_calling/p2_all_pf_sample_mapped.vcf > vcffilter55/all_variant_calling/p2_all_pf_sample_mapped_FILTERED55.vcf

for seq in `cat ../../data_mapped_vcf.txt`
do
vcffilter -f "DP > 5 & QUAL > 5" "${seq}" > vcffilter55/"${seq%.vcf}"_FILTERED55.vcf 
done

cd ../../
