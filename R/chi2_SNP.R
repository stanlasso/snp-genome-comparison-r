# Charger les packages nécessaires
#install.packages("vcfR")
library(vcfR)
library(genetics)
library(dplyr)
source("fonction_codif_genotype.R")


# Lire les fichiers VCF

aa_vcf_control = read.vcfR("bcf_snp_AA.vcf")
as_vcf_case = read.vcfR("bcf_snp_AS.vcf")
ss_vcf_case = read.vcfR("bcf_snp_SS.vcf")

# Fonction pour créer une table de contingence pour un SNP
create_contingency_table = function(geno_case, geno_control) {
  ref_case = sum(geno_case == "0",na.rm = T) * 2 + sum(geno_case == "1",na.rm = T)
  alt_case = sum(geno_case == "2",na.rm = T) * 2 + sum(geno_case == "1",na.rm = T)
  
  ref_control = sum(geno_control == "0",na.rm = T) * 2 + sum(geno_control == "1",na.rm = T)
  alt_control = sum(geno_control == "2",na.rm = T) * 2 + sum(geno_control == "1",na.rm = T)
  
  table = matrix(c(ref_case, alt_case, ref_control, alt_control), 
                  nrow = 2, byrow = TRUE)
  colnames(table) = c("Cas", "Témoins")
  rownames(table) = c("Ref", "Alt")
  
  return(table)
}


test_association = function(vcf_case, vcf_control) {
  # Extraire les génotypes des individus
  geno_case_all = extract.gt(vcf_case, element = "GT", as.numeric = FALSE)
  geno_control_all = extract.gt(vcf_control, element = "GT", as.numeric = FALSE)
  
  commun = rownames(geno_control_all) %in% rownames(geno_case_all)
  snp_commun = row.names(geno_control_all[commun, ])
  
  geno_case = geno_case_all[snp_commun, ] %>% t() %>% 
    #as.matrix()
    as.data.frame() %>% codage_geno_012()
    
  geno_control = geno_control_all[snp_commun, ] %>% t() %>% 
    #as.matrix()
    as.data.frame() %>% codage_geno_012()
  
  # Appliquer le test du Chi-2 pour chaque SNP
  
  results = data.frame(
    SNP = NA,
    statistics = NA,
    p_value = NA,
    test = NA
  )
  for (snp in colnames(geno_case)) {
    if (snp %in% colnames(geno_control)) {
      table = create_contingency_table(geno_case[, snp], geno_control[, snp])
      if (min(table, na.rm = T) < 5) {
        fisher = fisher.test(table)
        results = rbind(results,
                        c(snp, fisher$estimate, fisher$p.value, "fisher.test"))
      } else{
        chi2_test = chisq.test(table)
        results = rbind(results,
                        c(snp, chi2_test$statistic, chi2_test$p.value, "Chi2"))
      }
    }
  }
  return(results)
}

AS_AA = test_association(as_vcf_case,aa_vcf_control)
SS_AA = test_association(ss_vcf_case,aa_vcf_control)
SS_AS = test_association(ss_vcf_case,as_vcf_case)


AS_AA = AS_AA[-1,]
SS_AA = SS_AA[-1,]
SS_AS = SS_AS[-1,]
# write.csv2(AS_AA,file = "as_aa.csv")
# write.csv2(SS_AA,file = "ss_aa.csv")

as_bind = cbind(stringr::str_split(AS_AA$SNP,"_v3_",simplify = T),AS_AA)
ss_bind = cbind(stringr::str_split(SS_AA$SNP,"_v3_",simplify = T),SS_AA)
ss_as_bind = cbind(stringr::str_split(SS_AS$SNP,"_v3_",simplify = T),SS_AS)


names(as_bind)[1:2] = c("Chrom","Pos")
names(ss_bind)[1:2] = c("Chrom","Pos")
names(ss_as_bind)[1:2] = c("Chrom","Pos")

ss_bind = ss_bind %>% mutate(across(c(Pos,p_value), as.numeric)) %>% 
  filter(Chrom!="Pf3D7_MIT")
as_bind = as_bind %>% mutate(across(c(Pos,p_value), as.numeric)) %>% 
  filter(Chrom!="Pf3D7_MIT")
ss_as_bind = ss_as_bind %>% mutate(across(c(Pos,p_value), as.numeric)) %>% 
  filter(Chrom!="Pf3D7_MIT")

xlsx::write.xlsx(as_bind,"chi2_as_aa.xlsx", sheetName = "chi2_aa_vs_as")
xlsx::write.xlsx(ss_bind,"chi2_ss_aa.xlsx", sheetName = "chi2_aa_vs_ss")
xlsx::write.xlsx(ss_as_bind,"chi2_ss_as_bind.xlsx", sheetName = "chi2_ss_vs_as")

library(ggplot2)
gplot_chi2_aa_as = ggplot(as_bind,aes(x=Pos,y = p_value)) +
  geom_line() +
  facet_wrap(~Chrom, scales = "free", ncol = 4) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")

ggsave("gplot_chi2_aa_as.png", plot = gplot_chi2_aa_as, 
       width = 13, height = 6)
gplot_chi2_aa_as

gplot_chi2_aa_ss = ggplot(ss_bind,aes(x=Pos,y = p_value)) +
  geom_line() +
  facet_wrap(~Chrom, scales = "free", ncol = 4) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")

ggsave("gplot_chi2_aa_ss.png", plot = gplot_chi2_aa_ss, 
       width = 13, height = 6)
gplot_chi2_aa_ss


gplot_chi2_ss_as = ggplot(ss_as_bind,aes(x=Pos,y = p_value)) +
  geom_line() +
  facet_wrap(~Chrom, scales = "free", ncol = 4) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")

ggsave("gplot_chi2_ss_as.png", plot = gplot_chi2_ss_as, 
       width = 13, height = 6)
gplot_chi2_ss_as


as_bind = as_bind %>% arrange(p_value)
ss_bind = ss_bind %>% arrange(p_value)
ss_as_bind = ss_as_bind %>% arrange(p_value)

xlsx::write.xlsx(as_bind,"chi2_as_aa_triee.xlsx", sheetName = "chi2_aa_vs_as")
xlsx::write.xlsx(ss_bind,"chi2_ss_aa_triee.xlsx", sheetName = "chi2_aa_vs_ss")
xlsx::write.xlsx(ss_as_bind,"chi2_ss_as_triee.xlsx", sheetName = "chi2_ss_vs_as")

