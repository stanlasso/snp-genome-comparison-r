## ------------------------------------------------------------------------------------
#install.packages(c("vcfR", "tidyverse", "ggplot2"))
library(vcfR)
library(tidyverse)
library(ggplot2)


## ------------------------------------------------------------------------------------
# Définition des fichiers par groupe génotypique
vcf_directory = "/home/sherif/Bureau/yacouba_p2/Mapped/variant_calling_mapped/vcffilter55/Gene/fusion/bcftools_snp_ok/echantillon"
#vcf_files = list.files(vcf_directory, pattern = "*.vcf$", full.names = TRUE)

vcf_files = paste(vcf_directory,"/",c("bcf_snp_Ank086S142_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_Ank087S143_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK96S76_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK29S107_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ank123_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK67S102_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK33S100_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK28S106_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK55S101_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK54S108_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC45CredepositS190_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC41AredepositS193_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC41BredepositS194_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC45AredepositS188_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC45BredepositS189_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK113S77_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP11S103_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_ANK125S5_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP61S8_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP06S99_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP57S78_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP63S9_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP07S98_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP02S98_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP53S6_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP32S99_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_YOP55S7_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC31BredepositS192_sort_mapped_dup_FILTERED55.vcf", "bcf_snp_TC31AredepositS191_sort_mapped_dup_FILTERED55.vcf"), sep = "")
# Associer chaque fichier à son groupe (AA, AS, SS)
geno_types = rep(c("AA", "AS", "SS"),c(10,9,10))  
sample_names = basename(vcf_files)

# Création d'un tableau associant fichiers et génotypes
metadata = data.frame(Sample = sample_names, File = vcf_files, Genotype = geno_types)


## ------------------------------------------------------------------------------------
# Initialiser un data.frame pour stocker les résultats
results = data.frame()

# Boucle sur chaque fichier VCF
for (i in 1:nrow(metadata)) {
  
  # Lire le fichier VCF
  vcf = read.vcfR(metadata$File[i], verbose = FALSE)
  
    vcf_data = vcf@fix %>% as.data.frame()
  # Extraire DP et QUAL
  dp_values = extract.gt(vcf, element = "DP", as.numeric = TRUE)
  qual_values = as.numeric(vcf_data$QUAL)
  
  # Extraction du nombre total de SNPs
  num_snps = nrow(vcf_data)
  
  # Détection des types de SNPs
  ref_allele = vcf_data$REF
  alt_allele = vcf_data$ALT
  transitions = sum((ref_allele == "A" & alt_allele == "G") | (ref_allele == "G" & alt_allele == "A") |
                     (ref_allele == "C" & alt_allele == "T") | (ref_allele == "T" & alt_allele == "C"))
  transversions = num_snps - transitions
  
  # Ajouter au tableau des résultats
  results = rbind(results, data.frame(
    Sample = metadata$Sample[i],
    Genotype = metadata$Genotype[i],
    Mean_DP = mean(dp_values, na.rm = TRUE),
    Median_DP = median(dp_values, na.rm = TRUE),
    SD_DP = sd(dp_values, na.rm = TRUE),
    Mean_QUAL = mean(qual_values, na.rm = TRUE),
    Median_QUAL = median(qual_values, na.rm = TRUE),
    SD_QUAL = sd(qual_values, na.rm = TRUE),
    Num_SNPs = num_snps,
    Transitions = transitions,
    Transversions = transversions
  ))
}

# Sauvegarde des résultats
write.csv(results, "resultats_vcf_par_genotype.csv", row.names = FALSE)


## ------------------------------------------------------------------------------------
# Calcul des moyennes par génotype
summary_results = results %>%
  group_by(Genotype) %>%
  summarise(
    Moy_DP = mean(Mean_DP, na.rm = TRUE),
    SD_DP = sd(Mean_DP, na.rm = TRUE),
    Moy_QUAL = mean(Mean_QUAL, na.rm = TRUE),
    SD_QUAL = sd(Mean_QUAL, na.rm = TRUE),
    Total_SNPs = sum(Num_SNPs),
    Moy_Transitions = mean(Transitions, na.rm = TRUE),
    Moy_Transversions = mean(Transversions, na.rm = TRUE)
  )

# Affichage des résultats

write.csv(results, "resultats_vcf_stats_genotype.csv", row.names = FALSE)



## ------------------------------------------------------------------------------------
gplot_dp_gen = ggplot(results, aes(x = Genotype, y = Mean_DP, fill = Genotype)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution de la Couverture Moyenne (DP) par Génotype", x = "Génotype", y = "Profondeur Moyenne (DP)") + theme_classic()

ggsave("gplot_dp_gen.png",gplot_dp_gen)


## ------------------------------------------------------------------------------------
gplot_snp_gen = ggplot(results, aes(x = Genotype, y = Num_SNPs, fill = Genotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Nombre total de SNPs par Génotype", x = "Génotype", y = "Nombre de SNPs") + theme_classic()

ggsave("gplot_snp_gen.png",gplot_snp_gen)


## ------------------------------------------------------------------------------------
gplot_gen_trsit_trvers = ggplot(results, aes(x = Genotype, y = Transitions / Transversions, fill = Genotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Ratio Transitions / Transversions par Génotype", x = "Génotype", y = "Ratio")  + theme_classic()

ggsave("gplot_gen_trsit_trvers.png",gplot_gen_trsit_trvers)


## ------------------------------------------------------------------------------------
gplot_qual_gen = ggplot(results, aes(x = Genotype, y = Mean_QUAL, fill = Genotype)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution de la qualité Moyenne (QUAL) par Génotype", x = "Génotype", y = "Profondeur Moyenne (DP)") + theme_classic()

ggsave("gplot_qual_gen.png",gplot_qual_gen)


## ------------------------------------------------------------------------------------
#install.packages("VennDiagram")
library(VennDiagram)


## ------------------------------------------------------------------------------------
# Lecture des fichiers VCF et extraction des SNPs uniques pour chaque groupe
extract_snps = function(vcf_file) {
  vcf = read.vcfR(vcf_file, verbose = FALSE)
  vcf_data = vcf@fix %>% as.data.frame()
  return(unique(vcf_data$POS))  # Extraction des positions des SNPs
}

# Séparation des fichiers par génotype
aa_files = metadata$File[metadata$Genotype == "AA"]
as_files = metadata$File[metadata$Genotype == "AS"]
ss_files = metadata$File[metadata$Genotype == "SS"]

# Extraction des SNPs pour chaque groupe
snps_aa = unique(unlist(lapply(aa_files, extract_snps)))
snps_as = unique(unlist(lapply(as_files, extract_snps)))
snps_ss = unique(unlist(lapply(ss_files, extract_snps)))

# Création du diagramme de Venn
venn.plot = draw.triple.venn(
  area1 = length(snps_aa),
  area2 = length(snps_as),
  area3 = length(snps_ss),
  n12 = length(intersect(snps_aa, snps_as)),   # SNPs communs entre AA et AS
  n13 = length(intersect(snps_aa, snps_ss)),   # SNPs communs entre AA et SS
  n23 = length(intersect(snps_as, snps_ss)),   # SNPs communs entre AS et SS
  n123 = length(Reduce(intersect, list(snps_aa, snps_as, snps_ss))),  # SNPs communs aux trois groupes
  category = c("AA", "AS", "SS"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2
)

# Sauvegarde de l’image
pdf("venn_snp_distribution.pdf")
grid.draw(venn.plot)
dev.off()

