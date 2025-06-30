## ----setup, include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(echo = T, comment = "", warning = F)
library(vcfR) # Pour lire et manipuler les VCF en R
library(SNPRelate)  # GDS, PCA, FST, IBS, etc
library(readxl)
library(stringr)
library(mmod)
library(adegenet)
library(pegas)
library(ggplot2)
library(ggrepel)
library(vegan)   # Pour PERMANOVA (adonis, adonis2)
library(dplyr)   # Manipulations data.frame
load("Jost.RData")

#


## ----include=FALSE-------------------------------------------------------------------
# <!---BLOCK_LANDSCAPE_START--->


## ------------------------------------------------------------------------------------
# # Installation des packages utiles
# install.packages("vcfR")        # Pour lire et manipuler les VCF en R
# install.packages("gdsfmt")
# install.packages("vegan")       # Pour PERMANOVA (adonis, adonis2)
# install.packages("dplyr")       # Manipulations data.frame
# BiocManager::install("SNPRelate") # GDS, PCA, FST, IBS, etc


## ------------------------------------------------------------------------------------
vcf_file = "bcf_snp_merge30filter55vcf_moins_AC.vcf"    # ou "mon_fichier.vcf.gz"
vcf_data = read.vcfR(vcf_file,verbose = F)  # charge le VCF en mémoire

# vcf_data est un objet de classe vcfR contenant :
#  - vcf_data@meta : métadonnées du VCF
#  - vcf_data@fix  : colonnes fixes (CHROM, POS, REF, ALT, etc.)
#  - vcf_data@gt   : génotypes pour chaque échantillon



## ------------------------------------------------------------------------------------
library(readxl)
populations = read_excel("pop.data_Allele_P2.xlsx")

# str_replace_all(pop$ID,paste0("(", "S", ").*"), "\\1")

pop = str_replace_all(colnames(vcf_data@gt)[-1],paste0("S", ".*"), "")
pop = data.frame(ID=pop,ID2=colnames(vcf_data@gt)[-1]) %>% select(-c(ID))
pop.data = merge(pop,populations,by = "ID2",sort = F)
#write.csv2(pop.data,file="pop.data.csv")
sample_info = pop.data
names(sample_info)[1] = "sample.id"
sample_info$GENOTYPE = as.factor(sample_info$GENOTYPE)
all(colnames(vcf_data@gt)[-1] == pop.data$ID2)


## ------------------------------------------------------------------------------------
gds_file = "bcf_snp_merge30filter55vcf_moins_AC.gds"
snpgdsVCF2GDS(vcf_file, gds_file, method = "copy.num.of.ref")
#  method = "copy.num.of.ref" => code le génotype selon le nombre d’allèles de référence


## ------------------------------------------------------------------------------------
genofile = snpgdsOpen(gds_file)


## ------------------------------------------------------------------------------------
snpset = snpgdsSelectSNP(genofile,
                          #maf = 0.01,           # MAF>1%
                          #missing.rate = 0.05,  # <5% de missing
                          autosome.only=FALSE)

# snpset = snpgdsSelectSNP(genofile,
#                           autosome.only=FALSE)


# on peut créer un GDS filtré :
filtered_gds_file = "bcf_snp_merge30filter55vcf_moins_AC_filtered.gds"
name_genofile = "bcf_snp_merge30filter55vcf_moins_AC.gds"
snpgdsCreateGenoSet(src.fn = name_genofile, dest.fn = filtered_gds_file, snp.id = snpset)

filtered_gds = snpgdsOpen(filtered_gds_file)


## ------------------------------------------------------------------------------------
pca_res = snpgdsPCA(filtered_gds, autosome.only=FALSE) 
# autosome.only=TRUE si on veut ignorer X/Y


## ------------------------------------------------------------------------------------
pc_percent = pca_res$varprop * 100 
pc_percent %>% format(nsmall =  2,scientific = F)

pc = paste("PC",1:length(pc_percent),sep = "")
pc_percent_df = data.frame(
  PC = factor(pc,levels = pc),
  Percent = round(pc_percent,1),
  Couleur = ifelse(pc_percent>3.5, "Important","Moins important"))


## ------------------------------------------------------------------------------------
# Créer le graphe avec ggplot2
Screeplot = ggplot(pc_percent_df[1:23,], aes(x = PC, y = Percent, fill = Couleur)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Important" = "deeppink", "Moins important" = "blue")) +
  geom_text(aes(label = paste0(Percent, "%")), vjust = -0.5,hjust = 0.4 ,
            fontface = "bold",size = 3) +
  labs(
    title = "Variance explained by each component",
    y = "Variance explained (%)",
    x = ""
  ) +
  theme_classic() +
   theme(
    text = element_text(face = "bold", size = 12, family = "Times New Roman"),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", 
                               color = "black"),
    axis.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  ylim(0, max(pc_percent_df$Percent) + 1)
  # theme(
  #   axis.text.x = element_text(angle = 45, hjust = 1),
  #   legend.position = "none",
  #   plot.title = element_text(face = "bold", hjust = 0.5)
  # ) +
  

ggsave(filename = "Screeplot.png",dpi = 300, plot = Screeplot, 
       width = 8, height = 7,bg = "white")



## ------------------------------------------------------------------------------------
pca_tab = data.frame(
  sample.id = pca_res$sample.id,
  PC1 = pca_res$eigenvect[, 1],
  PC2 = pca_res$eigenvect[, 2],
  PC3 = pca_res$eigenvect[, 3],
  PC4 = pca_res$eigenvect[, 4],
  PC5 = pca_res$eigenvect[, 5],
  PC6 = pca_res$eigenvect[, 6]
)
head(pca_tab)


## ------------------------------------------------------------------------------------
pca_annot = pca_tab %>% 
  left_join(sample_info, by="sample.id")  # jointure sur sample.id

# Plot
# plot(pca_annot$PC1, pca_annot$PC2,
#      col = as.factor(pca_annot$GENOTYPE), pch=19,
#      xlab = paste0("PC1 (", round(pc_percent[1],2),"% )"),
#      ylab = paste0("PC2 (", round(pc_percent[2],2),"% )"),
#      main = "PCA - Couleur selon groupe Hb")
# legend("topright", legend=unique(pca_annot$GENOTYPE),
#        col=1:length(unique(pca_annot$GENOTYPE)), pch=19)



## ------------------------------------------------------------------------------------
# Charger les packages
#install.packages("plotly")
library(plotly)

# Visualisation 3D avec plotly
plot_ly(
  pca_annot,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~GENOTYPE,
  colors = c("deeppink", "blue", "forestgreen"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "Graphique 3D des trois premières composantes",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )



## ----eval=FALSE----------------------------------------------------------------------
# plot3D::scatter3D(
#   x = pca_annot$PC1,
#   y = pca_annot$PC2,
#   z = pca_annot$PC3,
#   #angle = angle,
#   #pch = shapes[as.numeric(pca_annot$GENOTYPE)],
#   col= pca_annot$GENOTYPE
# )



## ------------------------------------------------------------------------------------
gplo_pca12_avlab = ggplot(data = pca_annot) +
  geom_point(aes(x = PC1, y = PC2, colour = GENOTYPE)) +
  labs(x = paste("PC1", round(pc_percent[1],2),"%"),
       y = paste("PC2", round(pc_percent[2],2),"%"),
       #title = "PCA - Couleur selon groupe Hb"
       ) +
  geom_text_repel(aes(x = PC1, y = PC2, label = sample.id), vjust = -1,
                   max.overlaps = 30) +
   #geom_label(aes(x = PC1, y = PC2,label = sample.id), vjust = -1, alpha = 0.5) +
  theme_classic()+
   theme(
    # Texte général en gras
    text = element_text(face = "bold", size = 12, family = "Times New Roman"), # Titre en gras
    plot.title = element_text(face = "bold", hjust = 0.5),# Axes en gras
    axis.title = element_text(face = "bold"), # Légende en gras
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )
  # theme(text = element_text(family = "Times New Roman"))

ggsave(filename = "gplo_pca12_avlab.png",dpi = 300, plot = gplo_pca12_avlab)


## ------------------------------------------------------------------------------------
gplo_pca12 = ggplot(data = pca_annot) +
  geom_point(aes(x = PC1, y = PC2, colour = GENOTYPE)) +
  labs(x = paste("PC1", round(pc_percent[1],2),"%"),
       y = paste("PC2", round(pc_percent[2],2),"%"),
       #title = "PCA - Couleur selon groupe Hb"
       ) +
  theme_classic()+
   theme(
    # Texte général en gras
    text = element_text(face = "bold", size = 12, family = "Times New Roman"),
    # Titre en gras (peut être ajusté séparément)
    plot.title = element_text(face = "bold", hjust = 0.5),
    # Axes en gras
    axis.title = element_text(face = "bold"),
    # Légende en gras
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )
  # theme(text = element_text(family = "Times New Roman"))

ggsave(filename = "gplo_pca12.png",dpi = 300, plot = gplo_pca12)


## ----eval=FALSE----------------------------------------------------------------------
# # library(ggrepel)
# # ggplot(data = pca_annot) +
# #   geom_point(aes(x = PC1, y = PC2, colour = GENOTYPE)) +
# #   labs(
# #     x = paste("PC1", round(pc_percent[1], 2), "%"),
# #     y = paste("PC2", round(pc_percent[2], 2), "%"),
# #     title = "PCA - Couleur selon groupe Hb"
# #   ) +
# #   geom_label_repel(aes(x = PC1, y = PC2, label = sample.id), size = 2.5,
# #                    segment.square	= T) +
# #   theme_classic()+
# # theme(text = element_text(family = "Times New Roman"))
# 


## ------------------------------------------------------------------------------------
gplo_pca13 = ggplot(data = pca_annot) +
  geom_point(aes(x = PC1, y = PC3, colour = GENOTYPE)) +
  labs(x = paste("PC1", round(pc_percent[1],2),"%"),
       y = paste("PC3", round(pc_percent[3],2),"%"),
       # title = "PCA - Couleur selon groupe Hb"
       ) +
  theme_classic()+
   theme(
    # Texte général en gras
    text = element_text(face = "bold", size = 12, family = "Times New Roman"),
    # Titre en gras (peut être ajusté séparément)
    plot.title = element_text(face = "bold", hjust = 0.5),
    # Axes en gras
    axis.title = element_text(face = "bold"),
    # Légende en gras
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )
  # theme(text = element_text(family = "Times New Roman"))

ggsave(filename = "gplo_pca13.png",dpi = 300, plot = gplo_pca13)


## ------------------------------------------------------------------------------------
gplo_pca23 = ggplot(data = pca_annot) +
  geom_point(aes(x = PC2, y = PC3, colour = GENOTYPE)) +
  labs(x = paste("PC2", round(pc_percent[2],2),"%"),
       y = paste("PC3", round(pc_percent[3],2),"%"),
       title = "PCA - Couleur selon groupe Hb") +
  theme_classic()+
  theme(text = element_text(family = "Times New Roman"))

ggsave(filename = "gplo_pca23.png",dpi = 300, plot = gplo_pca23)


## ------------------------------------------------------------------------------------
# sample_info doit avoir un ordre aligné à pca_res$sample.id ou filtered_gds$sample.id
# Dans l'idéal, on a sample_info$sample.id, sample_info$pop
fst_res = snpgdsFst(
  gdsobj      = filtered_gds,
  sample.id   = sample_info$sample.id,
  population  = sample_info$GENOTYPE,        # vecteur de pop
  method      = "W&C84",                 # ou "W&H02"
  autosome.only=FALSE
)

fst_res$Fst       # Valeur globale
fst_res$MeanFst   # FST moyen sur l'ensemble des SNP


## ------------------------------------------------------------------------------------

ibs_res = snpgdsIBS(filtered_gds, autosome.only=FALSE)

# ibs_res$ibs est une matrice de similarité (n x n)
# On convertit en distance
dist_mat = 1 - ibs_res$ibs

rownames(dist_mat) = ibs_res$sample.id
colnames(dist_mat) = ibs_res$sample.id


## ------------------------------------------------------------------------------------
dist_mat_imp = dist_mat
dist_mat_imp[is.nan(dist_mat)] = 0

#dist_obj = as.dist(dist_mat)
dist_obj = as.dist(dist_mat_imp)

#dist_obj = as.dist(dist_mat)


## ------------------------------------------------------------------------------------

library(vegan)

# On doit aligner le vecteur group (HbAA, HbAS, HbSS) à l'ordre rownames(dist_mat)
metadata = data.frame(
  sample.id = ibs_res$sample.id,
  group = sample_info$GENOTYPE # ou sample_info$group si c'est la même correspondance
)

# reorder metadata si besoin, aligner rownames etc.
rownames(metadata) = metadata$sample.id

# PERMANOVA via adonis2 (ou adonis)
adonis_res = adonis2(dist_obj ~ group, data=metadata)
adonis_res


## ----include=FALSE-------------------------------------------------------------------
#<!---BLOCK_LANDSCAPE_STOP--->


## ------------------------------------------------------------------------------------
# install.packages("mmod")
# install.packages("adegenet")  # Pour gérer les données génétiques
# install.packages("pegas")     # Optionnel pour manipuler les données SNP
# library(mmod)
# library(adegenet)
# library(pegas)


## ------------------------------------------------------------------------------------
genind_data = vcfR2genind(vcfR::read.vcfR(vcf_file, verbose = F))
pop(genind_data) = pop.data$GENOTYPE


## ----eval=FALSE----------------------------------------------------------------------
# # Calculer Jost's D
# jost_d = D_Jost(genind_data)
# 


## ----eval=FALSE----------------------------------------------------------------------
# print(jost_d)


## --------------------------------------------------------------------------------
# Matrice pairwise de Jost's D
pairwise_jost_d = pairwise_D(genind_data)

# Afficher la matrice
print(pairwise_jost_d)



## ------------------------------------------------------------------------------------
print(pairwise_jost_d)


## ------------------------------------------------------------------------------------
heatmap(as.matrix(pairwise_jost_d), 
        main = "Pairwise Jost's D Heatmap", 
        xlab = "Population", 
        ylab = "Population", 
        col = heat.colors(10))


## ------------------------------------------------------------------------------------
dist_matrix = as.dist(pairwise_jost_d)
hc = hclust(dist_matrix, method = "average")
plot(hc, main = "Dendrogramme basé sur Jost's D", xlab = "", sub = "")

