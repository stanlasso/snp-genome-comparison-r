library(stringr)

list_file = system("ls *out62.txt", intern = T)
# noms = str_replace_all(list_file,paste0("S", ".*"), "") %>% 
#   str_replace_all("bcf_snp_","") %>% 
#   str_replace_all(paste0("_dup", ".*"), "")

colonnes = c("#Uploaded_variation", "Location", "Allele", "Gene", "Feature",
             "Feature_type", "Consequence", "cDNA_position", "CDS_position",
             "Protein_position", "Amino_acids", "Codons", "Existing_variation", 
             "IMPACT", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL",
             "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "CANONICAL", "MANE_SELECT",
             "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT",
             "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "SOURCE", "GENE_PHENO",
             "SIFT", "PolyPhen", "EXON", "INTRON", "DOMAINS", "miRNA", "HGVSc",
             "HGVSp", "HGVS_OFFSET", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF",
             "SAS_AF", "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF",
             "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF", 
             "gnomADe_NFE_AF", "gnomADe_OTH_AF", "gnomADe_SAS_AF", 
             "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF",
             "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", 
             "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_OTH_AF", 
             "gnomADg_SAS_AF", "MAX_AF", "MAX_AF_POPS", "CLIN_SIG", "SOMATIC",
             "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", 
             "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS", 
             "PlasmoDB-62_Pfalciparum3D7.gff.gz")

noms = list_file %>% str_replace_all("bcf_snp_","") %>% 
  #str_replace_all(paste0("S", ".*"), "") %>% 
  str_replace_all("_sort_mapped_dup_FILTERED55_out62.txt", "")
for (i in 1:length(list_file)) {
  ff=read.table(list_file[i])
  names(ff) = colonnes
  assign(noms[i],ff)
}


noms_AA = c("Ank086S142", "Ank087S143", "ANK96S76", "ANK29S107", "ank123", 
            "ANK67S102", "ANK33S100", "ANK28S106", "ANK55S101", "ANK54S108")
noms_AS = c("TC45CredepositS190", "TC41AredepositS193", "TC41BredepositS194", 
            "TC45AredepositS188", "TC45BredepositS189", "ANK113S77", "YOP11S103", 
            "ANK125S5", "YOP61S8")
noms_SS = c("YOP06S99", "YOP57S78", "YOP63S9", "YOP07S98", "YOP02S98", 
            "YOP53S6", "YOP32S99", "YOP55S7", "TC31BredepositS192", 
            "TC31AredepositS191")


library(dplyr)

# Liste des bases

global_genotype_function = function(sample_genotype, nom_col){
  bases = list()
  for (ech in sample_genotype) {
    dat = get(ech)
    bases[[ech]] = dat %>%  select(-c(STRAND))

  }
  # Ajouter l’indicateur à chaque base
  bases = lapply(names(bases), function(nom) {
    df = bases[[nom]]
    df[[nom]] = 1
    df
  })

  # Fusion successive sur toutes les colonnes communes (ID, Nom, Age)
  merged = Reduce(function(x, y) full_join(x, y, by = nom_col), bases)
  
  # Remplacer les NA dans les colonnes d’indicateur par 0
  merged[sample_genotype] = lapply(merged[sample_genotype],
                                                 function(x) ifelse(is.na(x), 0, x))
  
  merged = merged %>% 
    select(!(HGVS_OFFSET:`PlasmoDB-62_Pfalciparum3D7.gff.gz`)) %>% 
    select(!(CANONICAL:miRNA)) %>% 
    select(!c(SYMBOL_SOURCE, HGNC_ID, FLAGS))
  
  
  # Calcul de la prévalence
  merged$prevalence = rowSums(merged[sample_genotype]) / length(sample_genotype)
  #merged$prevalence = round(merged$prevalence,3)
  return(merged)
}


AA = global_genotype_function(sample_genotype = noms_AA,nom_col = colonnes[-16])
AS = global_genotype_function(sample_genotype = noms_AS,nom_col = colonnes[-16])
SS = global_genotype_function(sample_genotype = noms_SS,nom_col = colonnes[-16])

write.table(AA, file = "data/AA_echantillon_mutation_prevalence.txt",sep = "\t",
            na = "")
write.table(AS, file = "data/AS_echantillon_mutation_prevalence.txt",sep = "\t",
            na = "")
write.table(SS, file = "data/SS_echantillon_mutation_prevalence.txt",sep = "\t",
            na = "")
