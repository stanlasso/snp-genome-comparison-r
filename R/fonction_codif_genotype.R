
library(dplyr)

# fonction codification du genotype 0_1_1
# "0/0" ~ 0, 
# "1/0", "0/1" ~ 1,
# "1/1" ~ 1
codage_geno_011 = function(gt) {
  gt = gt %>%
    mutate(across(everything(), ~ case_when(
      . == "0/0" ~ 0, 
      . %in% c("1/0", "0/1") ~ 1,
      . == "1/1" ~ 1
    )))
  
  return(gt)
}



# fonction codification du genotype 0_1_2
# "0/0" ~ 0, 
# "1/0", "0/1" ~ 1,
# "1/1" ~ 2
codage_geno_012 = function(gt) {
  gt = gt %>%
    mutate(across(everything(), ~ case_when(
      . == "0/0" ~ 0, 
      . %in% c("1/0", "0/1") ~ 1,
      . == "1/1" ~ 2
    )))
  
  return(gt)
}


# fonction codification du genotype 1_2_3
# "0/0" ~ 1, 
# "1/0", "0/1" ~ 2,
# "1/1" ~ 3
codage_geno_123 = function(gt) {
  gt = gt %>%
    mutate(across(everything(), ~ case_when(
      . == "0/0" ~ 1, 
      . %in% c("1/0", "0/1") ~ 2,
      . == "1/1" ~ 3
    )))
  
  return(gt)
}
