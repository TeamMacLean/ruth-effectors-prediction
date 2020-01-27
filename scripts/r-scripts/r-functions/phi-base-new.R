library(tidyverse)

phi_base <- data.table::fread("~/Downloads/phi-base_current.csv", header = TRUE)

phi_small <- phi_base %>%
  dplyr::filter_all(any_vars(str_detect(., 'plant avirulence determinant')))

# write.csv(phi_small, "~/Downloads/phi-base-avirulenc_current.csv")


phi_small %>%
  dplyr::filter(!str_detect(`Mutant Phenotype`, 'plant avirulence determinant')) %>%
  dplyr::select(`Protein ID`)
