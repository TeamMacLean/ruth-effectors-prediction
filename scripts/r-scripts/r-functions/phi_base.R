library(tidyverse)

phi_base <- data.table::fread("data/phi-base_current.csv", header = TRUE)

# phi_small <- phi_base %>%
#   dplyr::filter(stringr::str_detect(V31, "Effector"))

phi_small <- phi_base %>%
  dplyr::filter(stringr::str_detect(V35, "(plant avirulence determinant)"))


phi_small <- phi_base %>%
  dplyr::filter_all(any_vars(str_detect(., 'plant avirulence determinant')))

phi_small %>%
  dplyr::filter(!str_detect(`Mutant Phenotype`, 'plant avirulence determinant')) %>%
  dplyr::select(`Gene`)

phi_proteinID <- phi_small %>%
  dplyr::select(`Protein ID`)


# find th eunique value

proteinID_unique <- unique(phi_proteinID) %>%
  dplyr::filter_all(any_vars(!str_detect(., 'no data found')))


