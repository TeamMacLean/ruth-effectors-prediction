library(tidyverse)

# read the mapping result from phi-base unique proteinIDs to uniprot -- retrieved in .csv
uniprot_raw <- data.table::fread("data/uniprot-results-mapped-raw.csv", fill = TRUE, sep = "\t")

# rename a column of protein IDs obtained from phi-base, and the entry IDs from
uniprot_raw <- uniprot_raw %>%
  rename(phi_ids = `yourlist:M201904266746803381A1F0E0DB47453E0216320D065C43S`, uniprot_ids = `Entry`)

# getting the ids which are different between phi-base and uniprot
diff_phi_uniprot_ids <- uniprot_raw %>%
  dplyr::filter(!(phi_ids %in% intersect(uniprot_raw[['phi_ids']], uniprot_raw[['uniprot_ids']]))) %>%
  select(`phi_ids`, `uniprot_ids`, `Organism`)

# matching the ids that can be mapped vs ids that have sequence available

uniprot_ids_with_unavailable_seq <- uniprot_raw %>%
  dplyr::filter(!(uniprot_ids %in% intersect(uniprot_effector_short_pathogen_name[['protein_id']], uniprot_raw[['uniprot_ids']])))

uniprot_ids_with_available_seq <- uniprot_raw %>%
  dplyr::filter((uniprot_ids %in% intersect(uniprot_effector_short_pathogen_name[['protein_id']], uniprot_raw[['uniprot_ids']])))

# matching the IDs

uniprot_effector_with_proteinID_correction <- uniprot_effector_short_pathogen_name %>%
  mutate('protein_id' = ifelse(protein_id == uniprot_ids_with_available_seq[['phi_ids']], protein_id, uniprot_ids_with_available_seq[['phi_ids']])) %>%
  select(`protein_id`, `pathogen_short`, `sequence`)

# since only 482 ids available then now we can filter

phi_proteinID_list <- phi_effector_proteinID_unique %>%
  dplyr::select(`Protein ID`) %>%
  unlist()

uniprot_proteinID_list <- uniprot_effector_with_proteinID_correction %>%
  dplyr::select(protein_id) %>%
  unique() %>%
  unlist()

phi_ids_available_in_uniprot <- phi_effector_proteinID_unique %>%
  dplyr::filter(`Protein ID` %in% dplyr::intersect(phi_proteinID_list, uniprot_proteinID_list))

#  getting the difference of the pathogen name and the ids

uniprot_pathogen_name_unique <- uniprot_effector_with_proteinID_correction %>%
  dplyr::select(pathogen_short) %>%
  unique() %>%
  unlist()

phi_pathogen_name_unique <- phi_effector_proteinID_unique %>%
  dplyr::select(`Pathogen species`) %>%
  unique() %>%
  unlist()

phi_uniprot_diff_pathogen_name <- uniprot_effector_with_proteinID_correction %>%
  dplyr::filter(!(pathogen_short %in% intersect(uniprot_pathogen_name_unique, phi_pathogen_name_unique)))


phi_uniprot_diff_pathogen_name  <- phi_uniprot_diff_pathogen_name %>%
  dplyr::select(`protein_id`, `pathogen_short`) %>%
  mutate(pathogen_short = ifelse(pathogen_short == "Salmonella typhi", "Salmonella typhimurium", pathogen_short)) %>%
  group_by(pathogen_short) %>%
  summarise(count = n(), protein_id = paste(protein_id, collapse = " "))


# getting all of the ids in

ids_diff_pathogen_name <- phi_uniprot_diff_pathogen_name  %>%
  dplyr::summarise(all_ids = paste(protein_id, collapse = " ")) %>%
  unlist() %>%
  strsplit(" ") %>%
  unlist() %>%
  unname()


# getting the true name for then phi-base
phi_correct_pathogen_name <- phi_ids_available_in_uniprot %>%
  dplyr::filter(`Protein ID` %in% ids_diff_pathogen_name) %>%
  dplyr::select(`Pathogen species`, `Protein ID`) %>%
  rename(pathogen_short = `Pathogen species`, protein_id = `Protein ID`) %>%
  group_by(protein_id, pathogen_short) %>%
  summarise()

effector_data <- uniprot_effector_with_proteinID_correction %>%
  left_join(., phi_correct_pathogen_name, by = "protein_id") %>%
  mutate(
    pathogen_short = ifelse(is.na(pathogen_short.y), pathogen_short.x, pathogen_short.y),
    name_src = ifelse(is.na(pathogen_short.y), "uniprot", "phi_base"),
    name_src = as.factor(name_src)
  ) %>%
  dplyr::select(-pathogen_short.x, -pathogen_short.y)

# getting the sample of non-effecor protein

#  count how many sample data of effector per
effector_data_count_each_organism_sample <- effector_data %>%
  group_by(pathogen_short) %>%
  summarise(count = n())

effector_data_count_each_organism_sample_list <- effector_data_count_each_organism_sample %>%
  dplyr::summarise(all_ids = paste(pathogen_short, collapse = "|")) %>%
  unlist() %>%
  strsplit("\\|") %>%
  unlist() %>%
  unname()

effector_rowids <- phi_base %>%
  tibble::rowid_to_column() %>%
  dplyr::filter_all(any_vars(str_detect(., "plant avirulence determinant"))) %>%
  select(rowid) %>%
  unlist() %>%
  unname()

noneffector <- phi_base %>%
  tibble::rowid_to_column() %>%
  filter(!(rowid %in% effector_rowids)) %>%
  rename(pathogen_short = `Pathogen species`) %>%
  select(-rowid)
