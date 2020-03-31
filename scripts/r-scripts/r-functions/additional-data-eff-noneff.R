# This scripts is used to get the additional effector and noneffector data based on the BLAST results previosuly obtained
# In order to make sure there is no identical protein within the .fatsa file, we can blast the data to themself

# The succesive order of the BLAST results of each the effector-noneffector data to itself:
# 1. Get the non-identical protein data
# 2. Take the sequences
# 3. Put together with the existing data

# read the result of BLAST
library(tidyverse)


# Function to concatenate the data
concat_ordered_cols <- function(data, var1, var2, sort = TRUE) {
  nrows <- nrow(data)
  unique_comb = rep(NA, nrows)

  var1_col <- data[[var1]] %>% as.character()
  var2_col <- data[[var2]] %>% as.character()

  if (sort) {
    for (i in 1:nrows) {
      unique_comb[[i]] <- sort(c(var1_col[i], var2_col[i])) %>% paste(collapse = "")
    }
  } else {
    for (i in 1:nrows) {
      unique_comb[[i]] <- c(var1_col[i], var2_col[i]) %>% paste(collapse = "")
    }
  }

  return(unique_comb)
}

# function to read the results and get the list of row index
blast_results <- function(result_path){

  # Read the results and turn it into dataframe
  # qseqid means Query Seq-id
  # qlen means Query sequence length
  # sseqid means Subject Seq-id
  # slen means Subject sequence length
  # length means Alignment length
  # nident means Number of identical matches
  df_results <- data.table::fread(result_path) %>%
    setNames(c("qseqid", "qlen", "sseqid", "slen", "length", "nident"))

  # Delete self-comparisons
  df_results <- df_results %>%
    dplyr::filter(qseqid != sseqid) %>%
    mutate(unique_id_comb = concat_ordered_cols(., "qseqid", "sseqid", sort = TRUE),
           unique_id_comb_unsorted = concat_ordered_cols(., "qseqid", "sseqid", sort = FALSE)
          ) %>%
    # Summarise nident per each identical unique_id_comb_unsorted value
    # the other variables are selected not to lose them, but they are not really the "grouping" variable,
    # as they are identical for each unique value of "unique_id_comb_unsorted"
    group_by(unique_id_comb_unsorted, unique_id_comb, qseqid, sseqid, qlen, slen) %>%
    summarise(nident = sum(nident)) %>%
    # Select only one example per combination
    group_by(unique_id_comb) %>%
    slice(1) %>%
    # Remove auxiliary columns
    ungroup() %>%
    select(-unique_id_comb, -unique_id_comb_unsorted) %>%
    rowwise() %>%
    mutate(
      percent_identical = (nident/max(qlen, slen))*100 # The percentage of identical sequence over the longer sequence
    ) %>%
    ungroup() %>%
    dplyr::filter(nident < 90)

  get_all_IDs <- c(df_results[["qseqid"]], df_results[["sseqid"]]) %>%
    unique()

  result_list <- list(
    df = df_results,
    uniq_pro_IDs = get_all_IDs
  )

  return(result_list)
}


blast_effector_bacteria <- blast_results("data/BLAST-data/additional-data/0001-blast-effector-bacteria.tsv")
blast_effector_fungi <- blast_results("data/BLAST-data/additional-data/0002-blast-effector-fungi.tsv")
blast_effector_oomycetes <- blast_results("data/BLAST-data/additional-data/0003-blast-effector-oomycetes.tsv")
blast_noneffector_bacteria <- blast_results("data/BLAST-data/additional-data/0004-blast-noneffector-bacteria.tsv")
blast_noneffector_fungi <- blast_results("data/BLAST-data/additional-data/0005-blast-noneffector-fungi.tsv")
blast_noneffector_oomycetes <- blast_results("data/BLAST-data/additional-data/0006-blast-noneffector-oomycetes.tsv")

# get the IDs of non-identical protein

eff_bacteria_IDs <- blast_effector_bacteria[["uniq_pro_IDs"]]
eff_fungi_IDs <- blast_effector_fungi[["uniq_pro_IDs"]]
eff_oomycetes_IDs <- blast_effector_oomycetes[["uniq_pro_IDs"]]
noneff_bacteria_IDs <- blast_noneffector_bacteria[["uniq_pro_IDs"]]
noneff_fungi_IDs <- blast_noneffector_fungi[["uniq_pro_IDs"]]
noneff_oomycetes_IDs <- blast_noneffector_oomycetes[["uniq_pro_IDs"]]

# get the sample with the  number that we need for the additional data

IDs_batch_entrez_effector <- c(sample(eff_bacteria_IDs, size = 27),
                               sample(eff_fungi_IDs, size = 11),
                               sample(eff_oomycetes_IDs, size = 6))

IDs_batch_entrez_noneffector <- c(sample(noneff_bacteria_IDs, size = 122),
                                  sample(noneff_fungi_IDs, size = 2),
                                  sample(noneff_oomycetes_IDs, size = 16))

write.table(IDs_batch_entrez_effector, "data/BLAST-data/additional-data/IDs_batch_entrez_effector.csv", row.names = FALSE, quote = FALSE)
write.table(IDs_batch_entrez_noneffector, "data/BLAST-data/additional-data/IDs_batch_entrez_noneffector.csv", row.names = FALSE, quote = FALSE)

saveRDS(IDs_batch_entrez_effector, "data/BLAST-data/additional-data/IDs_batch_entrez_effector.RDS")
saveRDS(IDs_batch_entrez_noneffector, "data/BLAST-data/additional-data/IDs_batch_entrez_noneffector.RDS")


