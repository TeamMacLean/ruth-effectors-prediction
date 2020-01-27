library(tidyverse)
library(magrittr)

# function to read the results and get the list of row index
blast_results <- function(result_path){

  # Read the results and turn it into dataframe
  # qseqid means Query Seq-id
  # qlen means Query sequence length
  # sseqid means Subject Seq-id
  # slen means Subject sequence length
  # length means Alignment length
  # nident means Number of identical matches
  # mismatch means Number of mismatches
  # positive means Number of positive-scoring matches
  df_results <- data.table::fread(result_path) %>%
    setNames(c("qseqid", "qlen", "sseqid", "slen", "length", "nident", "mismatch", "positive")) %>%
    rowwise() %>%
    mutate(
      percent_indentical = (nident/max(qlen, slen))*100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive/max(qlen, slen))*100 # The percentage of positive sequence over the longer sequence
    ) %>%
    ungroup()

  # Get the row indices of the source data for all of the identical percentage > 90%
  source_index_list_to_remove <- df_results %>%
    filter(percent_indentical > 90) %>%
    select(sseqid) %>%
    unique() %>%
    unlist()

  # Get the row indices of the query data for all of the identical percentage > 90%
  query_index_list_to_remove <- df_results %>%
    filter(percent_indentical > 90) %>%
    select(qseqid) %>%
    unique() %>%
    unlist()

  # Make list of all of the values
  list_results <- list(
    df = df_results,
    source_index = source_index_list_to_remove,
    query_index = query_index_list_to_remove
  )

  # Return the lists
  return(list_results)
}

# Read the data
ncbi_new_effector <- blast_results("data/BLAST-data/blast_effector1_effetor2.tsv")
df_effector1_x_effector2 <- ncbi_new_effector[["df"]]


# Function to concatenate the data
concat_ordered_cols <- function(data, var1, var2) {
  nrows <- nrow(data)
  unique_comb = rep(NA, nrows)

  var1_col <- as.character(data[[var1]])
  var2_col <- as.character(data[[var2]])

  for (i in 1:nrows) {
    unique_comb[[i]] <- sort(c(var1_col[i], var2_col[i])) %>% paste(collapse = "")
  }

  return(unique_comb)
}


# Delete self-comparisons
df_effector1_x_effector2 <- df_effector1_x_effector2 %>%
  filter(qseqid != sseqid)

df_effector1_x_effector2_with_comb <- df_effector1_x_effector2 %>%
  mutate(
    unique_id_comb = concat_ordered_cols(df_effector1_x_effector2, "qseqid", "sseqid")
  )

# Concatenate the data the function to the data
df_effector1_x_effector2_uniq <- df_effector1_x_effector2 %>%
  mutate(
    unique_id_comb = concat_ordered_cols(df_effector1_x_effector2, "qseqid", "sseqid")
  ) %>%
  group_by(unique_id_comb) %>%
  slice(1)

# Explore duplicated data
df_effector1_x_effector2 %>%
  mutate(
    unique_id_comb = concat_ordered_cols(df_effector1_x_effector2, "qseqid", "sseqid")
  ) %>%
  group_by(unique_id_comb) %>%
  summarise(count = n()) %>%
  filter(count > 2)



