# read the result of BLAST
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
      percent_identical = (nident/max(qlen, slen))*100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive/max(qlen, slen))*100 # The percentage of positive sequence over the longer sequence
      ) %>%
    ungroup()

  # Get the row indices of the source data for all of the identical percentage > 90%
  source_index_list_to_remove <- df_results %>%
    filter(percent_identical > 90) %>%
    select(sseqid) %>%
    unique() %>%
    unlist()

  # Get the row indices of the query data for all of the identical percentage > 90%
  query_index_list_to_remove <- df_results %>%
    filter(percent_identical > 90) %>%
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


# Function to remove from the list

remove_from_data_sets <- function(input_seq_path, label_seq_path, drop){

  # Read the data from the .csv file
  df_input <- data.table::fread(input_seq_path, header = FALSE) %>%
    setNames("sequence")
  df_label <- data.table::fread(label_seq_path, header = FALSE) %>%
    setNames("label")

  # Combine the input data and the label, then drop the rows based on the list
  df_new <- df_input %>%
    cbind(df_label) %>%
    filter(!row_number() %in% drop)

  # Get the information about the removed columns
  removed_rows <- df_label %>%
    filter(row_number() %in% drop) %>%
    select(label) %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("label", "freq")) %>%
    mutate(label = ifelse(label == 1, "effector", "noneffector"))

  # Create list for the results
  results_list <- list(
    df = df_new,
    removed = removed_rows
  )

  return(results_list)
}


# Function to get information from all data
get_info <- function(blast_train_x_val_path,
                     blast_train_x_test_path,
                     input_train_path,
                     label_train_path,
                     input_val_path,
                     label_val_path,
                     input_test_path,
                     label_test_path){

  # Getting the results from comparing validation dataset and training dataset
  blast_results_train_x_val <- blast_results(blast_train_x_val_path)
  df_train_x_val <- blast_results_train_x_val[["df"]]
  source_index_train_x_val <- blast_results_train_x_val[["source_index"]]
  query_index_train_x_val <- blast_results_train_x_val[["query_index"]]

  # Getting the results from comparing testing dataset and training dataset
  blast_results_train_x_test <- blast_results(blast_train_x_test_path)
  df_train_x_test <- blast_results_train_x_test[["df"]]
  source_index_train_x_test <- blast_results_train_x_test[["source_index"]]
  query_index_train_x_test <- blast_results_train_x_test[["query_index"]]

  # Remove all of the rows of the training data
  # Since there might be intersection between index on the training when comparing them with validation and testing sets,
  # we need to find the intesections

  intersec_training_index <- c(source_index_train_x_val, source_index_train_x_test) %>%
    unique() %>%
    unlist()

  results_after_removed_training <- remove_from_data_sets(input_train_path, label_train_path, intersec_training_index)
  removed_rows_training <- results_after_removed_training[["removed"]] %>%
    mutate(type = "training")

  results_after_removed_validation <- remove_from_data_sets(input_val_path, label_val_path, query_index_train_x_val)
  removed_rows_val <- results_after_removed_validation[["removed"]] %>%
    mutate(type = "validation")

  results_after_removed_testing <- remove_from_data_sets(input_test_path, label_test_path, query_index_train_x_test)
  removed_rows_test <- results_after_removed_testing[["removed"]] %>%
    mutate(type = "testing")

  all_removed_rows <- removed_rows_training %>%
    rbind(., removed_rows_val) %>%
    rbind(., removed_rows_test)

  all_ <- reshape2::dcast(all_removed_rows, label ~ type, value.var = "freq") %>%
    as.data.frame()

  results_list <- list(
    df_val = blast_results_train_x_val,
    df_test = blast_results_train_x_test,
    all = all_
  )

  return(results_list)
}


blast_train_x_val_path = "data/BLAST-data/blast_train_x_val.tsv"
blast_train_x_test_path = "data/BLAST-data/blast_train_x_test.tsv"
input_train_path = "data/BLAST-data/blast_train.csv"
label_train_path = "data/BLAST-data/blast_label_train.csv"
input_val_path = "data/BLAST-data/blast_val.csv"
label_val_path = "data/BLAST-data/blast_label_val.csv"
input_test_path = "data/BLAST-data/blast_test.csv"
label_test_path = "data/BLAST-data/blast_label_test.csv"


all_results <-  get_info(blast_train_x_val_path,
                                     blast_train_x_test_path,
                                     input_train_path,
                                     label_train_path,
                                     input_val_path,
                                     label_val_path,
                                     input_test_path,
                                     label_test_path)


summary <- all_results[["all"]]


# comparing effector data

ncbi_new_effector <- blast_results("data/BLAST-data/blast_effector1_effetor2.tsv")

df_effector1_x_effector2 <- ncbi_new_effector[["df"]]

diff_values <- df_effector1_x_effector2 %>%
  filter(., qseqid != sseqid)

same_values <- df_effector1_x_effector2 %>%
  filter(., qseqid == sseqid)

effector2 <- ncbi_new_effector[["query_index"]]

effector1 <- ncbi_new_effector[["source_index"]]

intersec_effector <- intersect(effector1, effector2)


