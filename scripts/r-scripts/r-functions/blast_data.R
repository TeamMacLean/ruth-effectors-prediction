# Function for getting the data blast-ed

library(tidyverse)
library(Biostrings)
library(seqRFLP)

get_fasta_from_df <- function(df, column_id, column_seq, label = NULL, fasta_name = NULL, dir_path = "data/secreted_data/split-blast/fasta_files") {
  #' Function to automatically change the dataframe to fasta data
  #'
  #' @param df dataframe. The dataframe we want to change to fasta data (need to be R dataframe not tibble)
  #' @param dir_path path string. String contains path where we will save the data

  if (is.null(fasta_name)) {
    df_name <- deparse(substitute(df))
  } else {
    df_name <- fasta_name
  }

  # Change the label to become ID name
  if (!is.null(label)) {
    df <- df %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        {{ column_id }} := stringr::str_c({{ column_id }}, label, sep = "_")
      ) %>%
      dplyr::select(
        {{ column_id }}, {{ column_seq }}
      )
  } else {
    df <- df %>%
      # dplyr::ungroup() %>%
      dplyr::select(
        {{ column_id }}, {{ column_seq }}
      )
  }

  data_fa <- df %>%
    as.data.frame() %>%
    seqRFLP::dataframe2fas(file = paste0(dir_path, "/", df_name, ".fasta"))
  message("The data frame has been saved in ", paste0(dir_path, "/", df_name, ".fasta"))
}

get_blast_data <- function(database_fasta_path, query_fasta_path, dir_path = "data/secreted_data/split-blast/blast_files") {
  #' Blast data using the
  #'
  #' @param database_fasta_path dataframe. The path of the data for the dataframe for the database
  #' @param query_fasta_path dataframe. The path of the data for the dataframe for the query
  #' @param dir_path dataframe. The path where the directory we want to save

  # function to get the actual name
  get_name <- function(path1, path2) {
    list_name <- c(path1, path2) %>%
      stringr::str_split("/") %>%
      purrr::map(
        .f = function(x) {
          x[[length(x)]] %>%
            stringr::str_remove_all(".fasta")
        }
      ) %>%
      unlist()

    return(list_name)
  }

  # List of data to blast
  db_name <- get_name(database_fasta_path, query_fasta_path)[1]
  query_name <- get_name(database_fasta_path, query_fasta_path)[2]
  # result_name <- here::here(dir_path, paste0(db_name, "_vs_", query_name, ".tsv"))
  result_name <- paste0(dir_path, "/", db_name, "_vs_", query_name, ".tsv")

  # Making the database
  system(paste("makeblastdb ", "-in ", database_fasta_path, "-dbtype ", "prot"), intern = TRUE)

  # Blast the database against the query name
  system(paste("blastp ", "-query ", query_fasta_path, "-db ", database_fasta_path, "-out ", result_name, " -outfmt ", "\"6  qseqid qlen sseqid slen length nident mismatch positive\""))
  # message("The BLAST results have been saved in ", result_name %>% stringr::str_remove_all("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/"))
}


# function to read the results and get the list of row index
blast_results <- function(result_path, percent_threshold) {
  #' Function to read all of the data
  #'
  #' @param result_path path string. Path of blast result done using function get_blast_data()
  #' @param percent_threshold integer. Certain percentage (%) of the threshold of identical protein

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
    dplyr::mutate(
      percent_indentical = (nident / max(qlen, slen)) * 100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive / max(qlen, slen)) * 100 # The percentage of positive sequence over the longer sequence
    )

  # Get the data frame where the percent identical > 90
  df_identical_protein <- df_results %>%
    filter(percent_indentical > percent_threshold)

  # Get the row indices of the subject data for all of the identical percentage > 90%
  subject_index_list_to_remove <- df_results %>%
    filter(percent_indentical > percent_threshold) %>%
    select(sseqid) %>%
    unique() %>%
    unlist()

  # Get the row indices of the query data for all of the identical percentage > 90%
  query_index_list_to_remove <- df_results %>%
    filter(percent_indentical > percent_threshold) %>%
    select(qseqid) %>%
    unique() %>%
    unlist()

  # Make list of all of the values
  list_results <- list(
    df = df_results,
    df_identical_protein = df_identical_protein,
    subject_index = subject_index_list_to_remove,
    query_index = query_index_list_to_remove
  )

  # Return the lists
  return(list_results)
}
