getting_sample_data <- function(tbl,
                                col_organism,
                                col_id,
                                col_seq,
                                organism_name,
                                num_sample) {

  # Getting the bag of sequence data only for the given organism name
  table_processed <- tbl %>%
    # Remove the duplicate data
    dplyr::filter({{ col_organism }} == organism_name) %>%
    dplyr::group_by({{ col_seq }}) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # Initialize list of data frames
  list_of_df <- list()
  max_iterations <- nrow(table_processed)
  message("The loop will run for a maximum of ", max_iterations, " iterations.")

  # Add first sequence to list
  list_of_df[["seq1"]] <- table_processed %>%
    dplyr::sample_n(size = 1, replace = FALSE)

  table_processed <- table_processed %>%
    dplyr::filter(
      {{ col_id }} != list_of_df[["seq1"]] %>% pull({{ col_id }})
    )

  # Create temporary folder
  temp_dir <- tempdir()

  # Setting the condition of the sample availability
  for (i in 2:max_iterations) {

    # Exit loop if enough samples have been found
    if (length(list_of_df) >= num_sample) {
      message("Found enough samples.")
      break
    }

    message("\n", "Trying to get sample number ", length(list_of_df) + 1, " out of ", num_sample)
    message("Global iteration ", i, " out of ", max_iterations)

    # Add to list of data frames
    list_of_df[[paste0("seq", i)]] <- table_processed %>%
      dplyr::sample_n(size = 1, replace = FALSE)

    # Update the table by removing the element that is taken already
    table_processed <- table_processed %>%
      dplyr::filter(
        {{ col_id }} != list_of_df[[paste0("seq", i)]] %>% pull({{ col_id }})
      )

    # FASTA from previous sequences
    get_fasta_from_df(
      df = list_of_df[seq(1, length(list_of_df) - 1, 1)] %>% purrr::reduce(rbind),
      column_id = {{ col_id }},
      column_seq = {{ col_seq }},
      fasta_name = "fasta_1",
      dir_path = temp_dir
    )

    # FASTA from new sequence
    get_fasta_from_df(
      df = list_of_df[[paste0("seq", i)]],
      column_id = {{ col_id }},
      column_seq = {{ col_seq }},
      fasta_name = "fasta_2",
      dir_path = temp_dir
    )

    # BLAST sequences
    get_blast_data(
      database_fasta_path = paste0(temp_dir, "/fasta_1.fasta"),
      query_fasta_path = paste0(temp_dir, "/fasta_2.fasta"),
      dir_path = temp_dir
    )

    res_list <- blast_results(
      result_path = paste0(temp_dir, "/fasta_1_vs_fasta_2.tsv"),
      percent_threshold = 95
    )

    percent <- res_list[["df"]] %>%
      dplyr::filter(qseqid != sseqid) %>%
      pull(percent_indentical)

    message("Checking percentages...")
    if (any(percent > 90)) {
      message("Percentages are: ", paste(round(percent, 2), collapse = ", "), ". Deleting from list_of_df.")
      list_of_df[[paste0("seq", i)]] <- NULL
    }

    # Delete temporary folder
    # unlink(temp_dir)

    # Notify about reaching the maximum possible iterations
    if (i == max_iterations) {
      message("It was not possible to find enough samples.")
    }
  }

  # Reduce lists to dataframes
  df_combined <- list_of_df %>%
    purrr::reduce(rbind)

  return(df_combined)
}
