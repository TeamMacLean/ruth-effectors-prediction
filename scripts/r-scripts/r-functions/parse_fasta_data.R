library(tidyverse)


parse_fasta_data_uniprot <- function(file_path) {
  # Read FASTA file
  fasta_data <- seqinr::read.fasta(file_path)
  # Number of entries
  num_data <- fasta_data %>% length()


  # Create empty data frame
  parsed_data <- data.frame(
    protein_id = rep(NA, num_data),
    pathogen = rep(NA, num_data),
    sequence = rep(NA, num_data)
  )

  for (i in 1:num_data) {
    # Read 'Annot' attribute and parse the string between 'OS=' and 'OX='
    pathogen <- fasta_data[[i]] %>%
      attr("Annot") %>%
      sub(".*OS= *(.*?) *OX=.*", "\\1", .)

    protein_id <- fasta_data[[i]] %>%
      attr("name") %>%
      strsplit("\\|") %>%
      unlist() %>%
      .[[2]]

    # Concatenate the vector of the sequence into a single string
    sequence <- fasta_data[[i]] %>%
      as.character() %>%
      toupper() %>%
      paste(collapse = "")

    # Input values into data frame
    parsed_data[i,] <- cbind(protein_id, pathogen, sequence)
  }

  return(parsed_data)
}

parse_fasta_data_ncbi <- function(file_path) {
  # Read FASTA file
  fasta_data <- seqinr::read.fasta(file_path)
  # Number of entries
  num_data <- fasta_data %>% length()


  # Create empty data frame
  parsed_data <- data.frame(
    protein_id = rep(NA, num_data),
    protein_fun = rep(NA, num_data),
    pathogen = rep(NA, num_data),
    sequence = rep(NA, num_data)
  )

  for (i in 1:num_data) {
    # Read 'Annot' attribute and parse the string between 'OS=' and 'OX='
    pathogen <- fasta_data[[i]] %>%
      attr("Annot") %>%
      sub(".*\\[ *(.*?) *\\].*", "\\1", .)

    protein_id <- fasta_data[[i]] %>%
      attr("name")

    protein_fun <- fasta_data[[i]] %>%
      attr("Annot") %>%
      stringr::str_remove(protein_id) %>%
      sub(".*> *(.*?) *\\[.*", "\\1", .)

    # Concatenate the vector of the sequence into a single string
    sequence <- fasta_data[[i]] %>%
      as.character() %>%
      toupper() %>%
      paste(collapse = "")

    # Input values into data frame
    parsed_data[i,] <- cbind(protein_id, protein_fun, pathogen, sequence)
  }

  return(parsed_data)
}


# Run function ---------------------------------------------

# Paths
uniprot_path     <- "~/Downloads/Cintaku/TSL/R/uniprot-data-mapped.fasta"
salmonella_path  <- "~/Downloads/Cintaku/TSL/R/salmonella.fasta"
xanthomonas_path <- "~/Downloads/Cintaku/TSL/R/xanthomonas.fasta"

# Parsed data
uniprot_parsed     <- parse_fasta_data_uniprot(uniprot_path)
salmonella_parsed  <- parse_fasta_data_ncbi(salmonella_path)
xanthomonas_parsed <- parse_fasta_data_ncbi(xanthomonas_path)

# Save data frames into CSV files
# write.csv(uniprot_parsed , "~/Downloads/uniprot-data-mapped.csv", row.names = FALSE)

