library(seqinr)
library(tidyverse)

# Read FASTA file
fasta_data <- seqinr::read.fasta("data/uniprot-data-mapped.fasta")

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
  parsed_data[i, ] <- cbind(protein_id, pathogen, sequence)
}

# Save data frame into CSV file
write.csv(parsed_data, "data/uniprot-data-mapped.csv", row.names = FALSE)
