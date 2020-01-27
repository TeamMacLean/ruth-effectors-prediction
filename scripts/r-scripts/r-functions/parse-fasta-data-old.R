library(seqinr)
library(tidyverse)

# Read FASTA file
fasta_data <- seqinr::read.fasta("data/uniprot-data-mapped.fasta")

# Number of entries
num_data <- fasta_data %>% length()


# Create empty data frame
parsed_data <- data.frame(
  pathogen = rep(NA, num_data),
  sequence = rep(NA, num_data)
)

for (i in 1:num_data) {
  # Read 'Annot' attribute and parse the string between 'OS=' and 'OX='
  pathogen <- fasta_data[[i]] %>%
    attr("Annot") %>%
    sub(".*OS= *(.*?) *OX=.*", "\\1", .)

  # Concatenate the vector of the sequence into a single string
  sequence <- fasta_data[[i]] %>%
    as.character() %>%
    toupper() %>%
    paste(collapse = "")

  # Input values into data frame
  parsed_data[i,] <- cbind(pathogen, sequence)
}

# Save data frame into CSV file
write.csv(parsed_data, "data/uniprot-data-mapped.csv", row.names = FALSE)

