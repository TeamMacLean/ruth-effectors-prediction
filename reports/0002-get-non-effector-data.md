Getting the Non-effector data
=============================

Getting the data from NCBI
--------------------------

| Bacteria                  | Number |     | Bacteria                    | Number |
|---------------------------|:------:|-----|-----------------------------|--------|
| Acinetobacter baumannii   |    1   |     | Mycobacterium tuberculosis  | 1      |
| Aeromonas hydrophila      |    2   |     | Pantoea stewartii           | 2      |
| Aeromonas salmonicida     |    1   |     | Pseudomonas aeruginosa      | 2      |
| Brucella abortus          |    1   |     | Pseudomonas cichorii        | 1      |
| Burkholderia glumae       |    1   |     | Pseudomonas savastanoi      | 14     |
| Burkholderia pseudomallei |   22   |     | Pseudomonas syringae        | 42     |
| Campylobacter jejuni      |    1   |     | Ralstonia solanacearum      | 46     |
| Citrobacter rodentium     |    1   |     | Salmonella enterica         | 64     |
| Clavibacter michiganensis |    2   |     | Shigella flexneri           | 2      |
| Coxiella burnetii         |    3   |     | Staphylococcus aureus       | 1      |
| Cystobacter fuscus        |    1   |     | Vibrio parahaemolyticus     | 1      |
| Edwardsiella ictaluri     |    8   |     | Xanthomonas axonopodis      | 8      |
| Erwinia amylovora         |    9   |     | Xanthomonas campestris      | 16     |
| Escherichia coli          |    6   |     | Xanthomonas oryzae          | 22     |
| Francisella tularensis    |    1   |     | Xylella fastidiosa          | 3      |
| Helicobacter pylori       |    1   |     | Yersinia enterocolitica     | 1      |
| Legionella pneumophila    |    7   |     | Yersinia pestis             | 1      |
| Listeria monocytogenes    |    2   |     | Yersinia pseudotuberculosis | 3      |

``` r
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
```

``` r
get_different_sequence_samples <- function(data, samples, seq_chars = 10) {
  data <- data %>%
    mutate(short_sequence = stringr::str_sub(sequence, 1, seq_chars)) %>%
    group_by(short_sequence) %>%
    slice(1) %>%
    ungroup() %>%
    select(-short_sequence) %>% 
    slice(1:samples)
  
  return(data)
}
```

``` r
# Run function ---------------------------------------------

# Paths
ncbi_others_path <- "../../data/ncbi-others.fasta"
ncbi_fungus_path <- "../../data/ncbi-fungus.fasta"
ncbi_oomycetes_path <- "../../data/ncbi-oomycetes.fasta"
ncbi_bacteria_path <- "../../data/ncbi-bacteria.fasta"
ncbi_bacteria_salmonella_path <- "../../data/ncbi-bacteria-salmonella.fasta"
ncbi_bacteria_xanthomonas_path <- "../../data/ncbi-bacteria-xanthomonas.fasta"

# Parsed data
ncbi_others_parsed <- parse_fasta_data_ncbi(ncbi_others_path)
ncbi_fungus_parsed <- parse_fasta_data_ncbi(ncbi_fungus_path)
ncbi_oomycetes_parsed <- parse_fasta_data_ncbi(ncbi_oomycetes_path)
ncbi_bacteria_parsed <- parse_fasta_data_ncbi(ncbi_bacteria_path)
ncbi_bacteria_xanthomonas_parsed <- parse_fasta_data_ncbi(ncbi_bacteria_xanthomonas_path)
ncbi_bacteria_salmonella_parsed <- parse_fasta_data_ncbi(ncbi_bacteria_salmonella_path)

# get the data from Salmonella and Xanthomonas
unique_bacteria_xanthomonas_parsed <- get_different_sequence_samples(ncbi_bacteria_xanthomonas_parsed, 22)
unique_bacteria_salmonella_parsed <- get_different_sequence_samples(ncbi_bacteria_salmonella_parsed, 64)

# merge all of data by rows
ncbi_noneffector_parsed <- ncbi_bacteria_parsed %>% 
  rbind(., unique_bacteria_xanthomonas_parsed) %>% 
  rbind(., unique_bacteria_salmonella_parsed) %>% 
  rbind(., ncbi_fungus_parsed) %>% 
  rbind(., ncbi_oomycetes_parsed) %>% 
  rbind(., ncbi_others_parsed)

# Save data frames into CSV files
write.csv(ncbi_noneffector_parsed, "../../data/ncbi_noneffector_parsed.csv", row.names = FALSE)
```

``` r
# add label on the data frame 
noneffector <- ncbi_noneffector_parsed %>% 
  select(sequence) %>% 
  mutate(label = as.factor(0))

write.csv(noneffector, "../../data/noneffector.csv", row.names = FALSE)
```

View the effector and noneffector data
--------------------------------------

``` r
effector <- data.table::fread("../../data/effector.csv", header = TRUE)
noneffector <- data.table::fread("../../data/noneffector.csv", header = TRUE)
```

``` r
effector %>% 
  head(10) %>%
  mutate(sequence = substr(sequence, 1, 30)) %>%
  knitr::kable()
```

| sequence                       |  label|
|:-------------------------------|------:|
| MKLSLLSVELALLIATTLPLCWAAALPVGL |      1|
| MHYTTLLLSTLLVGTALAQPTNPPAKTPKK |      1|
| MVQFKTIFLSTALAALFSTGSSSPATKNNV |      1|
| MKFLVLPLSLAFLQIGLVFSTPDRCRYTLC |      1|
| MKFNKTIPLYILAFFSTAVIAGGRKWTNKV |      1|
| MKCNNIILPFALVFFSTTVTAGGGWTNKQF |      1|
| MLFNAAAAAVFAPLLVMGNVLPRNAGNSPG |      1|
| MNFRALFAATVAALVGSTSATTCTTSQQTV |      1|
| MLFYSLFFFHTVAISAFTNIGTFSHPVYDY |      1|
| MRDEMWNTATEPIAIIGSGCKFPGGSTTPS |      1|

``` r
noneffector %>% 
  head(10) %>%
  mutate(sequence = substr(sequence, 1, 30)) %>%
  knitr::kable()
```

| sequence                       |  label|
|:-------------------------------|------:|
| MPVLTGTAEANSTISIFDGTTLLGTTTADA |      0|
| MLPVPPGSIDDRTVLYTYNVLNQKTSMTRV |      0|
| MSALRDAMSHAALGWVKPELDETLRQARNE |      0|
| MNSKLYKLIFCRRLGCLIAVGEFTRSYGRA |      0|
| MISGAPSQDSLLPDNRHAADYQQLRERLIQ |      0|
| MQEAKVDGGYVKFVNLFKLSDLEQALLPYK |      0|
| MWWSKQKSRITEGQRTLAAPMIMSLEPRML |      0|
| MSSQEFQSIIAGLQHRQITPAEARQRLQRL |      0|
| MHHPSSAQPESASLEAVGDVSAARPASSVP |      0|
| MEIAMAAIWAQVLGIQRVGRQDNFFELGGH |      0|
