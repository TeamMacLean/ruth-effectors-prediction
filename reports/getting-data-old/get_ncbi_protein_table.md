Downloading protein data from the whole genome from NCBI
========================================================

Main function
-------------

``` r
download_protein_table <- function(species, strain, verbose = TRUE) {
  # URL of the queries: 
  #   https://www.ncbi.nlm.nih.gov/genome/?term=<QUERY>
  # URL of the Genome Assembly and Annotation report: 
  #   https://www.ncbi.nlm.nih.gov/genome/genomes/<GNOME_ID>
  # URL of the protein data:
  #   https://www.ncbi.nlm.nih.gov/genome/genomes/<GNOME_ID>?genome_assembly_id=<ASSEMBLY_ID>
  
  # Prepare species and strain for querying
  species_query <- species %>%
    stringr::str_replace_all(" ", "+")
  strain_query <- strain %>%
    stringr::str_replace_all(" ", "+") %>%
    stringr::str_replace_all("\\&", "\\\\\\&") %>%
    stringr::str_replace_all("\\(", "\\\\\\(") %>%
    stringr::str_replace_all("\\)", "\\\\\\)")

  # Create download directory if needed
  protein_dir <- "proteins"
  if (!dir.exists(protein_dir)) {
    system(paste("mkdir", protein_dir))
  }

  # Table filename
  table_name <- paste0("protein_table_", species_query %>% tolower(), "_", strain_query %>% tolower(), ".tsv")
  table_path <- paste0(protein_dir, "/", table_name)
  
  # Standard output header
    if (verbose) cat("\n--------------------------------------------------")
    if (verbose) cat(paste("\nSpecies:", species, "|" ,"Strain:", strain, "\n"))

  if (file.exists(table_path)) {
    if (verbose) cat("The file already exists\n")
    return(1)
  } else {
    # Base URL and queries
    base_url <- "https://www.ncbi.nlm.nih.gov"
    search_query <- "/genome/?term="
    protein_query <- paste(species_query, strain_query, sep = "+")

    # Temporary file
    dest_file <- tempfile()

    # Download HTML file
    query_url <- paste0(base_url, search_query, protein_query)

    query_url %>%
      paste0("wget ", ., " --output-document=", dest_file, " --quiet") %>%
      system()

    # Read file
    raw_html_file <- dest_file %>%
      readLines(warn = FALSE)

    # Get genome ID
    genome_id <- raw_html_file %>%
      grep("Genome Assembly and Annotation report", ., value = TRUE) %>%
      .[1] %>%
      stringr::str_split("\"") %>%
      unlist() %>%
      grep("genomes", ., value = TRUE) %>%
        stringr::str_replace("\\?", "") %>%
        stringr::str_split("/") %>%
        unlist() %>%
        .[length(.)]

    # Check if the data can be obtained
    if (length(genome_id) == 0) {
      if (verbose) cat("The species couldn't be found\n")
      return(0)
    } else {
      # Get genome assembly ID
      assembly_id <- raw_html_file %>%
        grep("genome_assembly_id", ., value = TRUE) %>%
        .[1] %>%
        stringr::str_split("\"") %>%
        unlist() %>%
        grep("genome_assembly_id", ., value = TRUE) %>%
        .[1] %>%
        stringr::str_split("genome_assembly_id") %>%
        unlist() %>%
        .[2] %>%
        stringr::str_replace("=", "")

      # URL and queries for downloading table
      table_base_url <- "https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi\\?action\\=GetFeatures4Grid\\&amp\\;download\\=on\\&amp\\;type\\=Proteins\\&amp\\;genome_id\\="
      table_query <- "\\&mode\\=2\\&is_locus\\=1\\&is_locus_tag\\=1\\&optcols\\=1,1,0,0,0"
      table_url <- paste0(table_base_url, genome_id, "\\&genome_assembly_id\\=", assembly_id, table_query)

      # Download table
      if (verbose) cat(paste0("Downloading (genome: ", genome_id, ", assembly: ", assembly_id, ")...\n"))

      table_url %>%
        paste0("wget ", ., " --output-document=", table_path, " --quiet") %>%
        system()
      
      if (verbose) cat(paste("Saved table as", table_name, "\n"))

      return(1)
    }
  }
}
```

Usage of the function
---------------------

### Simple example

``` r
species <- "Aeromonas salmonicida"
strain <- "A449"

download_protein_table(species, strain)
```

### Download data using the function

``` r
effector_table_raw <- data.table::fread("effector_sequence_with_metadata_uniq_to_download.csv") %>% 
  `colnames<-`(c("num", "species", "strain"))

effector_table_results <- effector_table_raw %>%
  dplyr::slice(1:10) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(status = download_protein_table(species, strain, verbose = FALSE))

effector_table_results
```

    ## Source: local data frame [10 x 4]
    ## Groups: <by row>
    ## 
    ## # A tibble: 10 x 4
    ##      num species                 strain             status
    ##    <int> <chr>                   <chr>               <dbl>
    ##  1     1 Acinetobacter baumannii ATCC 17978              1
    ##  2     2 Aeromonas hydrophila    SSU                     0
    ##  3     3 Aeromonas salmonicida   A449                    1
    ##  4     4 Beauveria bassiana      ARSEF2860               1
    ##  5     5 Blumeria graminis       f. sp. Hordei           1
    ##  6     6 Blumeria graminis       f. sp. hordei DH14      1
    ##  7     7 Botrytis cinerea        B05.10                  1
    ##  8     8 Botrytis cinerea        no data found           0
    ##  9     9 Brucella abortus        2308                    1
    ## 10    10 Burkholderia glumae     106619                  1

Check results after running function
------------------------------------

This allows us to check the status of the downloaded files and make sure
they are not empty.

``` r
effector_table_check <- effector_table_raw %>%
  # Recreate table path from species and strain
  dplyr::mutate(
    species_query = species %>%
      stringr::str_replace_all(" ", "+") %>%
      tolower(),
    strain_query = strain %>%
      stringr::str_replace_all(" ", "+") %>%
      stringr::str_replace_all("\\&", "\\\\\\&") %>%
      stringr::str_replace_all("\\(", "\\\\\\(") %>%
      stringr::str_replace_all("\\)", "\\\\\\)") %>%
      tolower(),
    table_path = paste0("protein_table_", species_query %>% tolower(), "_", strain_query %>% tolower(), ".tsv")
  ) %>%
  select(num, species, strain, table_path) %>%
  # Check if file exists
  dplyr::mutate(
    status = ifelse(file.exists(paste0("proteins/", table_path)), 1, 0),
    table_path = ifelse(status == 1, table_path, NA)
  ) %>%
  # Check contents of the file
  dplyr::rowwise() %>%
  dplyr::mutate(
    lines = ifelse(
      status == 1,
      system(paste0("wc -l ", "proteins/", table_path), intern = TRUE) %>%
        stringr::str_split(" ") %>%
        unlist() %>%
        .[1],
      NA
    ),
    # Substract header line
    lines = as.numeric(lines) - 1,
    # Mark status = 2 for empty downloaded files
    status = ifelse(status == 1 & lines == 0, 2, status)
  ) %>% 
  select(num, species, strain, status, table_path, lines)
```

### Filter files that are empty

These files are empty because although they genome ID and the genome
assembly ID are defined and can be queried, there is no protein table on
the database. This results in downloading a file with just the header,
but no data.

``` r
file_ids_empty <- effector_table_check %>% 
  dplyr::filter(status == 2) %>% 
  select(-status)

file_ids_empty
```

    ## Source: local data frame [0 x 5]
    ## Groups: <by row>
    ## 
    ## # A tibble: 0 x 5
    ## # … with 5 variables: num <int>, species <chr>, strain <chr>,
    ## #   table_path <chr>, lines <dbl>

In total, there are `nrow(file_ids_empty)` empty files.

### Filter unobtainable proteins

Likewise, we can check those species that couldn’t be queried.

``` r
file_ids_missing <- effector_table_check %>% 
  dplyr::filter(status == 0) %>% 
  select(-status)

file_ids_missing
```

    ## Source: local data frame [90 x 5]
    ## Groups: <by row>
    ## 
    ## # A tibble: 90 x 5
    ##      num species               strain                      table_path lines
    ##    <int> <chr>                 <chr>                       <chr>      <dbl>
    ##  1     2 Aeromonas hydrophila  SSU                         <NA>          NA
    ##  2     8 Botrytis cinerea      no data found               <NA>          NA
    ##  3    15 Cercospora apii       no data found               <NA>          NA
    ##  4    16 Cercospora beticola   no data found               <NA>          NA
    ##  5    18 Clavibacter michigan… Cm15-2.0 sm                 <NA>          NA
    ##  6    19 Colletotrichum orbic… 104-T                       <NA>          NA
    ##  7    20 Coxiella burnetii     Nine Mile (RSA493)          <NA>          NA
    ##  8    21 Coxiella burnetii     Nine Mile Phase II Clone 4… <NA>          NA
    ##  9    22 Coxiella burnetii     NMII                        <NA>          NA
    ## 10    25 Dothistroma septospo… no data found               <NA>          NA
    ## # … with 80 more rows

In total, there are `nrow(file_ids_missing)` unobtainable species.

Get list of protein products
----------------------------

With this command, we can get all protein products in a single file
using `bash`. This file can be later used in Batch Entrez
(<a href="https://www.ncbi.nlm.nih.gov/sites/batchentrez" class="uri">https://www.ncbi.nlm.nih.gov/sites/batchentrez</a>)
to get the fasta sequence data.

``` bash
tail -n+2 proteins/protein_table_*.tsv | cut -d$'\t' -f9 | uniq | grep -v 'protein_table' > proteins/proteins_list.csv
```

In total, there are unique proteins.

### Downloading the FASTA data

Given the amount of data, this is nearly impossible to process within
the Batch Entrez website, as the `proteins_list.csv` would need to be
split into files of 1000 lines, as that’s the maximum entries Batch
Entrez will process at a time. Even if we decided to do that, the FASTA
data can only be downloaded 200 entries at a time. This would result in
manually “downloading” (i.e. copying and pasting the “FASTA text”
results) roughly 0 times.

Therefore, this will need to be done with BLAST+.

For this we need to download the NCBI non-redundant protein database
(<a href="ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz" class="uri">ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz</a>)
and load it into our BLAST+ installation.

Once we have the database (it’s still downloading), we should be able to
make a bash script that parses the individual protein products like
shown bellow.

``` bash
blastdbcmd -db nr -entry WP_042791309.1 -outfmt "%f"
```
