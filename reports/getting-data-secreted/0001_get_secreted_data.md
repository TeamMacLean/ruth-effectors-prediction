Getting secreted data
=====================

Introduction
------------

Ideally, the heatmap plot made in jupyter lab book can show which
positions of the sequences activated the most by convolutional layers,
however after grouping the data based on the pathogen organisms (it
still has not shows similar pattern for each organism). Hypothesis: This
may due to the effector data which is not secreted by the organism.

Aim
---

In this report, several steps will be conducted in order to get the
protein data secreted from each pathogen organisms (to improve the
perfomance of the model in recognising the pattern).

Execution
---------

Several steps will be done in order to achieve the aims above get the
all of protein sequences from the same genomic sequence, and use signalP
5.0 to see the probability of those protein to have signal peptide,
which means the protein is secreted by each species.

``` r
effector_data <- data.table::fread("../../../data/getting-data-new/binary-class-data/effector_data.csv")
```

``` r
# Load phi base data
all_phi_base <- data.table::fread("../../../data/phi-base-current-data/phi-base_without_column_desc.csv")
```

``` r
list_protein_id <- effector_data %>%
  pull(ProteinID)
```

``` r
# Getting the information about the Gene names and Gene IDs, and the protein source
phi_base_prot <- all_phi_base %>%
  dplyr::filter(ProteinID %in% list_protein_id) %>%
  group_by(ProteinID, ProteinIDsource) %>%
  slice(1) %>%
  dplyr::select(ProteinIDsource, ProteinID, GeneIDsource, GeneID, Gene)

# Getting the dictinct value of each data
phi_base_prot %>%
  dplyr::distinct(Gene) %>%
  nrow()
```

    ## [1] 402

From the result, it seems that the they have unique name. NOw, using the
unique name we have above, we want to get the genome name.

The example of the genome I can get from gene is `KGO57615` ==&gt; input
in Gene, and then I can get the info of the Genome, and using that
Genome ID, retrive it in Genome NCBI. Get the whole protein sequences.
In this case, using the gene `LysM1` will not result in a specific one.

``` r
connected.to.internet <- function() {
  if (curl::curl_fetch_memory("www.google.com")$status_code == 200) {
    return(TRUE)
  } else {
    stop(
      "It seems that you are not connected to the internet. Could you please check?",
      call. = FALSE
    )
  }
}

connected.to.internet()
```

    ## [1] TRUE

``` r
# get_proteins_in_genome <- function(gene){
#   # URL of the queries:
#   #   https://www.ncbi.nlm.nih.gov/genome/?term=<QUERY>
#   # URL to get the Genome ID if the gene ID is known
#   #   https://www.ncbi.nlm.nih.gov/gene/?term=<GENE_ID>
#   # URL to get all protein data in all genome
#   #   https://www.ncbi.nlm.nih.gov/nuccore/<GNOME_ID>
#
#   # gene_name <- paste0(gene)
#   # Get the gene information
#   # Base URL and queries
#   base_url <- "https://www.ncbi.nlm.nih.gov"
#   search_gene_query <- "/gene?term="
#
#   # Download HTML file
#   search_gene_url <- paste0(base_url, search_gene_query, gene)
#   search_gene_html_file <- textreadr::read_html(search_gene_url)
#
#   # Get the genomic sequence ID
#   genomic_sequence_id <- search_gene_html_file %>%
#   stringr::str_detect("Genomic Sequence") %>%
#   which() %>%
#   `+`(1) %>%
#   search_gene_html_file[.]
#
#   # Get the all of sequence from the whole genomic sequence
#   search_whole_genome_query <- "/nuccore/"
#   search_whole_genome_url <- paste0(base_url, search_whole_genome_query, genomic_sequence_id)
#   search_whole_genome_html_file <- textreadr::read_html(search_whole_genome_url)
#
#   print(search_whole_genome_url)
#   # return(search_whole_genome_html_file)
# }
```

Function
--------

``` r
get_signalp_pred <- function(fasta_filename, verbose = FALSE) {
  signalp_path <- "/Users/kristian/Documents/Workspace/Software/signalp/bin"
  verbose_string <- tolower(deparse(substitute(verbose)))
  fasta_path <- fasta_filename

  command_string <- paste0(
    "export PATH=", signalp_path, ":$PATH ;",
    " signalp -fasta ", fasta_path, " -org euk -format short ",
    "-verbose=", verbose_string, " -stdout"
  )

  data_lines <- command_string %>%
    system(intern = TRUE)

  data <- data_lines %>%
    .[-c(1, 2)] %>%
    as.list() %>%
    data.frame() %>%
    t() %>%
    as_tibble() %>%
    tidyr::separate(col = V1, into = c("ID", "Prediction", "SP(Sec/SPI)", "OTHER", "CS Position"), sep = "\t")

  return(data)
}
```

``` r
# Function if protein ID has protein ID
check_if_genomic_seq_exist <- function(gene_term) {
  if (nchar(gene_term) > 0) {
    # Get the Gene ID from NCBI
    ncbi_gene_id <- entrez_search(db = "gene", term = gene_term)$ids

    # Return boolean whether the sequence exists or not
    seq_exists <- (length(ncbi_gene_id) != 0)

    return(seq_exists)
  } else {
    return(NA)
  }
}
```

``` r
df_check_genomic_seq <- phi_base_prot %>%
  ungroup() %>%
  rowwise() %>%
  dplyr::mutate("ncbi_genomic_seq_exist" = check_if_genomic_seq_exist(GeneID))
```

``` r
phi_base_prot %>%
  ungroup() %>%
  slice(1:20) %>%
  rowwise() %>%
  dplyr::mutate("ncbi_genomic_seq_exist" = check_if_genomic_seq_exist(GeneID))
```

    ## Source: local data frame [20 x 6]
    ## Groups: <by row>
    ## 
    ## # A tibble: 20 x 6
    ##    ProteinIDsource ProteinID  GeneIDsource  GeneID  Gene   ncbi_genomic_se…
    ##    <chr>           <chr>      <chr>         <chr>   <chr>  <lgl>           
    ##  1 Uniprot         A0A023UJQ9 Genbank       AHY021… Avr5(… FALSE           
    ##  2 Uniprot         A0A059UDR8 Genbank       KJ5712… BEC10… FALSE           
    ##  3 Uniprot         A0A0A0S3X0 Genbank       AIW074… AvrLm2 FALSE           
    ##  4 Uniprot         A0A0A2ILW0 Genbank       KGO576… LysM1  TRUE            
    ##  5 Uniprot         A0A0A2JB06 Genbank       KGO519… LysM4  TRUE            
    ##  6 Uniprot         A0A0A2JVC2 Genbank       KGO587… LysM2  TRUE            
    ##  7 Uniprot         A0A0A2JY92 Genbank       KGO597… LysM3  TRUE            
    ##  8 Uniprot         A0A0A7DLN4 Genbank       AIP869… Iug6   FALSE           
    ##  9 Uniprot         A0A0A7DM22 Genbank       AIP869… Iug9   FALSE           
    ## 10 Uniprot         A0A0D1C5E3 Genbank       XP_011… cce1 … TRUE            
    ## 11 Uniprot         A0A0D1E5A8 Genbank       KIS708… apB73  TRUE            
    ## 12 Uniprot         A0A0E3GIM3 Genbank       AJR197… HgGLA… FALSE           
    ## 13 Uniprot         A0A0L0VYD3 Ensembl Geno… PSTG_0… PST02… FALSE           
    ## 14 Uniprot         A0A0M5K865 Genbank       ALD519… PsCRN… FALSE           
    ## 15 Uniprot         A0A0N8SZV2 Genbank       KPY754… HopAI1 FALSE           
    ## 16 Uniprot         A0A0Q0BGR4 Genbank       KPY886… HopBB1 FALSE           
    ## 17 Uniprot         A0A169T708 Genbank       BAU983… AVR1   FALSE           
    ## 18 Uniprot         A0A1B0RFQ0 ""            ""      Avr4   NA              
    ## 19 Uniprot         A0A2P1C6A6 Genbank       AVJ356… Pst_8… FALSE           
    ## 20 Uniprot         A0A2R2Z552 Genbank       AUD400… Pc107… FALSE

``` r
get_all_prot_seqs <- function(gene_term) {

  # Get the GeneID given gene information form Phibase
  ncbi_gene_id <- entrez_search(db = "gene", term = gene_term)$ids

  if (length(ncbi_gene_id) == 0) {
    return(NULL)
  } else {

    # Get the genomic sequence ID
    ncbi_genomic_id <- entrez_link(dbfrom = "gene", id = ncbi_gene_id, db = "all")$links$gene_nuccore_pos

    # Get all of the protein sequence from given genomic sequence ID in NCBI
    ncbi_all_protein_id <- entrez_link(dbfrom = "nuccore", id = ncbi_genomic_id, db = "all")$links$nuccore_protein

    # Get all of the protein sequences
    prot_fasta <- entrez_fetch(db = "protein", id = ncbi_all_protein_id, rettype = "fasta")

    # Save in temp file
    temp <- tempfile()
    write(prot_fasta, temp)
    # parsed_fasta <- ape::read.FASTA(temp, type="AA")

    # Get the SignalP prediction from the fasta file we get
    pred <- get_signalp_pred(temp)

    # Get the information about fasta file if it is secreted
    pred_signalP_exist <- pred %>%
      dplyr::filter(Prediction == "SP(Sec/SPI)") %>%
      dplyr::mutate(Gene = gene_term)
  }
  return(pred_signalP_exist)
}
```

``` r
entrez_search(db = "gene", term = "KGO57615")$ids
```

    ## [1] "27674751"

``` r
get_all_prot_seqs("AHY02126")
```

    ## NULL
