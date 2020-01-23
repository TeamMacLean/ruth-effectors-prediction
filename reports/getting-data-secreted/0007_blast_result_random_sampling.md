BLAST the result of random sampling data to the original effector data
======================================================================

Background
----------

In this report, the process to show how the process to blast the result
of random sampling from the SignalP prediction is shown.

Load libraries
--------------

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(seqRFLP)
library(stringr)
```

Load data
---------

``` r
effector_data_info <- readRDS("effector_data_info.RDS")

effector_data <- data.table::fread("../../../data/getting-data-new/binary-class-data/effector_data.csv")

effector_seq <- left_join(effector_data_info %>% 
  dplyr::select(c(ProteinID, class)), 
                effector_data %>% 
                  dplyr::select(c(ProteinID, Sequence)), 
                by = "ProteinID")

effector_seq %>% 
  saveRDS("eff_seq_with_class.RDS")
```

### Function to filter the data based on the class

``` r
get_seq_each_class <- function(df_effector, class){
  
  df_seq <- df_effector %>% 
    dplyr::filter(class == class) %>%
    dplyr::select(c(ProteinID, Sequence))
  
  return(df_seq)
}
```

### Get effector data for each class

``` r
effector_seq_with_class <- readRDS("../../../data/secreted_data/data_processed_after_signalp/eff_seq_with_class.RDS")

effector_seq_fungi <- get_seq_each_class(effector_seq_with_class, "fungi")
effector_seq_bacteria <- get_seq_each_class(effector_seq_with_class, "bacteria")
effector_seq_oomycete <- get_seq_each_class(effector_seq_with_class, "oomycete")
```

### Get the non-effector candidate data from sampling the SignalP prediction data

``` r
non_effector_seq_fungi <- readRDS("../../../data/secreted_data/updated_signalp_results/fungi_non_eff_secreted_data.RDS") %>%
  dplyr::select(c(ID, sequence))
```

``` r
non_effector_seq_bacteria <- readRDS("../../../data/secreted_data/updated_signalp_results/bacteria_non_eff_secreted_data.RDS") %>% 
  dplyr::select(c(ID, sequence))
```

``` r
non_effector_seq_oomycete <- readRDS("../../../data/secreted_data/updated_signalp_results/oomycete_non_eff_secreted_data.RDS") %>% 
  dplyr::select(c(ID, sequence))
```

Dataframe to Fasta data
-----------------------

Since teh source data is R dataframe, we need to convert it to `.fasta`
data. And in order to do so, we can use the function `dataframe2fas()`
which is available in `library(seqRFLP)`.

### Function to convert dataframe to fasta data

``` r
# Function to create fasta data

get_fasta <- function(df) {
  df_name <- deparse(substitute(df))
  data_fa <- df %>%
    as.data.frame() %>% 
    seqRFLP::dataframe2fas(file = paste0("../../../data/secreted_data/data_processed_after_signalp/", df_name, ".fasta"))
}
```

### Get fasta data for each pathogen organism

``` r
# Getting the seq fasta 

# Fungi
get_fasta(effector_seq_fungi)
get_fasta(non_effector_seq_fungi)

# Bacteria
get_fasta(effector_seq_bacteria)
get_fasta(non_effector_seq_bacteria)

# Oomycete
get_fasta(effector_seq_oomycete)
get_fasta(non_effector_seq_oomycete)
```

Blast the data
--------------

### Function to blast data from R

``` r
get_blast_data <- function(database_fasta_path, query_fasta_path) {

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
  result_name <- paste0("../../../data/secreted_data/data_processed_after_signalp/", db_name, "_vs_", query_name, ".tsv")


  # Making the database
  system(paste("makeblastdb ", "-in ", database_fasta_path, "-dbtype ", "prot"), intern = TRUE)

  # Blast the database against the query name
  system(paste("blastp ", "-query ", query_fasta_path, "-db ", database_fasta_path, "-out ", result_name, " -outfmt ", "\"6  qseqid qlen sseqid slen length nident mismatch positive\""))
}
```

``` r
# Check teh current directory 

getwd()
```

    ## [1] "/Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/r-scripts/getting-secreted-data"

``` r
# Fungi 
get_blast_data("../../../data/secreted_data/data_processed_after_signalp/effector_seq_fungi.fasta", "../../../data/secreted_data/data_processed_after_signalp/non_effector_seq_fungi.fasta")

# Bacteria
get_blast_data("../../../data/secreted_data/data_processed_after_signalp/effector_seq_bacteria.fasta", "../../../data/secreted_data/data_processed_after_signalp/non_effector_seq_bacteria.fasta")

# Oomycete
get_blast_data("../../../data/secreted_data/data_processed_after_signalp/effector_seq_oomycete.fasta", "../../../data/secreted_data/data_processed_after_signalp/non_effector_seq_oomycete.fasta")
```

Read the BLAST result
---------------------

### Function to convert the result of the blast to dataframe

``` r
# function to read the results and get the list of row index
blast_results <- function(result_path) {

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
      percent_indentical = (nident / max(qlen, slen)) * 100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive / max(qlen, slen)) * 100 # The percentage of positive sequence over the longer sequence
    )

  # Get the data frame where the percent identical > 90
  df_identical_protein <- df_results %>%
    filter(percent_indentical > 90)

  # Get the row indices of the subject data for all of the identical percentage > 90%
  subject_index_list_to_remove <- df_results %>%
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
    df_identical_protein = df_identical_protein,
    subject_index = subject_index_list_to_remove,
    query_index = query_index_list_to_remove
  )

  # Return the lists
  return(list_results)
}
```

### Check the identical protein to remove

``` r
blast_results("../../../data/secreted_data/data_processed_after_signalp/effector_seq_fungi_vs_non_effector_seq_fungi.tsv")[2] %>% 
  knitr::kable()
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 20%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">qseqid</th>
<th style="text-align: right;">qlen</th>
<th style="text-align: left;">sseqid</th>
<th style="text-align: right;">slen</th>
<th style="text-align: right;">length</th>
<th style="text-align: right;">nident</th>
<th style="text-align: right;">mismatch</th>
<th style="text-align: right;">positive</th>
<th style="text-align: right;">percent_indentical</th>
<th style="text-align: right;">percent_positive</th>
</tr>
</thead>
<tbody>
</tbody>
</table>

</td>
</tr>
</tbody>
</table>
``` r
blast_results("../../../data/secreted_data/data_processed_after_signalp/effector_seq_bacteria_vs_non_effector_seq_bacteria.tsv")[2] %>% 
  knitr::kable()
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 20%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">qseqid</th>
<th style="text-align: right;">qlen</th>
<th style="text-align: left;">sseqid</th>
<th style="text-align: right;">slen</th>
<th style="text-align: right;">length</th>
<th style="text-align: right;">nident</th>
<th style="text-align: right;">mismatch</th>
<th style="text-align: right;">positive</th>
<th style="text-align: right;">percent_indentical</th>
<th style="text-align: right;">percent_positive</th>
</tr>
</thead>
<tbody>
</tbody>
</table>

</td>
</tr>
</tbody>
</table>
``` r
blast_results("../../../data/secreted_data/data_processed_after_signalp/effector_seq_oomycete_vs_non_effector_seq_oomycete.tsv")[2] %>% 
  knitr::kable()
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 20%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">qseqid</th>
<th style="text-align: right;">qlen</th>
<th style="text-align: left;">sseqid</th>
<th style="text-align: right;">slen</th>
<th style="text-align: right;">length</th>
<th style="text-align: right;">nident</th>
<th style="text-align: right;">mismatch</th>
<th style="text-align: right;">positive</th>
<th style="text-align: right;">percent_indentical</th>
<th style="text-align: right;">percent_positive</th>
</tr>
</thead>
<tbody>
</tbody>
</table>

</td>
</tr>
</tbody>
</table>
### Summary / Conclusion

Based on teh results above, as we can see that there is no identical
protein between effector and non-effector for each pathogen data
(bacteria, fungi, and oomycete). This can be used to make sure that we
do not randomly choose effector data in negative datasets. Also to
minimize the probablity to choose effector data in negative datasets.
