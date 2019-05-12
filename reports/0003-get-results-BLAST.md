Getting the identical protein using BLAST
=========================================

Introduction
------------

After getting all of data sets (training, validation / development, and
testing), we need to make sure that there is no identical protein
between them, otherwise it would not help the deep learning model to
learn the data. In order to identify the identical protein sequence
between the data sets, we can use BLAST.

Using BLAST
-----------

We can use terminal to run BLAST as follows

1.  Defining database

Now let us BLAST both validation or development datasets against the
entire training datasets. First, we need to tell BLAST that the protein
sequences in training datasets are the database. That’s done by calling
‘makeblastdb’ (creating database using local file):

``` bash
$ makeblastdb -in blast_train.fasta -dbtype prot
```

1.  Check the query to the source data base

Next, we call BLAST to do the search.

``` bash
$ blastp -query blast_val.fasta -db blast_train.fasta -out blast_train_x_val.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

``` bash
blastp -query blast_test.fasta -db blast_train.fasta -out blast_train_x_test.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

Here, `-outfmt` indicates how we specify the output format. The option
`6` will give us the output in tabular, with the supported format
specifiers:

-   qseqid means Query Seq-id
-   qlen means Query sequence length
-   sseqid means Subject Seq-id
-   slen means Subject sequence length
-   qstart means Start of alignment in query
-   length means Alignment length
-   nident means Number of identical matches
-   mismatch means Number of mismatches
-   positive means Number of positive-scoring matches

Results Analysis using R
------------------------

### Defining the function

``` r
# read the result of BLAST
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.1       ✔ purrr   0.3.2  
    ## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
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
      percent_indentical = (nident/max(qlen, slen))*100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive/max(qlen, slen))*100 # The percentage of positive sequence over the longer sequence
      )

  # Get the row indices of the source data for all of the identical percentage > 90%
  source_index_list_to_remove <- df_results %>%
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
    all_val = blast_results_train_x_val,
    all_test = blast_results_train_x_test,
    all = all_
  )

  return(results_list)
}
```

### Use the function defined

``` r
# define all the path of the results and the datasets in .csv
blast_train_x_val_path = "../../data/BLAST-data/blast_train_x_val.tsv"
blast_train_x_test_path = "../../data/BLAST-data/blast_train_x_test.tsv"
input_train_path = "../../data/BLAST-data/blast_train.csv"
label_train_path = "../../data/BLAST-data/blast_label_train.csv"
input_val_path = "../../data/BLAST-data/blast_val.csv"
label_val_path = "../../data/BLAST-data/blast_label_val.csv"
input_test_path = "../../data/BLAST-data/blast_test.csv"
label_test_path = "../../data/BLAST-data/blast_label_test.csv"


all_results <-  get_info(blast_train_x_val_path,
                                     blast_train_x_test_path,
                                     input_train_path,
                                     label_train_path,
                                     input_val_path,
                                     label_val_path,
                                     input_test_path,
                                     label_test_path)
```

``` r
df_val_x_train <- all_results[["all_val"]][["df"]]
 
df_val_x_train %>% 
  head(10) %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_indentical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|--------------------:|------------------:|
|       1|  1027|     154|  2037|      39|      15|        24|        23|            0.7363770|          1.1291114|
|       1|  1027|      32|  2210|      55|      21|        33|        30|            0.9502262|          1.3574661|
|       1|  1027|      32|  2210|      32|      14|        17|        20|            0.6334842|          0.9049774|
|       1|  1027|       5|  1047|      53|      16|        34|        26|            1.5281757|          2.4832856|
|       1|  1027|     361|  2212|      55|      21|        33|        29|            0.9493671|          1.3110307|
|       1|  1027|     361|  2212|      32|      14|        17|        20|            0.6329114|          0.9041591|
|       1|  1027|     115|  2266|      94|      24|        68|        40|            1.0591350|          1.7652251|
|       1|  1027|     401|  1548|      92|      29|        54|        43|            1.8733850|          2.7777778|
|       1|  1027|     116|  1635|      56|      20|        33|        28|            1.2232416|          1.7125382|
|       1|  1027|      58|  2469|      94|      23|        69|        39|            0.9315512|          1.5795869|

``` r
df_test_x_train <- all_results[["all_test"]][["df"]]

df_test_x_train %>% 
  head(10) %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_indentical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|--------------------:|------------------:|
|       1|  1926|     226|  1926|    1926|    1926|         0|      1926|          100.0000000|        100.0000000|
|       1|  1926|     189|   266|      47|      16|        30|        24|            0.8307373|          1.2461059|
|       1|  1926|     335|  1838|      34|      15|        18|        20|            0.7788162|          1.0384216|
|       1|  1926|     372|  4034|      94|      27|        51|        41|            0.6693109|          1.0163609|
|       1|  1926|     360|  1945|      28|       9|        19|        13|            0.4627249|          0.6683805|
|       2|   331|     523|   352|     320|     234|        81|       262|           66.4772727|         74.4318182|
|       2|   331|     370|   352|     320|     234|        81|       262|           66.4772727|         74.4318182|
|       2|   331|     453|   321|     323|     147|       158|       197|           44.4108761|         59.5166163|
|       2|   331|     546|   323|     318|     123|       165|       168|           37.1601208|         50.7552870|
|       2|   331|      66|  2124|      99|      28|        57|        45|            1.3182674|          2.1186441|

Lastly, we can get the result of the summary of total number of sequence
that has to be removed from each data sets:

``` r
summary <- all_results[["all"]]

summary %>% 
  knitr::kable()
```

| label       |  testing|  training|  validation|
|:------------|--------:|---------:|-----------:|
| effector    |       11|        21|           6|
| noneffector |       28|        76|          36|
