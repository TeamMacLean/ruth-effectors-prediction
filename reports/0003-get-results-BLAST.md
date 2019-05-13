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

-   Defining database

Now let us BLAST both validation or development datasets against the
entire training datasets. First, we need to tell BLAST that the protein
sequences in training datasets are the database. That’s done by calling
‘makeblastdb’ (creating database using local file):

``` bash
$ makeblastdb -in blast_train.fasta -dbtype prot
```

-   Check the query to the source data base

Next, we call BLAST to do the search.

``` bash
$ blastp -query blast_val.fasta -db blast_train.fasta -out blast_train_x_val.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

``` bash
$ blastp -query blast_test.fasta -db blast_train.fasta -out blast_train_x_test.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
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

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.1       ✔ purrr   0.3.2  
    ## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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
  
  # Get the data frame where the percent identical > 90
  df_identical_protein <- df_results %>% 
    filter(percent_indentical > 90)

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
    df_identical_protein = df_identical_protein,
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
    removed_training = results_after_removed_training, 
    removed_validation = results_after_removed_validation, 
    removed_testing = results_after_removed_testing,
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

#### Results for comparing the validation data set and training dataset

Here is teh result of the data that need to be removed since they are
considered as identical (90%):

``` r
df_iden_train_x_val <- all_results[["all_val"]][["df_identical_protein"]]

df_iden_train_x_val %>% 
  select(-c(mismatch, positive, percent_positive)) %>%
  rename("Query Seq-ID" = qseqid,
         "Query Length" = qlen,
         "Subject Seq-ID" = sseqid,
         "Subject Length" = slen,
         "Alignment Length" = length,
         "No Identical Matches" = nident,
         "% Identical Matches" =  percent_indentical
          ) %>%
  knitr::kable()
```

|  Query Seq-ID|  Query Length|  Subject Seq-ID|  Subject Length|  Alignment Length|  No Identical Matches|  % Identical Matches|
|-------------:|-------------:|---------------:|---------------:|-----------------:|---------------------:|--------------------:|
|             2|          2022|             109|            2015|              2024|                  2005|             99.15925|
|             2|          2022|             412|            2001|              2022|                  1998|             98.81306|
|             2|          2022|             184|            2023|              2047|                  1965|             97.13297|
|             3|          2057|             528|            2057|              2057|                  2048|             99.56247|
|             3|          2057|             106|            2035|              2035|                  2034|             98.88187|
|             6|          2290|             548|            2301|              2301|                  2097|             91.13429|
|            16|          1548|             222|            1548|              1546|                  1438|             92.89406|
|            26|          2048|             273|            2047|              2048|                  2040|             99.60938|
|            26|          2048|              84|            2048|              2048|                  2041|             99.65820|
|            27|          2339|              51|            2339|              2339|                  2339|            100.00000|
|            27|          2339|             542|            2401|              2401|                  2177|             90.67055|
|            29|          2046|             530|            2049|              2049|                  2042|             99.65837|
|            29|          2046|             332|            2050|              2050|                  2038|             99.41463|
|            29|          2046|             480|            2053|              2053|                  2039|             99.31807|
|            29|          2046|             249|            2053|              2053|                  2038|             99.26936|
|            29|          2046|             117|            2055|              2055|                  2037|             99.12409|
|            29|          2046|              31|            2067|              2067|                  2034|             98.40348|
|            33|          1140|             180|            1140|              1140|                  1140|            100.00000|
|            33|          1140|             385|            1141|              1140|                  1098|             96.23138|
|            34|          2419|             145|            2421|              2419|                  2188|             90.37588|
|            37|          2053|             480|            2053|              2053|                  2047|             99.70774|
|            37|          2053|             249|            2053|              2053|                  2044|             99.56162|
|            37|          2053|             117|            2055|              2055|                  2041|             99.31873|
|            37|          2053|             332|            2050|              2053|                  2039|             99.31807|
|            37|          2053|             530|            2049|              2053|                  2036|             99.17194|
|            37|          2053|              31|            2067|              2067|                  2038|             98.59700|
|            40|          2043|             390|            2043|              2043|                  2043|            100.00000|
|            40|          2043|             505|            2037|              2044|                  2009|             98.33578|
|            40|          2043|             540|            2039|              2045|                  1991|             97.45472|
|            42|          2339|              34|            2314|              2314|                  2287|             97.77683|
|            42|          2339|             234|            2303|              2303|                  2266|             96.87901|
|            43|          2003|             275|            2003|              2003|                  2003|            100.00000|
|            53|          2053|             249|            2053|              2053|                  2049|             99.80516|
|            53|          2053|             480|            2053|              2053|                  2042|             99.46420|
|            53|          2053|             332|            2050|              2053|                  2042|             99.46420|
|            53|          2053|             117|            2055|              2055|                  2040|             99.27007|
|            53|          2053|             530|            2049|              2053|                  2039|             99.31807|
|            53|          2053|              31|            2067|              2067|                  2039|             98.64538|
|            54|           561|             268|             561|               561|                   561|            100.00000|
|            70|          2393|             331|            2393|              2393|                  2384|             99.62390|
|            75|          2073|             520|            2073|              2072|                  1969|             94.98312|
|            75|          2073|             338|            2073|              2072|                  1971|             95.07959|
|            75|          2073|             573|            2073|              2072|                  1989|             95.94790|
|            77|          2276|              44|            2276|              2276|                  2276|            100.00000|
|            82|          1095|               3|            1096|              1096|                   994|             90.69343|
|            91|          2393|             331|            2393|              2393|                  2384|             99.62390|
|            93|          2154|             100|            2154|              2154|                  2085|             96.79666|
|            93|          2154|             420|            2154|              2154|                  1948|             90.43640|
|            93|          2154|             285|            2152|              2154|                  1944|             90.25070|
|            97|           352|             523|             352|               352|                   350|             99.43182|
|            97|           352|             370|             352|               352|                   350|             99.43182|
|            98|          2056|             184|            2023|              2066|                  1981|             96.35214|
|            98|          2056|             412|            2001|              2057|                  1963|             95.47665|
|           103|          2455|             465|            2475|              2454|                  2450|             98.98990|
|           103|          2455|             177|            2386|              2384|                  2304|             93.84929|
|           107|           311|             352|             311|               311|                   308|             99.03537|
|           107|           311|             443|             311|               311|                   306|             98.39228|
|           108|          2053|             480|            2053|              2053|                  2042|             99.46420|
|           108|          2053|             249|            2053|              2053|                  2039|             99.31807|
|           108|          2053|             530|            2049|              2053|                  2041|             99.41549|
|           108|          2053|              31|            2067|              2067|                  2043|             98.83890|
|           108|          2053|             117|            2055|              2055|                  2038|             99.17275|
|           108|          2053|             332|            2050|              2053|                  2038|             99.26936|
|           109|          2049|             530|            2049|              2049|                  2042|             99.65837|
|           109|          2049|             249|            2053|              2053|                  2040|             99.36678|
|           109|          2049|             480|            2053|              2053|                  2039|             99.31807|
|           109|          2049|             117|            2055|              2055|                  2039|             99.22141|
|           109|          2049|             332|            2050|              2053|                  2039|             99.46341|
|           109|          2049|              31|            2067|              2067|                  2040|             98.69376|
|           117|          2035|             423|            2041|              2035|                  1991|             97.55022|
|           117|          2035|             508|            2035|              2035|                  1991|             97.83784|
|           117|          2035|             367|            2041|              2035|                  1982|             97.10926|
|           124|          2260|             340|            2260|              2260|                  2256|             99.82301|
|           124|          2260|             417|            2260|              2260|                  2231|             98.71681|
|           125|          2310|             151|            2412|              2312|                  2284|             94.69320|
|           125|          2310|             564|            2395|              2310|                  2278|             95.11482|
|           125|          2310|             247|            2237|              2236|                  2179|             94.32900|
|           130|           163|             208|             163|               163|                   162|             99.38650|
|           132|          2006|             437|            2006|              2006|                  1995|             99.45165|
|           132|          2006|             427|            2006|              2006|                  1995|             99.45165|
|           132|          2006|             353|            2006|              2006|                  1995|             99.45165|
|           132|          2006|             221|            2006|              2006|                  1986|             99.00299|
|           132|          2006|             455|            2003|              2006|                  1987|             99.05284|
|           132|          2006|             192|            2006|              2005|                  1980|             98.70389|
|           132|          2006|             307|            2006|              2006|                  1981|             98.75374|
|           134|          1752|             330|            1752|              1752|                  1752|            100.00000|
|           135|          2290|              92|            2290|              2290|                  2286|             99.82533|
|           135|          2290|             125|            2437|              2290|                  2254|             92.49077|
|           154|          2493|             465|            2475|              2475|                  2436|             97.71360|
|           154|          2493|             177|            2386|              2384|                  2292|             91.93742|
|           162|          2276|              44|            2276|              2276|                  2275|             99.95606|
|           170|          2046|             117|            2055|              2055|                  2037|             99.12409|
|           170|          2046|             249|            2053|              2053|                  2034|             99.07453|
|           170|          2046|             480|            2053|              2053|                  2035|             99.12323|
|           170|          2046|             332|            2050|              2053|                  2035|             99.26829|
|           170|          2046|              31|            2067|              2067|                  2034|             98.40348|
|           170|          2046|             530|            2049|              2053|                  2030|             99.07272|
|           171|          2442|             191|            2425|              2442|                  2332|             95.49550|
|           177|          2152|             100|            2154|              2154|                  2024|             93.96472|
|           178|          2224|             361|            2212|              2212|                  2212|             99.46043|
|           178|          2224|              32|            2210|              2193|                  2053|             92.31115|
|           182|          1028|               3|            1096|              1096|                  1005|             91.69708|
|           192|          2057|             528|            2057|              2057|                  2030|             98.68741|
|           192|          2057|             106|            2035|              2035|                  2008|             97.61789|

And the following table shows the identical protein between training
data sets and testing datasets.

``` r
df_iden_train_x_test <- all_results[["all_test"]][["df_identical_protein"]]

df_iden_train_x_test %>% 
  select(-c(mismatch, positive, percent_positive)) %>% 
  rename("Query Seq-ID" = qseqid,
         "Query Length" = qlen,
         "Subject Seq-ID" = sseqid,
         "Subject Length" = slen,
         "Alignment Length" = length,
         "No Identical Matches" = nident,
         "% Identical Matches" =  percent_indentical
          ) %>% 
  knitr::kable()
```

|  Query Seq-ID|  Query Length|  Subject Seq-ID|  Subject Length|  Alignment Length|  No Identical Matches|  % Identical Matches|
|-------------:|-------------:|---------------:|---------------:|-----------------:|---------------------:|--------------------:|
|             1|          1926|             226|            1926|              1926|                  1926|            100.00000|
|             4|          1126|              81|            1163|              1163|                  1108|             95.27085|
|             4|          1126|             346|            1164|              1164|                  1094|             93.98625|
|             4|          1126|               3|            1096|              1130|                  1063|             94.40497|
|            14|          2059|             350|            2059|              2059|                  2059|            100.00000|
|            15|          2260|             340|            2260|              2260|                  2256|             99.82301|
|            15|          2260|             417|            2260|              2260|                  2231|             98.71681|
|            16|           125|               8|             126|               124|                   115|             91.26984|
|            26|          2273|              60|            2273|              2282|                  2255|             99.20810|
|            26|          2273|              96|            2323|              2332|                  2253|             96.98666|
|            31|          2210|              32|            2210|              2210|                  2210|            100.00000|
|            31|          2210|             361|            2212|              2212|                  2061|             93.17360|
|            34|          2257|             151|            2412|              2259|                  2221|             92.08126|
|            34|          2257|             564|            2395|              2257|                  2213|             92.40084|
|            34|          2257|             247|            2237|              2236|                  2171|             96.18963|
|            37|          2326|              25|            2345|              2328|                  2275|             97.01493|
|            37|          2326|             495|            2332|              2332|                  2274|             97.51286|
|            41|          2031|              85|            2005|              2005|                  2004|             98.67061|
|            44|            85|             197|              85|                85|                    83|             97.64706|
|            47|          2255|             148|            2255|              2255|                  2231|             98.93570|
|            48|          2014|             476|            2014|              2018|                  1821|             90.41708|
|            49|          2154|             100|            2154|              2154|                  2022|             93.87187|
|            49|          2154|             420|            2154|              2154|                  1942|             90.15785|
|            49|          2154|             285|            2152|              2154|                  1941|             90.11142|
|            53|          1548|             222|            1548|              1548|                  1545|             99.80620|
|            55|          1096|               3|            1096|              1096|                  1076|             98.17518|
|            57|          2406|             548|            2301|              2301|                  2265|             94.13965|
|            59|          2004|             159|            2004|              2004|                  1931|             96.35729|
|            63|           543|             471|             543|               543|                   512|             94.29098|
|            64|          2015|             109|            2015|              2019|                  1979|             98.21340|
|            64|          2015|             412|            2001|              2015|                  1975|             98.01489|
|            64|          2015|             184|            2023|              2030|                  1971|             97.42956|
|            68|          2380|             559|            2380|              2380|                  2342|             98.40336|
|            68|          2380|             311|            2380|              2380|                  2342|             98.40336|
|            68|          2380|             213|            2388|              2388|                  2239|             93.76047|
|            78|          2301|             548|            2301|              2301|                  2266|             98.47892|
|            83|           720|              72|             720|               720|                   694|             96.38889|
|            89|          2256|             151|            2412|              2258|                  2220|             92.03980|
|            89|          2256|             564|            2395|              2256|                  2212|             92.35908|
|            89|          2256|             247|            2237|              2236|                  2171|             96.23227|
|            92|          2255|             148|            2255|              2255|                  2234|             99.06874|
|            93|          2066|             187|            2066|              2066|                  2042|             98.83833|
|           106|          1858|             103|            1858|              1858|                  1858|            100.00000|
|           106|          1858|              53|            1859|              1859|                  1850|             99.51587|
|           106|          1858|             570|            1852|              1856|                  1830|             98.49300|
|           106|          1858|             507|            1852|              1856|                  1830|             98.49300|
|           108|           138|             252|             138|               138|                   138|            100.00000|
|           112|           147|             371|             147|               147|                   144|             97.95918|
|           121|          2006|             437|            2006|              2006|                  2006|            100.00000|
|           121|          2006|             427|            2006|              2006|                  2006|            100.00000|
|           121|          2006|             353|            2006|              2006|                  2006|            100.00000|
|           121|          2006|             455|            2003|              2006|                  1990|             99.20239|
|           121|          2006|             221|            2006|              2006|                  1981|             98.75374|
|           121|          2006|             307|            2006|              2006|                  1980|             98.70389|
|           121|          2006|             192|            2006|              2005|                  1975|             98.45464|
|           134|           240|             334|             240|               240|                   239|             99.58333|
|           135|          2405|             542|            2401|              2405|                  2253|             93.67983|
|           135|          2405|              51|            2339|              2409|                  2193|             91.18503|
|           156|          2049|              85|            2005|              2005|                  2002|             97.70620|
|           171|          2124|             456|            2018|              2121|                  2011|             94.67985|
|           172|           685|             279|             685|               685|                   684|             99.85401|
|           180|          2049|              85|            2005|              2005|                  2004|             97.80381|
|           181|           223|             127|             224|               224|                   218|             97.32143|
|           181|           223|              49|             224|               224|                   218|             97.32143|
|           181|           223|             509|             224|               224|                   217|             96.87500|
|           181|           223|             301|             224|               224|                   215|             95.98214|
|           185|          2179|              34|            2314|              2179|                  2154|             93.08557|
|           185|          2179|             234|            2303|              2179|                  2147|             93.22623|
|           191|          1548|             222|            1548|              1548|                  1544|             99.74160|

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
