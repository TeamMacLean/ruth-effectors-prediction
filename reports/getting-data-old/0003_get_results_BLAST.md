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

-   Check the query to the subject data base

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
      percent_identical = (nident/max(qlen, slen))*100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive/max(qlen, slen))*100 # The percentage of positive sequence over the longer sequence
      )
  
  # Get the data frame where the percent identical > 90
  df_identical_protein <- df_results %>% 
    filter(percent_identical > 90)

  # Get the row indices of the subject data for all of the identical percentage > 90%
  subject_index_list_to_remove <- df_results %>%
    filter(percent_identical > 90) %>%
    select(sseqid) %>%
    unique() %>%
    unlist()

  # Get the row indices of the query data for all of the identical percentage > 90%
  query_index_list_to_remove <- df_results %>%
    filter(percent_identical > 90) %>%
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
    dplyr::filter(!(row_number() %in% drop))
  
  # Combine the input data and the label, then take the data that should be removed
  removed_rows <- df_input %>%
    cbind(df_label) %>%
    dplyr::filter(row_number() %in% drop)

  # Get the information about the removed columns
  removed_rows_freq <- df_label %>%
    filter(row_number() %in% drop) %>%
    select(label) %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("label", "freq")) %>%
    mutate(label = ifelse(label == 1, "effector", "noneffector"))

  # Create list for the results
  results_list <- list(
    df = df_new,
    removed_freq = removed_rows_freq, 
    removed_rows = removed_rows
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
  subject_index_train_x_val <- blast_results_train_x_val[["subject_index"]]
  query_index_train_x_val <- blast_results_train_x_val[["query_index"]]

  # Getting the results from comparing testing dataset and training dataset
  blast_results_train_x_test <- blast_results(blast_train_x_test_path)
  df_train_x_test <- blast_results_train_x_test[["df"]]
  subject_index_train_x_test <- blast_results_train_x_test[["subject_index"]]
  query_index_train_x_test <- blast_results_train_x_test[["query_index"]]

  # Remove all of the rows of the training data
  # Since there might be intersection between index on the training when comparing them with validation and testing sets,
  # we need to find the intesections

  intersec_training_index <- c(subject_index_train_x_val, subject_index_train_x_test) %>%
    unique() %>%
    unlist()

  results_after_removed_training <- remove_from_data_sets(input_train_path, label_train_path, intersec_training_index)
  removed_rows_training_freq <- results_after_removed_training[["removed_freq"]] %>%
    mutate(type = "training")
  removed_rows_training <- results_after_removed_training[["removed_rows"]]

  results_after_removed_validation <- remove_from_data_sets(input_val_path, label_val_path, query_index_train_x_val)
  removed_rows_val_freq <- results_after_removed_validation[["removed_freq"]] %>%
    mutate(type = "validation")
  removed_rows_val <- results_after_removed_validation[["removed_rows"]]
  
  
  results_after_removed_testing <- remove_from_data_sets(input_test_path, label_test_path, query_index_train_x_test)
  removed_rows_test_freq <- results_after_removed_testing[["removed_freq"]] %>%
    mutate(type = "testing")
  removed_rows_test <- results_after_removed_testing[["removed_rows"]]

  all_removed_rows_freq <- removed_rows_training_freq %>%
    rbind(., removed_rows_val_freq) %>%
    rbind(., removed_rows_test_freq)
  
  all_removed_rows <- removed_rows_training %>%
    rbind(., removed_rows_val) %>%
    rbind(., removed_rows_test)
  
  all_ <- all_removed_rows_freq %>%
    reshape2::melt(id.var = "label")

  # all_ <- reshape2::dcast(all_removed_rows_freq, label ~ type, value.var = "freq") %>%
    # as.data.frame()

  results_list <- list(
    all_val = blast_results_train_x_val,
    all_test = blast_results_train_x_test,
    removed_training = results_after_removed_training, 
    removed_validation = results_after_removed_validation, 
    removed_testing = results_after_removed_testing,
    all_removed_rows_freq = all_removed_rows_freq, 
    all_removed_rows = all_removed_rows,
    all = all_
  )

  return(results_list)
}
```

### Use the function defined

``` r
# define all the path of the results and the datasets in .csv
blast_train_x_val_path = "../../../data/getting-data-old/BLAST-data/blast_train_x_val.tsv"
blast_train_x_test_path = "../../../data/getting-data-old/BLAST-data/blast_train_x_test.tsv"
input_train_path = "../../../data/getting-data-old/BLAST-data/blast_train.csv"
label_train_path = "../../../data/getting-data-old/BLAST-data/blast_label_train.csv"
input_val_path = "../../../data/getting-data-old/BLAST-data/blast_val.csv"
label_val_path = "../../../data/getting-data-old/BLAST-data/blast_label_val.csv"
input_test_path = "../../../data/getting-data-old/BLAST-data/blast_test.csv"
label_test_path = "../../../data/getting-data-old/BLAST-data/blast_label_test.csv"


all_results <-  get_info(blast_train_x_val_path,
                         blast_train_x_test_path,
                         input_train_path,
                         label_train_path,
                         input_val_path,
                         label_val_path,
                         input_test_path,
                         label_test_path)

all_removed_rowss_freq <- all_results[["all_removed_rows_freq"]]

# write.csv(all_removed_rowss_freq, "../../../data/getting-data-old/all_removed_rows_freq.csv", row.names = FALSE)
```

#### Results for comparing the validation data set and training dataset

Here is the result of the data that need to be removed since they are
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
         "% Identical Matches" =  percent_identical
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

Below result is one of the example for pairwise for each sequence:

``` bash
Query= 2

Length=2022
                                                                      Score        E
Sequences producing significant alignments:                          (Bits)     Value

109                                                                   3934       0.0   
412                                                                   3915       0.0   
184                                                                   3823       0.0   
381                                                                   579        2e-173
221                                                                   576        4e-173
192                                                                   574        3e-172
455                                                                   570        5e-171
307                                                                   570        1e-170
437                                                                   569        2e-170
427                                                                   569        2e-170
353                                                                   569        2e-170
528                                                                   566        2e-169
106                                                                   565        5e-169
55                                                                    550        1e-163
520                                                                   540        4e-160
338                                                                   534        6e-158
573                                                                   524        1e-154
508                                                                   512        1e-150
423                                                                   512        1e-150
367                                                                   510        4e-150
524                                                                   478        2e-139
56                                                                    478        2e-139
420                                                                   262        2e-71 
285                                                                   262        3e-71 
38                                                                    261        4e-71 
100                                                                   249        2e-67 
253                                                                   239        4e-64 
58                                                                    234        9e-63 
115                                                                   232        3e-62 
116                                                                   231        6e-62 
235                                                                   230        1e-61 
47                                                                    230        1e-61 
175                                                                   230        1e-61 
568                                                                   225        5e-60 
120                                                                   223        2e-59 
364                                                                   221        1e-58 
44                                                                    202        3e-53 
284                                                                   190        3e-49 
241                                                                   188        9e-49 
548                                                                   187        2e-48 
494                                                                   186        5e-48 
533                                                                   178        8e-46 
398                                                                   178        9e-46 
512                                                                   157        4e-39 
372                                                                   151        1e-37 
532                                                                   146        5e-36 
478                                                                   136        4e-33 
213                                                                   117        3e-27 
559                                                                   110        4e-25 
311                                                                   110        4e-25 
417                                                                   70.5       6e-13 
340                                                                   68.6       3e-12 
193                                                                   55.8       2e-08 
451                                                                   48.9       2e-06 
43                                                                    39.7       0.002 
541                                                                   33.9       0.073 
444                                                                   32.7       0.15  
176                                                                   31.2       0.47  
521                                                                   29.6       1.5   
187                                                                   29.3       1.8   
431                                                                   27.3       6.9   
85                                                                    26.9       9.3   


>109
Length=2015

 Score = 3934 bits (10202),  Expect = 0.0, Method: Compositional matrix adjust.
 Identities = 2005/2024 (99%), Positives = 2005/2024 (99%), Gaps = 11/2024 (1%)

Query  1     MSFEDVLETCRERQIELWCDGEQLRYRAPAGALDAPLAARINAQRAAFVRFLRDGAWRAQ  60
             MSFEDVLETCRERQIELWCDGEQLRYRAPAGALDAPLAARINAQRAAFVRFLRDGAWRAQ
Sbjct  1     MSFEDVLETCRERQIELWCDGEQLRYRAPAGALDAPLAARINAQRAAFVRFLRDGAWRAQ  60

Query  61    PSRTHERFALTPVQAAYVLGRHAAFEFGGSACHLYVEYRESADLDVRRFEAAWNACVARH  120
             PSRTHERFALTPVQAAYVLGRHAAFEFGGSACHLYVEYRESADLDVRRFEAAWNACVARH
Sbjct  61    PSRTHERFALTPVQAAYVLGRHAAFEFGGSACHLYVEYRESADLDVRRFEAAWNACVARH  120

Query  121   PMLRAIVEDNAWQRILPDVPWQTLAVHDLRDANAAAFDAHVRRVRERLDHAVHALDRWPV  180
             PMLRAIVEDNAWQRILPDVPWQTLAVHDLRDANAAAFDAHVRRVRERLDHAVHALDRWPV
Sbjct  121   PMLRAIVEDNAWQRILPDVPWQTLAVHDLRDANAAAFDAHVRRVRERLDHAVHALDRWPV  180

Query  181   LQPEVSIGPDGAVLHLSVDFTLIDYASLQLLLAEWRRRHDDPQWKPAPLDVTFRDYVVNE  240
             LQPEVSIGPDGAVLHLSVDFTLIDYASLQLLLAEWRRRHDDPQWKPAPLDVTFRDYVVNE
Sbjct  181   LQPEVSIGPDGAVLHLSVDFTLIDYASLQLLLAEWRRRHDDPQWKPAPLDVTFRDYVVNE  240

Query  241   ARAGRRAQHARDAAWWLARIDGLPGRPDLPVLPRPPEGRADARPRFTHRHARLDRARWDA  300
             ARAGRRAQHARDAAWWLARIDGLPGRPDLPVLPRPPEGRADARPRFTHRHARLDRARWD 
Sbjct  241   ARAGRRAQHARDAAWWLARIDGLPGRPDLPVLPRPPEGRADARPRFTHRHARLDRARWDV  300

Query  301   LVAFASRFGLSPAGVALAAFAEVVGRWSQSPAFCLNLTVLNRPPVHPQIDAVLGDFTALS  360
             LVAFASRFGLSPAGVALAAFAEVVGRWSQSPAFCLNLTVLNRPPVHPQIDAVLGDFTALS
Sbjct  301   LVAFASRFGLSPAGVALAAFAEVVGRWSQSPAFCLNLTVLNRPPVHPQIDAVLGDFTALS  360

Query  361   LLAVDVAQGRDFAERARAIGAQMFDDLDHRAFTGVEVMRELARRRGKDAALMPVVFTSGI  420
             LLAVDVAQG DFA RARAIGAQMFDDLDHRAFTGVEVMRELARRRGKDAALMPVVFTSGI
Sbjct  361   LLAVDVAQGSDFAARARAIGAQMFDDLDHRAFTGVEVMRELARRRGKDAALMPVVFTSGI  420

Query  421   GSVGRLLGEHGEHGARAAPPCYMISQTPQVWLDCQVTDQFGGLEIGWDVRDDLFPAGMPA  480
             GSVGRLLGEHGEHGARAAPPCYMISQTPQVWLDCQVTDQFGGLEIGWDVRDDLFPAGMPA
Sbjct  421   GSVGRLLGEHGEHGARAAPPCYMISQTPQVWLDCQVTDQFGGLEIGWDVRDDLFPAGMPA  480

Query  481   AMFDAYVALLARLASDAAWWARPGDIVLPSGPPVACAHPDADRHLAAGFAAQARRTPDAT  540
             AMFDAYVALLARLASDAAWWARPGDIVLPSGPPVACAHPDADRHLAAGFAAQARRTPDAT
Sbjct  481   AMFDAYVALLARLASDAAWWARPGDIVLPSGPPVACAHPDADRHLAAGFAAQARRTPDAT  540

Query  541   VVIDAAGAHTYRDVAQRAAAVRAALERAGVAPGDKVAVRMPKGANQLVAVLGIVQAGAAY  600
             VVIDAAGAHTYRDVAQRAAAVRAALERAGVAPGDKVAVRMPKGANQLVAVLGIVQAGAAY
Sbjct  541   VVIDAAGAHTYRDVAQRAAAVRAALERAGVAPGDKVAVRMPKGANQLVAVLGIVQAGAAY  600

Query  601   VPIDYRQPALRRRAILRNAQVGAIVTERALDGEPDGCARIDVDALAPDPRWPPRDAHPLE  660
             VPIDYRQPALRRRAILRNAQVGAIVTERALDGEPDGCARIDVDALAPDPRWPPRDAHPLE
Sbjct  601   VPIDYRQPALRRRAILRNAQVGAIVTERALDGEPDGCARIDVDALAPDPRWPPRDAHPLE  660

Query  661   GDALAYVIYTSGSTGEPKGVMVSHAAVCNTLADINARHAVGAGDAVLGLAELSFDLSVYD  720
             GDALAYVIYTSGSTGEPKGVMVSHAAVCNTLADINARHAVGAGDAVLGLAELSFDLSVYD
Sbjct  661   GDALAYVIYTSGSTGEPKGVMVSHAAVCNTLADINARHAVGAGDAVLGLAELSFDLSVYD  720

Query  721   LFGATARGARVVLPDPARGNDPSHWAELIARHGVTLWNSVPAQGQMLIDYLESEPARAMP  780
             LFGATARGARVVLPDPARGNDPSHWAELIARHGVTLWNSVPAQGQMLIDYLESEPARAMP
Sbjct  721   LFGATARGARVVLPDPARGNDPSHWAELIARHGVTLWNSVPAQGQMLIDYLESEPARAMP  780

Query  781   GPRCVMWSGDWIPVSLPTRWWRRWPDSRLFSLGGATEASIWSIEHPIRPEDTRLVSIPYG  840
             GPRCVMWSGDWIPVSLPTRWWRRWPDSRLFSLGGATEASIWSIEHPIRPEDTRLVSIPYG
Sbjct  781   GPRCVMWSGDWIPVSLPTRWWRRWPDSRLFSLGGATEASIWSIEHPIRPEDTRLVSIPYG  840

Query  841   RALTGQTVDVLDALGRPCPPGVRGEIHIGGVGLATGYANDPARTAERFVRHADGRRLYRT  900
             RALTGQTVDVLDALGRPCPPGVRGEIHIGGVGLATGYANDPARTAERFVRHADGRRLYRT
Sbjct  841   RALTGQTVDVLDALGRPCPPGVRGEIHIGGVGLATGYANDPARTAERFVRHADGRRLYRT  900

Query  901   GDLGRRRADGSLEFLGRQDDQVKIRGYRIELAEIDAALSAHPRVASAATIVLGDAAQRRL  960
             GDLGRRRADGSLEFLGRQDDQVKIRGYRIELAEIDAALSAHPRVASAATIVLGDAAQRRL
Sbjct  901   GDLGRRRADGSLEFLGRQDDQVKIRGYRIELAEIDAALSAHPRVASAATIVLGDAAQRRL  960

Query  961   ASFVTLHGAAPDPRRRDAQLLAIAQRVRDALAAERWPARAEIRRSVGQLDAACVASLASW  1020
             ASFVTLHGAAPDPRRRDAQL AIAQRVRDALAAERWPARAEIRRSVGQLDAACVASLASW
Sbjct  961   ASFVTLHGAAPDPRRRDAQLHAIAQRVRDALAAERWPARAEIRRSVGQLDAACVASLASW  1020

Query  1021  IVRAGGLTRDASVDFATLAERLRVPAARHGLLRHWLALLEARGVLCTDTAGPVAGVDVTE  1080
             IVRAGGLTRDASVDFATLAERLRVPAARHGLLRHWLALLEARGVLCTDTAGPVAGVDVTE
Sbjct  1021  IVRAGGLTRDASVDFATLAERLRVPAARHGLLRHWLALLEARGVLCTDTAGPVAGVDVTE  1080

Query  1081  SENSAASRGAADADAFEAFEAFEAFEAFEAFEAFEAFEAFEAFEAFEAPEAPEAPEAPEA  1140
             SENSAASRGAADAD      AFEAFEAFEAFEAFEAFE FEAFEA EAPEAPEAPEAPEA
Sbjct  1081  SENSAASRGAADAD------AFEAFEAFEAFEAFEAFEVFEAFEAPEAPEAPEAPEAPEA  1134

Query  1141  PEAPEAPEAPEAPEAPEAPEATEATEATEATEAPEATEATEATEAIETIETIETIETIET  1200
             PEAPEAPEAPEAPEAPEAPEATEATEATEATEAPEA EATEATEA   IETIETIETIET
Sbjct  1135  PEAPEAPEAPEAPEAPEAPEATEATEATEATEAPEAPEATEATEA---IETIETIETIET  1191

Query  1201  IETIETIETIETPNAADVTKATKSAGTTDTTTSPAGWRVARAHAARATDAAGAVPLSDAD  1260
             IETIETIETIETPNAADVTKATKSAGTTDTTTSPAGWRVARAHAARATDAAGAVPLSDAD
Sbjct  1192  IETIETIETIETPNAADVTKATKSAGTTDTTTSPAGWRVARAHAARATDAAGAVPLSDAD  1251

Query  1261  ACWDAFARDASPALWPAVLIDYFRASAACLGEQVDDRVSPAGLMFPQGAPHVAEAMYSDG  1320
             ACWDAFARDASPALWPAVLIDYFRASAACLGEQVDDRVSPAGLMFPQGAPHVAEAMYSDG
Sbjct  1252  ACWDAFARDASPALWPAVLIDYFRASAACLGEQVDDRVSPAGLMFPQGAPHVAEAMYSDG  1311

Query  1321  MHARALHRAMAEAVSGIVAREPRRAWRILEIGAGTAAATRAIVDALAPHAESGVRIDYLF  1380
             MHARALHRAMAEAVSGIVAREPRRAWRILEIGAGTAAATRAIVDALAPHAESGVRIDYLF
Sbjct  1312  MHARALHRAMAEAVSGIVAREPRRAWRILEIGAGTAAATRAIVDALAPHAESGVRIDYLF  1371

Query  1381  SDVSSYFLAAARERFAAYPWVRFIRFDMNAPLDAQGVAPHSVDLIASSGALNNARDTVAL  1440
             SDVSSYFLAAARERFAAYPWVRFIRFDMNAPLDAQGVAPHSVDLIASSGALNNARDTVAL
Sbjct  1372  SDVSSYFLAAARERFAAYPWVRFIRFDMNAPLDAQGVAPHSVDLIASSGALNNARDTVAL  1431

Query  1441  VAGLRALSSADAWWVIQELTAEHPEISISQGLMMEPPRDARAASHRLFVHRAQWLEWLRA  1500
             VAGLRALSSADAWWVIQELTAEHPEISISQGLMMEPPRDARAASHRLFVHRAQWLEWLRA
Sbjct  1432  VAGLRALSSADAWWVIQELTAEHPEISISQGLMMEPPRDARAASHRLFVHRAQWLEWLRA  1491

Query  1501  AGDRALGCVEPDSPLDALGYDILLARVKADAARVDADELLAFAAERVPRYMVPSQLCALE  1560
             AGDRALGCVEPDSPLDALGYDILLARVKADAARVDADELLAFAAERVPRYMVPSQLCALE
Sbjct  1492  AGDRALGCVEPDSPLDALGYDILLARVKADAARVDADELLAFAAERVPRYMVPSQLCALE  1551

Query  1561  RLPVTANGKIDRRALAQIAGVHEPRAAHGRAVRPCDAEPLLARLIGLWEAVLDTQGVAAD  1620
             RLPVTANGKIDRRALAQIAGVHEPRAAHGRAVR CDAEPLLARLIGLWEAVLDTQGVAAD
Sbjct  1552  RLPVTANGKIDRRALAQIAGVHEPRAAHGRAVRACDAEPLLARLIGLWEAVLDTQGVAAD  1611

Query  1621  QDFFAAGGDSLLIAQLVSRVRAEEPLARAHPFDRLLRWALAQPTPAGLAQCLRDAKAQAQ  1680
             QDFFAAGGDSLLIAQLVSRVRAEEPLARAHPFDRLLRWALAQPTPAGLAQCLRDAKAQAQ
Sbjct  1612  QDFFAAGGDSLLIAQLVSRVRAEEPLARAHPFDRLLRWALAQPTPAGLAQCLRDAKAQAQ  1671

Query  1681  SSTAAPSR--TATATATATATATATAATAAARAAPAASIVHAPAPRGRIDVRPLRLAPGV  1738
             SSTAAPSR  TATATATATATATATAATAAARAAPAASIVHAPAPRGRIDVRPLRLAPGV
Sbjct  1672  SSTAAPSRTATATATATATATATATAATAAARAAPAASIVHAPAPRGRIDVRPLRLAPGV  1731

Query  1739  GVPRVVVHEGLGTVHAYRPIVPALARLGPVLGFAVRDAQDYLDLPAQHLSATLGARYAYA  1798
             GVPRVVVHEGLGTVHAYRPIVPALARLGPVLGFAVRDAQDYLDLPAQHLSATLGARYAYA
Sbjct  1732  GVPRVVVHEGLGTVHAYRPIVPALARLGPVLGFAVRDAQDYLDLPAQHLSATLGARYAYA  1791

Query  1799  LSREGIGEVDVLGYCSGGLIALEMAKTLVQLGVAVRSLDIVSSYRIPYLIEDERLVLFNF  1858
             LSREGIGEVDVLGYCSGGLIALEMAKTLVQLGVAVRSLDIVSSYRIPYLIEDERLVLFNF
Sbjct  1792  LSREGIGEVDVLGYCSGGLIALEMAKTLVQLGVAVRSLDIVSSYRIPYLIEDERLVLFNF  1851

Query  1859  AATLGLPLDALGFPAAHLLADALADALKADPARVAPGSIQAQLEAFGDRCAPLDILRRRV  1918
             AATLGLPLDALGFPAAHLLADALADALKADPARVAPGSIQAQLEAFGDRCAPLDILRRRV
Sbjct  1852  AATLGLPLDALGFPAAHLLADALADALKADPARVAPGSIQAQLEAFGDRCAPLDILRRRV  1911

Query  1919  LRVSAGLPADGDEPHPLVDERERLYRLFMHSVHASHWAAHAPYAGPLRLFVPERCNPLIP  1978
             LRVSAGLPADGDEPHPLVDERERLYRLFMHSVHASHWAAHAPYAGPLRLFVPERCNPLIP
Sbjct  1912  LRVSAGLPADGDEPHPLVDERERLYRLFMHSVHASHWAAHAPYAGPLRLFVPERCNPLIP  1971

Query  1979  QQRAALFDYWRDQALGGIALVDVPGGHFDCLTAAFVDTHLKEAR  2022
             QQRAALFDYWRDQALGGIALVDVPGGHFDCLTAAFVDTHLKEAR
Sbjct  1972  QQRAALFDYWRDQALGGIALVDVPGGHFDCLTAAFVDTHLKEAR  2015
```

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
         "% Identical Matches" =  percent_identical
          ) %>% 
  head(10) %>% 
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

Lastly, we can get the result of the summary of total number of sequence
that has to be removed from each data sets:

``` r
summary <- all_results[["all"]]

summary %>% 
  knitr::kable()

# save summary to .csv

write.csv(summary, "../../../data/getting-data-old/summary_count_data_removed.csv", row.names = FALSE
```

``` r
all_removed_rows <- all_results[["all_removed_rows"]]

removed_all_effector <- all_removed_rows %>% 
  filter(., label == 1)

removed_all_noneffector <- all_removed_rows %>% 
  filter(., label == 0)

removed_all_noneffector %>% 
  nrow()
```

    ## [1] 140

View how the data looks like

``` r
bacteria <- c("Acinetobacter baumannii", 
"Aeromonas hydrophila", 
"Aeromonas salmonicida", 
"Brucella abortus", 
"Burkholderia glumae", 
"Burkholderia pseudomallei", 
"Campylobacter jejuni", 
"Citrobacter rodentium", 
"Clavibacter michiganensis", 
"Coxiella burnetii", 
"Cystobacter fuscus", 
"Edwardsiella ictaluri", 
"Erwinia amylovora", 
"Escherichia coli", 
"Francisella tularensis", 
"Helicobacter pylori", 
"Legionella pneumophila",  
"Listeria monocytogenes", 
"Mycobacterium tuberculosis", 
"Pantoea stewartii", 
"Pseudomonas aeruginosa", 
"Pseudomonas cichorii", 
"Pseudomonas savastanoi", 
"Pseudomonas syringae", 
"Ralstonia solanacearum", 
"Salmonella enterica", 
"Shigella flexneri", 
"Staphylococcus aureus", 
"Vibrio parahaemolyticus", 
"Xanthomonas axonopodis", 
"Xanthomonas campestris", 
"Xanthomonas oryzae", 
"Xylella fastidiosa", 
"Yersinia enterocolitica", 
"Yersinia pestis", 
"Yersinia pseudotuberculosis")

fungi <- c("Beauveria bassiana",
"Blumeria graminis",
"Botrytis cinerea",
"Cercospora apii",
"Cercospora beticola",
"Colletotrichum orbiculare",
"Dothistroma septosporum",
"Fusarium oxysporum",
"Leptosphaeria maculans",
"Magnaporthe oryzae",
"Melampsora lini",
"Parastagonospora nodorum",
"Passalora fulva",
"Penicillium expansum",
"Pseudocercospora fuligena",
"Puccinia striiformis",
"Rhynchosporium commune",
"Verticillium dahliae",
"Ustilago maydis",
"Zymoseptoria tritici")

oomycetes <- c("Hyaloperonospora arabidopsidis", 
"Phytophthora cactorum", 
"Phytophthora capsici", 
"Phytophthora infestans", 
"Phytophthora parasitica", 
"Phytophthora sojae", 
"Pythium aphanidermatum", 
"Plasmopara halstedii", 
"Phytophthora megakarya")
```

``` r
effector_data <- data.table::fread("../../../data/getting-data-old/effector_with_IDs_organism.csv")

effector_removed_with_species <- removed_all_effector %>%
  left_join(effector_data, by = "sequence") %>% 
  select(sequence, label, pathogen_short) %>% 
  mutate(category = ifelse(pathogen_short %in% bacteria, "bacteria", 
                    ifelse(pathogen_short %in% fungi, "fungus", 
                    ifelse(pathogen_short %in% oomycetes, "oomycetes", "others")))) %>% 
  group_by(category) %>% 
  summarise(count = n())

noneffector_data <- data.table::fread("../../../data/getting-data-old/noneffector_with_IDs_organism.csv")

noneffector_removed_with_species <- removed_all_noneffector %>%
  left_join(noneffector_data, by = "sequence") %>% 
  unique() %>% 
  select(sequence, label, pathogen_short) %>%
  mutate(category = ifelse(pathogen_short %in% bacteria, "bacteria",
                    ifelse(pathogen_short %in% fungi, "fungus",
                    ifelse(pathogen_short %in% oomycetes, "oomycetes", "others")))) %>% 
  group_by(category) %>% 
  summarise(count = n())
```

``` r
effector_removed_with_species %>% 
  knitr::kable()
```

| category  |  count|
|:----------|------:|
| bacteria  |     27|
| fungus    |     11|
| oomycetes |      6|

``` r
noneffector_removed_with_species %>% 
  knitr::kable()
```

| category  |  count|
|:----------|------:|
| bacteria  |    122|
| fungus    |      2|
| oomycetes |     16|

Get a new additional effector and non-effector data
---------------------------------------------------

``` r
df_training_new_after_removed <- all_results[["removed_training"]][["df"]]
df_validation_new_after_removed <- all_results[["removed_validation"]][["df"]]
df_testing_new_after_removed <- all_results[["removed_testing"]][["df"]]

# save them to .csv format
write.csv(df_training_new_after_removed, "../../../data/getting-data-old/df_training_new_after_removed.csv", col.names = TRUE)
write.csv(df_testing_new_after_removed, "../../../data/getting-data-old/df_testing_new_after_removed.csv", col.names = TRUE)
write.csv(df_validation_new_after_removed, "../../../data/getting-data-old/df_validation_new_after_removed.csv", col.names = TRUE)
```

According to the blast results we obtained previously, there are some
data that need to be removed since they are identical with other protein
sequences.

### Load the data for each datasets that has been removed

``` r
df_train <- data.table::fread("../../../data/getting-data-old/df_training_new_after_removed.csv", drop = "V1")
df_val <- data.table::fread("../../../data/getting-data-old/df_validation_new_after_removed.csv", drop = "V1")
df_test <- data.table::fread("../../../data/getting-data-old/df_testing_new_after_removed.csv", drop = "V1")
```

### Load fasta data

After some process that has been done previously and resulting `.fasta`
file data, then now we can read the fasta file and add them to the data
frame without any identical protein sequence data.

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
# path 

add_effector_path <- "../../../data/getting-data-old/BLAST-data/additional-data-fasta/batch_entrez_effector.fasta"
add_noneffector_path <- "../../../data/getting-data-old/BLAST-data/additional-data-fasta/batch_entrez_noneffector.fasta"

add_effector_parsed <- parse_fasta_data_ncbi(add_effector_path)
add_noneffector_parsed <- parse_fasta_data_ncbi(add_noneffector_path)

# 
add_noneffector <- add_noneffector_parsed %>% 
  select(sequence) %>% 
  mutate(label = as.factor(0))

add_effector <- add_effector_parsed %>% 
  select(sequence) %>% 
  mutate(label = as.factor(1))

all_removed_rows_freq <- data.table::fread("../../../data/getting-data-old/all_removed_rows_freq.csv")

all_removed_rows_freq %>%
  knitr::kable()
```

``` r
# Add the data of non-effector and effector to the training data

df_train <- df_train %>%
  rbind(., add_noneffector[1:76, ]) %>% 
  rbind(., add_effector[1:21, ])

df_val <- df_val %>% 
  rbind(., add_noneffector[77:112, ]) %>% 
  rbind(., add_effector[22:27, ])

df_test <- df_test %>% 
  rbind(., add_noneffector[113:140, ]) %>% 
  rbind(., add_effector[28:38, ])

write.csv(df_train, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-data.csv", col.names = TRUE)
write.csv(df_val, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-data.csv", col.names = TRUE)
write.csv(df_test, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-data.csv", col.names = TRUE)
```

### View the new data

``` r
library(tidyverse)
df_train <- data.table::fread("../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-data.csv")

df_train_sequence <- df_train %>% 
  select(sequence) 

# write.table(df_train_sequence, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-sequence.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_train_label <- df_train %>% 
  select(label)

# write.table(df_train_label, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-label.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_val <- data.table::fread("../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-data.csv")

df_val_sequence <- df_val %>% 
  select(sequence)

# write.table(df_val_sequence, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-sequence.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_val_label <- df_val %>% 
  select(label)

# write.table(df_val_label, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-label.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)


df_test <- data.table::fread("../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-data.csv")

df_test_sequence <- df_test %>% 
  select(sequence)
# write.table(df_test_sequence, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-sequence.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_test_label <- df_test %>% 
  select(label)
# write.table(df_test_label, "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-label.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

#### Training data

``` r
nrow(df_train)
```

    ## [1] 578

``` r
nrow(df_val)
```

    ## [1] 193

``` r
nrow(df_test)
```

    ## [1] 193
