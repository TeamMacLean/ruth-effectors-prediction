Checking the final result of the data using BLAST
=================================================

In order to make sure that there is no identical protein between the
protein datasets, we need again to blast the tetsing and validation
datasets against training datasets.

However, since we only have the sequence in .csv format, we need to
change the file to fasta

``` python
import sys

#File input
fileInput = open(sys.argv[1], "r")

#File output
fileOutput = open(sys.argv[2], "w")

#Seq count
count = 1 ;

#Loop through each line in the input file
print "Converting to FASTA..."
for strLine in fileInput:

    #Strip the endline character from each input line
    strLine = strLine.rstrip("\n")

    #Output the header
    fileOutput.write(">" + str(count) + "\n")
    fileOutput.write(strLine + "\n")

    count = count + 1
print ("Done.")

#Close the input and output file
fileInput.close()
fileOutput.close()
```

And we can call the scripts converter above in terminal using command
line

``` bash
python ConvertFASTA.py 0003-new-data-sets/training-sequence.csv 0003-new-data-sets/training.fasta
```

The above syntax is the example when we convert our training sequence
data from .csv to .fatsa file.

Check the data using BLAST
--------------------------

After converting those files to .fatsa, now we can directly BLAST the
fasta file using `blastp`.

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

Identify the results
--------------------

``` r
add_blast_train_x_val <- blast_results("../../../data/getting-data-old/BLAST-data/0003-new-data-sets/blast_train_x_val.tsv")
add_blast_train_x_test <- blast_results("../../../data/getting-data-old/BLAST-data/0003-new-data-sets/blast_train_x_test.tsv")
```

``` r
add_blast_train_x_val[["df_identical_protein"]] %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
|      39|  1795|     559|  1767|    1800|    1714|        48|      1728|            95.48747|           96.26741|
|     157|  2334|     508|  2347|    2338|    2317|        17|      2328|            98.72177|           99.19046|
|     157|  2334|     500|  2347|    2338|    2315|        19|      2326|            98.63656|           99.10524|
|     157|  2334|     538|  2347|    2338|    2304|        30|      2321|            98.16787|           98.89220|
|     157|  2334|     533|  2347|    2338|    2304|        30|      2321|            98.16787|           98.89220|
|     157|  2334|     501|  2347|    2338|    2304|        30|      2321|            98.16787|           98.89220|
|     157|  2334|     514|  2347|    2338|    2300|        34|      2315|            97.99744|           98.63656|
|     157|  2334|     536|  2347|    2338|    2300|        34|      2320|            97.99744|           98.84960|
|     157|  2334|     526|  2347|    2338|    2300|        34|      2320|            97.99744|           98.84960|
|     157|  2334|     504|  2347|    2338|    2300|        34|      2320|            97.99744|           98.84960|
|     157|  2334|     530|  2345|    2338|    2283|        49|      2310|            97.35608|           98.50746|
|     157|  2334|     524|  2345|    2338|    2281|        51|      2309|            97.27079|           98.46482|

``` r
add_blast_train_x_test[["df_identical_protein"]] %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
|      96|  1795|     559|  1767|    1800|    1717|        45|      1731|            95.65460|           96.43454|
|     191|  1388|     345|  1388|    1388|    1258|       130|      1316|            90.63401|           94.81268|
|     192|  1411|     566|  1411|    1411|    1283|       128|      1343|            90.92842|           95.18072|
|     192|  1411|     565|  1411|    1411|    1283|       128|      1343|            90.92842|           95.18072|
|     192|  1411|     564|  1411|    1411|    1283|       128|      1343|            90.92842|           95.18072|
|     192|  1411|     563|  1411|    1411|    1281|       130|      1339|            90.78668|           94.89724|
|     192|  1411|     562|  1411|    1411|    1281|       130|      1339|            90.78668|           94.89724|
|     192|  1411|     561|  1411|    1411|    1281|       130|      1339|            90.78668|           94.89724|
|     193|  1388|     566|  1411|    1388|    1388|         0|      1388|            98.36995|           98.36995|
|     193|  1388|     565|  1411|    1388|    1388|         0|      1388|            98.36995|           98.36995|
|     193|  1388|     564|  1411|    1388|    1388|         0|      1388|            98.36995|           98.36995|
|     193|  1388|     563|  1411|    1388|    1367|        21|      1370|            96.88164|           97.09426|
|     193|  1388|     562|  1411|    1388|    1367|        21|      1370|            96.88164|           97.09426|
|     193|  1388|     561|  1411|    1388|    1367|        21|      1370|            96.88164|           97.09426|
|     193|  1388|     576|  1411|    1388|    1286|       102|      1329|            91.14103|           94.18852|
|     193|  1388|     578|  1388|    1388|    1279|       109|      1323|            92.14697|           95.31700|
|     193|  1388|     558|  1388|    1388|    1276|       112|      1326|            91.93084|           95.53314|

``` r
#  checking all of the results for identical protein idetification

blast_train_x_val_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/blast_train_x_val.tsv"
blast_train_x_test_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/blast_train_x_test.tsv"
input_train_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-sequence.csv"
label_train_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-label.csv"
input_val_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-sequence.csv"
label_val_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-label.csv"
input_test_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-sequence.csv"
label_test_path = "../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-label.csv"

check_add_data <- get_info(blast_train_x_val_path,
                     blast_train_x_test_path,
                     input_train_path,
                     label_train_path,
                     input_val_path,
                     label_val_path,
                     input_test_path,
                     label_test_path)
```
