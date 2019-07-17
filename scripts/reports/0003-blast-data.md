BLAST different datasets of data
================================

Introduction
------------

In order to make sure that we do not have identical sequence among the
datasets, we need to BLAST the datasets. However, first of all, we need
to split the data into training, development, and test datasets.

Split the data into training, development, and test datasets
------------------------------------------------------------

As alternative of scikit-learn that we use in Python, we can use caret
package.

``` r
library(caret)
```

    ## Loading required package: lattice

    ## Loading required package: ggplot2

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ tibble  2.1.3     ✔ purrr   0.3.2
    ## ✔ tidyr   0.8.3     ✔ dplyr   0.8.3
    ## ✔ readr   1.3.1     ✔ stringr 1.4.0
    ## ✔ tibble  2.1.3     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ purrr::lift()   masks caret::lift()

``` r
# Read CSV of the effector data
effector_data <- data.table::fread("effector_data.csv") %>% 
  dplyr::select(-V1) %>% 
  mutate(data = 1) %>% 
  dplyr::select(Sequence, data) %>% 
  rename(sequence = Sequence)

# Read CSV of the noneffector data 
non_effector_data <- data.table::fread("non_effector_data.csv") %>% 
  dplyr::select(-rowid) %>%
  mutate(data = 0) %>% 
  dplyr::select(sequence, data)
```

``` r
#  Combine both data into same dataframe rowwise()

data_to_split <- effector_data %>% 
  rbind(non_effector_data)
```

``` r
set.seed(100)  # For reproducibility

# Create index for testing and training data
training_id <- createDataPartition(y = data_to_split$data, 
                               p = 0.6, list = FALSE)

# subset power_plant data to training
training <- data_to_split[training_id,]


# subset the rest to test
rest <- data_to_split[-training_id,]

# Splitting the rest into two different class of datasets

val_id <- createDataPartition(y = rest$data, 
                               p = 0.5, list = FALSE)

validation <- rest[val_id,]

testing <- rest[-val_id,]
```

``` r
rbind("Training set" = nrow(training)/nrow(data_to_split),
      "Validation set" = nrow(validation)/nrow(data_to_split), 
      "Testing set" = nrow(testing)/nrow(data_to_split)) %>% 
       round(2)*100
```

    ##                [,1]
    ## Training set     60
    ## Validation set   20
    ## Testing set      20

Take the label into separated data frame

``` r
#  Training data
training_seq <- training %>% 
  dplyr::select(sequence)

training_label <- training %>% 
  dplyr::select(data)

# Validation data 
validation_seq <- validation %>% 
  dplyr::select(sequence)

validation_label <- validation %>% 
  dplyr::select(data) 

# Testing data 
testing_seq <- testing %>% 
  dplyr::select(sequence)

testing_label <- testing  %>% 
  dplyr::select(data)
```

``` r
# Save the data to CSV format

# Sequence Input
write_csv(training_seq, "data-sets/training_input.csv", col_names = FALSE)
write_csv(validation_seq, "data-sets/validation_input.csv", col_names = FALSE)
write_csv(testing_seq, "data-sets/testing_input.csv", col_names = FALSE)

# Label data
write_csv(training_label, "data-sets/training_label.csv", col_names = FALSE)
write_csv(validation_label, "data-sets/validation_label.csv", col_names = FALSE)
write_csv(testing_label, "data-sets/testing_label.csv", col_names = FALSE)
```

Change the datasets to .fasta file
----------------------------------

To change the data into fasta file, we can run a python scripts called:

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
python ConvertFASTA.py training_input.csv training_input.fasta
```

The above syntax is the example when we convert our training sequence
data from .csv to .fatsa file.

BLAST the datasets
------------------

After getting the datasets in fasta file, now we can get BLAST the
datasets as follows:

### Making the database

Now let us BLAST two datasets of sequence (validation and testing)
against the entire training datasets. First, we need to tell BLAST that
are (a) a database, and (b) a protein database. That’s done by calling
‘makeblastdb’:

``` bash
makeblastdb -in training_input.csv -dbtype prot
```

### BLAST the validation dataset against training datasets

``` bash
blastp -query validation_input.fasta -db training_input.fasta -out val_vs_training.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

### BLAST the validation dataset against training datasets

``` bash
blastp -query testing_input.fasta -db training_input.fasta -out testing_vs_training.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

Read the BLAST results
----------------------

In order to read the BLAST results we obtained from running `blastp`, we
can use the function below:

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
```

``` r
validation_VS_training_path <- "data-sets/val_vs_training.tsv"
testing_VS_training_path <- "data-sets/testing_vs_training.tsv"

validation_VS_training_results <- blast_results(validation_VS_training_path)
testing_VS_training_results <- blast_results(testing_VS_training_path)
```

View the results

``` r
validation_VS_training_results[["df_identical_protein"]] %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_indentical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|--------------------:|------------------:|
|       8|   113|     282|   113|     113|     111|         2|       111|             98.23009|           98.23009|
|       8|   113|     347|   113|     113|     111|         2|       112|             98.23009|           99.11504|
|       8|   113|     346|   113|     113|     110|         3|       111|             97.34513|           98.23009|
|       8|   113|     281|   113|     113|     109|         4|       110|             96.46018|           97.34513|
|      36|   311|     157|   311|     311|     308|         3|       309|             99.03537|           99.35691|
|      36|   311|     131|   311|     311|     306|         5|       307|             98.39228|           98.71383|
|      49|  1028|     208|  1096|    1096|    1005|        23|      1015|             91.69708|           92.60949|
|      59|  1163|     136|  1164|    1164|    1124|        39|      1141|             96.56357|           98.02406|
|      59|  1163|     209|  1126|    1164|    1100|        25|      1111|             94.58298|           95.52880|
|      75|   309|     260|   296|     296|     296|         0|       296|             95.79288|           95.79288|
|     101|   160|     343|   160|     160|     159|         1|       160|             99.37500|          100.00000|
|     104|   113|     282|   113|     113|     112|         1|       112|             99.11504|           99.11504|
|     104|   113|     347|   113|     113|     112|         1|       113|             99.11504|          100.00000|
|     104|   113|     346|   113|     113|     111|         2|       112|             98.23009|           99.11504|
|     104|   113|     281|   113|     113|     110|         3|       111|             97.34513|           98.23009|

``` r
testing_VS_training_results[["df_identical_protein"]] %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_indentical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|--------------------:|------------------:|
|      18|   113|     282|   113|     113|     112|         1|       113|             99.11504|          100.00000|
|      18|   113|     347|   113|     113|     112|         1|       112|             99.11504|           99.11504|
|      18|   113|     281|   113|     113|     112|         1|       112|             99.11504|           99.11504|
|      18|   113|     346|   113|     113|     111|         2|       111|             98.23009|           98.23009|
|      40|   720|     120|   720|     720|     694|        26|       711|             96.38889|           98.75000|
|      58|   187|     175|   187|     187|     176|        11|       180|             94.11765|           96.25668|
|      68|  1795|     184|  1795|    1795|    1792|         3|      1792|             99.83287|           99.83287|
|      72|  1096|     208|  1096|    1096|    1076|        20|      1086|             98.17518|           99.08759|
|      72|  1096|     209|  1126|    1129|    1071|        22|      1083|             95.11545|           96.18117|
|      72|  1096|     103|  1095|    1096|     991|       104|      1025|             90.41971|           93.52190|
|      88|   309|     260|   296|     296|     296|         0|       296|             95.79288|           95.79288|
|      89|   309|     260|   296|     296|     296|         0|       296|             95.79288|           95.79288|
|      91|   799|     269|   863|     799|     799|         0|       799|             92.58401|           92.58401|
|      98|   267|     309|   267|     267|     263|         4|       264|             98.50187|           98.87640|
|     116|   446|     355|   446|     446|     445|         1|       446|             99.77578|          100.00000|

#### Validation against Training Results

``` r
# Get all index of the training data that are identical to the validation data
training_id_from_val <- validation_VS_training_results[["subject_index"]]
```

``` r
# View all of the index that needs to be removed 
id_to_be_removed_validation <- validation_VS_training_results[["query_index"]]
```

#### Validation against Testing Results

``` r
# View all of the index that needs to be removed 
training_id_from_test <- testing_VS_training_results[["subject_index"]]
```

``` r
# View all of the index that needs to be removed 
id_to_be_removed_testing <- testing_VS_training_results[["query_index"]]
```

``` r
# View all of the index that needs to be removed 
training_id_from_val <- validation_VS_training_results[["subject_index"]]
```

``` r
# View all of the index that needs to be removed 
training_id_from_test <- testing_VS_training_results[["subject_index"]]
```

``` r
id_to_be_removed_training <- c(training_id_from_val, training_id_from_test) %>%
    unique() %>%
    unlist() 
```

``` r
#  MAking a dataframe from all results 

after_blast_summary <- data.frame("Total" = c(nrow(training_seq), nrow(validation_seq), nrow(testing_seq))) %>% 
  mutate(Removed = c(length(id_to_be_removed_training), length(id_to_be_removed_validation), length(id_to_be_removed_testing))) %>% 
  mutate(Datasets = c("Training", "Validation", "Testing")) %>% 
  select(Datasets, everything())
```

``` r
after_blast_summary %>% 
  knitr::kable()
```

| Datasets   |  Total|  Removed|
|:-----------|------:|--------:|
| Training   |    480|       18|
| Validation |    160|        7|
| Testing    |    160|       10|
