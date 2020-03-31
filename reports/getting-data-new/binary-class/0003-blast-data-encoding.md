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

non_effector_data
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
```

``` r
validation_VS_training_path <- here::here("data/getting-data-new/binary-class-data/data-sets", "val_vs_training.tsv")
testing_VS_training_path <- here::here("data/getting-data-new/binary-class-data/data-sets", "testing_vs_training.tsv")

validation_VS_training_results <- blast_results(validation_VS_training_path)
testing_VS_training_results <- blast_results(testing_VS_training_path)
```

View the results

``` r
validation_VS_training_results[["df_identical_protein"]] %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
|       8|   113|     282|   113|     113|     111|         2|       111|            98.23009|           98.23009|
|       8|   113|     347|   113|     113|     111|         2|       112|            98.23009|           99.11504|
|       8|   113|     346|   113|     113|     110|         3|       111|            97.34513|           98.23009|
|       8|   113|     281|   113|     113|     109|         4|       110|            96.46018|           97.34513|
|      36|   311|     157|   311|     311|     308|         3|       309|            99.03537|           99.35691|
|      36|   311|     131|   311|     311|     306|         5|       307|            98.39228|           98.71383|
|      49|  1028|     208|  1096|    1096|    1005|        23|      1015|            91.69708|           92.60949|
|      59|  1163|     136|  1164|    1164|    1124|        39|      1141|            96.56357|           98.02406|
|      59|  1163|     209|  1126|    1164|    1100|        25|      1111|            94.58298|           95.52880|
|      75|   309|     260|   296|     296|     296|         0|       296|            95.79288|           95.79288|
|     101|   160|     343|   160|     160|     159|         1|       160|            99.37500|          100.00000|
|     104|   113|     282|   113|     113|     112|         1|       112|            99.11504|           99.11504|
|     104|   113|     347|   113|     113|     112|         1|       113|            99.11504|          100.00000|
|     104|   113|     346|   113|     113|     111|         2|       112|            98.23009|           99.11504|
|     104|   113|     281|   113|     113|     110|         3|       111|            97.34513|           98.23009|

``` r
testing_VS_training_results[["df_identical_protein"]] %>% 
  knitr::kable()
```

|  qseqid|  qlen|  sseqid|  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|-------:|-----:|-------:|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
|      18|   113|     282|   113|     113|     112|         1|       113|            99.11504|          100.00000|
|      18|   113|     347|   113|     113|     112|         1|       112|            99.11504|           99.11504|
|      18|   113|     281|   113|     113|     112|         1|       112|            99.11504|           99.11504|
|      18|   113|     346|   113|     113|     111|         2|       111|            98.23009|           98.23009|
|      40|   720|     120|   720|     720|     694|        26|       711|            96.38889|           98.75000|
|      58|   187|     175|   187|     187|     176|        11|       180|            94.11765|           96.25668|
|      68|  1795|     184|  1795|    1795|    1792|         3|      1792|            99.83287|           99.83287|
|      72|  1096|     208|  1096|    1096|    1076|        20|      1086|            98.17518|           99.08759|
|      72|  1096|     209|  1126|    1129|    1071|        22|      1083|            95.11545|           96.18117|
|      72|  1096|     103|  1095|    1096|     991|       104|      1025|            90.41971|           93.52190|
|      88|   309|     260|   296|     296|     296|         0|       296|            95.79288|           95.79288|
|      89|   309|     260|   296|     296|     296|         0|       296|            95.79288|           95.79288|
|      91|   799|     269|   863|     799|     799|         0|       799|            92.58401|           92.58401|
|      98|   267|     309|   267|     267|     263|         4|       264|            98.50187|           98.87640|
|     116|   446|     355|   446|     446|     445|         1|       446|            99.77578|          100.00000|

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

training_id_from_test
```

    ##  sseqid1  sseqid2  sseqid3  sseqid4  sseqid5  sseqid6  sseqid7  sseqid8 
    ##      282      347      281      346      120      175      184      208 
    ##  sseqid9 sseqid10 sseqid11 sseqid12 sseqid13 sseqid14 
    ##      209      103      260      269      309      355

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
# Take the intersect data from BLAST results training vs val and training vs testing
id_to_be_removed_training <- c(training_id_from_val, training_id_from_test) %>%
    unique() %>%
    unlist() 
```

Load all of the data

``` r
# Header = FALSE since the first column of the data is not the header

# Training data
training_seq <- here::here("data/getting-data-new/binary-class-data/data-sets", "training_input.csv") %>%
  data.table::fread(header = FALSE) %>% 
  rename(sequence = V1)

training_label <- here::here("data/getting-data-new/binary-class-data/data-sets", "training_label.csv") %>%
  data.table::fread(header = FALSE) %>% 
  rename(label = V1)

# Validation data
validation_seq <- here::here("data/getting-data-new/binary-class-data/data-sets", "validation_input.csv") %>%
  data.table::fread(header = FALSE) %>% 
  rename(sequence = V1)

validation_label <- here::here("data/getting-data-new/binary-class-data/data-sets", "validation_label.csv") %>%
  data.table::fread(header = FALSE) %>% 
  rename(label = V1)

# Testing data
testing_seq <- here::here("data/getting-data-new/binary-class-data/data-sets", "testing_input.csv") %>%
  data.table::fread(header = FALSE) %>% 
  rename(sequence = V1)

testing_label <- here::here("data/getting-data-new/binary-class-data/data-sets", "testing_label.csv") %>%
  data.table::fread(header = FALSE) %>% 
  rename(label = V1)
```

``` r
#  Making a dataframe from all results 
after_blast_summary <- data.frame(
  "Total" = c(nrow(training_seq), nrow(validation_seq), nrow(testing_seq))
) %>%
  mutate(
    "To be removed" = c(length(id_to_be_removed_training), length(id_to_be_removed_validation), length(id_to_be_removed_testing)),
    Datasets = c("Training", "Validation", "Testing")
  ) %>%
  select(Datasets, everything())
```

``` r
# Print the summary of the datasets (Total data and data needed to be removed)
after_blast_summary %>% 
  knitr::kable()
```

| Datasets   |  Total|  To be removed|
|:-----------|------:|--------------:|
| Training   |    480|             18|
| Validation |    160|              7|
| Testing    |    160|             10|

Since the percentage of data that need to be removed is considered to be
quite small then we do not need to replace with the new data. Now we
just need to remove the data that are identical by using the index we
already had.

``` r
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
# Get the new dataframe after removing the identical sequence data
training <- remove_from_data_sets("data-sets/training_input.csv", "data-sets/training_label.csv", id_to_be_removed_training)[["df"]]

validation <- remove_from_data_sets("data-sets/validation_input.csv", "data-sets/validation_label.csv", id_to_be_removed_validation)[["df"]]

testing <- remove_from_data_sets("data-sets/testing_input.csv", "data-sets/testing_label.csv", id_to_be_removed_testing)[["df"]]

# Save the dataframe to CSV
# write_csv(training, "data-sets/new_training_data.csv")
# write_csv(validation, "data-sets/new_validation_data.csv")
# write_csv(testing, "data-sets/new_testing_data.csv")
```

Encoding the data using Pandas
------------------------------

``` r
# Load keras library
library(keras)
```

Import the library that we need to reprocess the data.

``` python
import pandas as pd
import numpy as np
```

``` python
import os
cwd = os.getcwd()
print(cwd)
```

### Load the data

``` python
training_data = pd.read_csv('data-sets/new_training_data.csv', index_col = False)
validation_data = pd.read_csv('data-sets/new_validation_data.csv', index_col = False)
testing_data = pd.read_csv('data-sets/new_testing_data.csv', index_col = False)
```

``` python
# Define the input and the label of data 

# Training datasets
input_train = training_data[["sequence"]]
label_train = training_data[["label"]]

# Validation datasets
input_val = validation_data[["sequence"]]
label_val = validation_data[["label"]]

# Testing data 
input_test= testing_data[["sequence"]]
label_test = testing_data[["label"]]
```

``` python
training_data[["sequence"]]
```

``` python
from collections import Counter
field_length_train = input_train.sequence.astype(str).map(len) 
field_length_val = input_val.sequence.astype(str).map(len)
field_length_test = input_test.sequence.astype(str).map(len) 

print(max(field_length_train)) 
print(max(field_length_val)) 
print(max(field_length_test))
```

``` python
import matplotlib.pyplot as plt
plt.clf()
plt.hist(field_length_train, bins=100)
plt.ylabel('Count')
plt.xlabel('Length')
plt.title('Histogram of Length of Sequences')
plt.show()
```

### One hot encoding

``` python
def get_key(mydict, element):
    key = list(mydict.keys())[list(mydict.values()).index(element)]
    return(key)

amino = ['R', 'K', 'D', 'E', 'Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'W', 'A', 'I', 'L', 'M', 'F', 'V', 'P', 'G']
token_index = dict(zip(range(1, (len(amino)+1)), amino))

max_length = 4034
def get_encoding(mydata, max_length):
    results = np.zeros((len(mydata), max_length, max(token_index.keys())))
    for i, sample in enumerate(mydata):
        for j, character in enumerate(sample):
            if character in token_index.values():
                index = get_key(token_index, character) - 1
                results[i, j, index] = 1. 
            else:
                results[i, j, :] = results[i, j, :]
    return results
```

``` python
# Change the data to list
x_train = input_train.sequence.tolist()
x_val = input_val.sequence.tolist()
x_test = input_test.sequence.tolist()
```

``` python
# Encoding by calling the function get_encoding()
one_hot_train = get_encoding(x_train, max_length)
one_hot_val = get_encoding(x_val, max_length)
one_hot_test = get_encoding(x_test, max_length)
```

``` python
input_train.info()
# print(one_hot_test[1:2, :20, :20])
```

#### Change the label into list data format

``` python
# Change the data into 
y_train = label_train.label.tolist()
y_val = label_val.label.tolist()
y_test = label_test.label.tolist()
```

#### Save all of the data

``` python
# Save the input data
np.save('data-sets/x_train.npy', one_hot_train)
np.save('data-sets/x_val.npy', one_hot_val)
np.save('data-sets/x_test.npy', one_hot_test)

# Save the label data 
np.save('data-sets/y_train.npy', y_train)
np.save('data-sets/y_val.npy', y_val)
np.save('data-sets/y_test.npy', y_test)
```

### Integer Encoding

``` python
def get_value(mydict, element):
    key = mydict.get(element)
    return(key)
    
dic = {'A':1, 'B':22, 'U':23,'J':24,'Z':25,'O':26,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20,'X':21}
max_len = 4034

def encoding_each_sequence(input):
    results = np.zeros((len(input), max_len))
    for i, sample in enumerate(input):
        for j, character in enumerate(sample):
            results[i, j] = get_value(dic, character)
    return results
```

``` python
# Encoding by calling the function get_encoding()
integer_train = encoding_each_sequence(x_train)
integer_val = encoding_each_sequence(x_val)
integer_test = encoding_each_sequence(x_test)
```

``` python
print(integer_train.shape)
```

Save the results

``` python
np.save('data-sets/x_train_int.npy', integer_train)
np.save('data-sets/x_val_int.npy', integer_val)
np.save('data-sets/x_test_int.npy', integer_test)
```
