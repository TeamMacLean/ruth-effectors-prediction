Split and Encode data after all sampled and having no identical protein
=======================================================================

Background
----------

After having all data sampled and making sure there is no identical
protein (more than 90% identical), then the data is ready to be
processed (splitted and encoded)

Functions
---------

### Load libraries

``` r
library(docstring)
library(tidyverse)
# library(taxize)
library(caret)
# reticulate::use_condaenv(condaenv = "tensorflow2", conda = "/anaconda3/bin/conda")

# Get the source of the function used to split the data
source(here::here("scripts/r-scripts/r-functions", "split_datasets.R"))
```

``` r
# Load keras library
library(keras)
```

Import the library that we need to reprocess the data.

``` python
import pandas as pd
import numpy as np
```

### Define functions

#### Get sequences for each class

``` r
# Funtion to get each pathogen class given a full table of effector data
get_seq_each_class <- function(df_effector, class_var) {
  df_seq <- df_effector %>%
    dplyr::filter(class == class_var) %>%
    dplyr::select(Sequence)

  return(df_seq)
}
```

#### One hot encoding funtion

``` python
# Funtion to get the index of each character
def get_key(mydict, element):
    key = list(mydict.keys())[list(mydict.values()).index(element)]
    return(key)

# List of amino acids to encode
amino = ['R', 'K', 'D', 'E', 'Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'W', 'A', 'I', 'L', 'M', 'F', 'V', 'P', 'G']
token_index = dict(zip(range(1, (len(amino)+1)), amino))


max_length = 934 # Max sequence on the validation data
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

### Load the data

``` r
effector_final_after_blast <- readRDS("../../../data/secreted_data/data_processed_after_signalp/effector_final_after_blast.RDS")

effector_seq_fungi <- get_seq_each_class(effector_final_after_blast, class_var = "fungi")
effector_seq_bacteria <- get_seq_each_class(effector_final_after_blast, class_var = "bacteria")
effector_seq_oomycete <- get_seq_each_class(effector_final_after_blast, class_var = "oomycete")
```

### Get the non-effector from randomly sampling data

``` r
non_effector_seq_fungi <- readRDS("../../../data/secreted_data/ready_to_process/sampled_data_without_identical/fungi_sampled_table_good.RDS") %>% 
  dplyr::select(sequence) %>% 
  `colnames<-`("Sequence")
```

``` r
non_effector_seq_bacteria <- readRDS("../../../data/secreted_data/ready_to_process/sampled_data_without_identical/bacteria_sampled_table_good.RDS")  %>% 
  dplyr::select(sequence) %>% 
  `colnames<-`("Sequence")
```

``` r
non_effector_seq_oomycete <- readRDS("../../../data/secreted_data/ready_to_process/sampled_data_without_identical/oomycete_sampled_table_good.RDS")  %>% 
  dplyr::select(sequence) %>% 
  `colnames<-`("Sequence")
```

Splitting data
--------------

### For each pathogen organism

### Fungi

``` r
# Combine and labeled data
fungi_full_datasets <- get_data_labeled_binary(effector_seq_fungi, non_effector_seq_fungi)

# Splitted data
fungi_splitted <- get_data_splitted(fungi_full_datasets, p1 = 0.6, p2 = 0.2, test_dataset = TRUE)

fungi_splitted[[1]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/fungi_training.csv")

fungi_splitted[[2]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/fungi_validation.csv")

fungi_splitted[[3]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/fungi_testing.csv")
```

### Oomycete

``` r
oomycete_full_datasets <- get_data_labeled_binary(effector_seq_oomycete, non_effector_seq_oomycete)

oomycete_splitted <- get_data_splitted(oomycete_full_datasets, p1 = 0.6, p2 = 0.2, test_dataset = TRUE)

oomycete_splitted[[1]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/oomycete_training.csv")

oomycete_splitted[[2]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/oomycete_validation.csv")

oomycete_splitted[[3]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/oomycete_testing.csv")
```

### Bacteria

``` r
bacteria_full_datasets <- get_data_labeled_binary(effector_seq_bacteria, non_effector_seq_bacteria)

bacteria_splitted <- get_data_splitted(bacteria_full_datasets, p1 = 0.6, p2 = 0.2, test_dataset = TRUE)

bacteria_splitted[[1]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/bacteria_training.csv")

bacteria_splitted[[2]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/bacteria_validation.csv")

bacteria_splitted[[3]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/bacteria_testing.csv")
```

### Data all together

``` r
# Get both effector and non-effector ready
effector_final <- effector_final_after_blast %>% 
  dplyr::select(Sequence)

non_effector_final <- rbind(non_effector_seq_bacteria, non_effector_seq_fungi, non_effector_seq_oomycete) 

# Combine and label all of the effector and non-effector data 
full_data <- get_data_labeled_binary(effector_final, non_effector_final)

# Split datasets
full_data_splitted <- get_data_splitted(full_data, p1 = 0.6, p2 = 0.2, test_dataset = TRUE)

full_data_splitted[[1]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/all_training.csv")

full_data_splitted[[2]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/all_validation.csv")

full_data_splitted[[3]] %>% 
  data.table::fwrite("../../../data/secreted_data/ready_to_process/splitted-data/all_testing.csv")
```

Encode Data
-----------

### Oomycete

``` python
# Oomycete
training_oomycete = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/oomycete_training.csv", index_col = False)
validation_oomycete = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/oomycete_validation.csv", index_col = False)
testing_oomycete = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/oomycete_testing.csv", index_col = False)
```

``` python
training_oomycete.head(2)
```

    ##                                             Sequence  label
    ## 0  MVKLYCAVVGVAGSAFSVRVDESDTVDDLKDAIKAKKPNDFKDIDA...      1
    ## 1  MRLTNTLVVAVAAILLASENAFSAATDADQATVSKLAAAEFDTLVD...      1

``` python
# Define the input and the label of data 

# Training datasets
input_train_oomycete = training_oomycete[["Sequence"]]
label_train_oomycete = training_oomycete[["label"]]

# Validation datasets
input_val_oomycete = validation_oomycete[["Sequence"]]
label_val_oomycete = validation_oomycete[["label"]]

# Testing data 
input_test_oomycete = testing_oomycete[["Sequence"]]
label_test_oomycete = testing_oomycete[["label"]]
```

``` python
# To get the information about the length of the data

from collections import Counter
field_length_train_oomycete = input_train_oomycete.Sequence.astype(str).map(len) 
field_length_val_oomycete = input_val_oomycete.Sequence.astype(str).map(len)
field_length_test_oomycete = input_test_oomycete.Sequence.astype(str).map(len) 

print(max(field_length_train_oomycete)) 
```

    ## 820

``` python
print(max(field_length_val_oomycete)) 
```

    ## 695

``` python
print(max(field_length_test_oomycete))
```

    ## 536

``` python
# Change the data to list
x_train_oomycete = input_train_oomycete.Sequence.tolist()
x_val_oomycete = input_val_oomycete.Sequence.tolist()
x_test_oomycete = input_test_oomycete.Sequence.tolist()
```

``` python
# Define the maximum list
max_length = 820

# Encoding by calling the function get_encoding()
one_hot_train_oomycete = get_encoding(x_train_oomycete, max_length)
one_hot_val_oomycete = get_encoding(x_val_oomycete, max_length)
one_hot_test_oomycete = get_encoding(x_test_oomycete, max_length)
```

``` python
# View the encoding results
input_train_oomycete.info()
```

    ## <class 'pandas.core.frame.DataFrame'>
    ## RangeIndex: 102 entries, 0 to 101
    ## Data columns (total 1 columns):
    ## Sequence    102 non-null object
    ## dtypes: object(1)
    ## memory usage: 896.0+ bytes

``` python
print(one_hot_test_oomycete[1:2, :20, :20])
```

    ## [[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.]
    ##   [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]]]

#### Change the label into list data format

``` python
# Change the data into 
y_train_oomycete = label_train_oomycete.label.tolist()
y_val_oomycete = label_val_oomycete.label.tolist()
y_test_oomycete = label_test_oomycete.label.tolist()
```

``` python
# View the label data

print(y_train_oomycete)
```

    ## [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#### Save all of the Oomycete all of the data encoded

``` python
# Save the input data
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_train_oomycete.npy', one_hot_train_oomycete)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_val_oomycete.npy', one_hot_val_oomycete)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_test_oomycete.npy', one_hot_test_oomycete)

# Save the label data 
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_train_oomycete.npy', y_train_oomycete)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_val_oomycete.npy', y_val_oomycete)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_test_oomycete.npy', y_test_oomycete)
```

### Fungi

``` python
# fungi
training_fungi = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/fungi_training.csv", index_col = False)
validation_fungi = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/fungi_validation.csv", index_col = False)
testing_fungi = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/fungi_testing.csv", index_col = False)
```

``` python
training_fungi.head(2)
```

    ##                                             Sequence  label
    ## 0  MQSVLLLTVLTQSFIATASPLVERSTPLSFAEKRPQKVSYDWTTPY...      1
    ## 1  MRLANFLFYLAPMIVSSLAFDFVPLSGELDFSQEMVFINLTQQQFS...      1

``` python
# Define the input and the label of data 

# Training datasets
input_train_fungi = training_fungi[["Sequence"]]
label_train_fungi = training_fungi[["label"]]

# Validation datasets
input_val_fungi = validation_fungi[["Sequence"]]
label_val_fungi = validation_fungi[["label"]]

# Testing data 
input_test_fungi = testing_fungi[["Sequence"]]
label_test_fungi = testing_fungi[["label"]]
```

``` python
# To get the information about the length of the data

from collections import Counter
field_length_train_fungi = input_train_fungi.Sequence.astype(str).map(len) 
field_length_val_fungi = input_val_fungi.Sequence.astype(str).map(len)
field_length_test_fungi = input_test_fungi.Sequence.astype(str).map(len) 

print(max(field_length_train_fungi)) 
```

    ## 709

``` python
print(max(field_length_val_fungi)) 
```

    ## 4034

``` python
print(max(field_length_test_fungi))
```

    ## 690

``` python
# Change the data to list
x_train_fungi = input_train_fungi.Sequence.tolist()
x_val_fungi = input_val_fungi.Sequence.tolist()
x_test_fungi = input_test_fungi.Sequence.tolist()
```

``` python
# Define the maximum list
max_length = 4034

# Encoding by calling the function get_encoding()
one_hot_train_fungi = get_encoding(x_train_fungi, max_length)
one_hot_val_fungi = get_encoding(x_val_fungi, max_length)
one_hot_test_fungi = get_encoding(x_test_fungi, max_length)
```

``` python
# View the encoding results
input_train_fungi.info()
```

    ## <class 'pandas.core.frame.DataFrame'>
    ## RangeIndex: 118 entries, 0 to 117
    ## Data columns (total 1 columns):
    ## Sequence    118 non-null object
    ## dtypes: object(1)
    ## memory usage: 1.0+ KB

``` python
print(one_hot_test_fungi[1:2, :20, :20])
```

    ## [[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]]]

#### Change the label into list data format

``` python
# Change the data into 
y_train_fungi = label_train_fungi.label.tolist()
y_val_fungi = label_val_fungi.label.tolist()
y_test_fungi = label_test_fungi.label.tolist()
```

``` python
# View the label data

print(y_train_fungi)
```

    ## [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#### Save all of the fungi all of the data encoded

``` python
# Save the input data
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_train_fungi.npy', one_hot_train_fungi)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_val_fungi.npy', one_hot_val_fungi)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_test_fungi.npy', one_hot_test_fungi)

# Save the label data 
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_train_fungi.npy', y_train_fungi)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_val_fungi.npy', y_val_fungi)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_test_fungi.npy', y_test_fungi)
```

### Bacteria

``` python
# bacteria
training_bacteria = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/bacteria_training.csv", index_col = False)
validation_bacteria = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/bacteria_validation.csv", index_col = False)
testing_bacteria = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/bacteria_testing.csv", index_col = False)
```

``` python
training_bacteria.head(2)
```

    ##                                             Sequence  label
    ## 0  MGNICGTSGSHYVYSPPVSPRHVSGSSTPVHSVGGQGLTSVYQLSA...      1
    ## 1  MGNVCVGGSRMSHQVYSPDRADTPPRSERNTPDRRQRAAGDAERTQ...      1

``` python
# Define the input and the label of data 

# Training datasets
input_train_bacteria = training_bacteria[["Sequence"]]
label_train_bacteria = training_bacteria[["label"]]

# Validation datasets
input_val_bacteria = validation_bacteria[["Sequence"]]
label_val_bacteria = validation_bacteria[["label"]]

# Testing data 
input_test_bacteria = testing_bacteria[["Sequence"]]
label_test_bacteria = testing_bacteria[["label"]]
```

``` python
# To get the information about the length of the data

from collections import Counter
field_length_train_bacteria = input_train_bacteria.Sequence.astype(str).map(len) 
field_length_val_bacteria = input_val_bacteria.Sequence.astype(str).map(len)
field_length_test_bacteria = input_test_bacteria.Sequence.astype(str).map(len) 

print(max(field_length_train_bacteria)) 
```

    ## 2338

``` python
print(max(field_length_val_bacteria)) 
```

    ## 2574

``` python
print(max(field_length_test_bacteria))
```

    ## 1835

``` python
# Change the data to list
x_train_bacteria = input_train_bacteria.Sequence.tolist()
x_val_bacteria = input_val_bacteria.Sequence.tolist()
x_test_bacteria = input_test_bacteria.Sequence.tolist()
```

``` python
# Define the maximum list
max_length = 2574

# Encoding by calling the function get_encoding()
one_hot_train_bacteria = get_encoding(x_train_bacteria, max_length)
one_hot_val_bacteria = get_encoding(x_val_bacteria, max_length)
one_hot_test_bacteria = get_encoding(x_test_bacteria, max_length)
```

``` python
# View the encoding results
input_train_bacteria.info()
```

    ## <class 'pandas.core.frame.DataFrame'>
    ## RangeIndex: 228 entries, 0 to 227
    ## Data columns (total 1 columns):
    ## Sequence    228 non-null object
    ## dtypes: object(1)
    ## memory usage: 1.9+ KB

``` python
print(one_hot_test_bacteria[1:2, :20, :20])
```

    ## [[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
    ##   [0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
    ##   [0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]]

#### Change the label into list data format

``` python
# Change the data into 
y_train_bacteria = label_train_bacteria.label.tolist()
y_val_bacteria = label_val_bacteria.label.tolist()
y_test_bacteria = label_test_bacteria.label.tolist()
```

``` python
# View the label data

print(len(y_train_bacteria))
```

    ## 228

#### Save all of the bacteria all of the data encoded

``` python
# Save the input data
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_train_bacteria.npy', one_hot_train_bacteria)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_val_bacteria.npy', one_hot_val_bacteria)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_test_bacteria.npy', one_hot_test_bacteria)

# Save the label data 
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_train_bacteria.npy', y_train_bacteria)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_val_bacteria.npy', y_val_bacteria)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_test_bacteria.npy', y_test_bacteria)
```

### Data all together

``` python
# all
training_all = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/all_training.csv", index_col = False)
validation_all = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/all_validation.csv", index_col = False)
testing_all = pd.read_csv("../../../data/secreted_data/ready_to_process/splitted-data/all_testing.csv", index_col = False)
```

``` python
training_all.head(2)
```

    ##                                             Sequence  label
    ## 0  MQSVLLLTVLTQSFIATASPLVERSTPLSFAEKRPQKVSYDWTTPY...      1
    ## 1  MRLANFLFYLAPMIVSSLAFDFVPLSGELDFSQEMVFINLTQQQFS...      1

``` python
# Define the input and the label of data 

# Training datasets
input_train_all = training_all[["Sequence"]]
label_train_all = training_all[["label"]]

# Validation datasets
input_val_all = validation_all[["Sequence"]]
label_val_all = validation_all[["label"]]

# Testing data 
input_test_all = testing_all[["Sequence"]]
label_test_all = testing_all[["label"]]
```

``` python
# To get the information about the length of the data

from collections import Counter
field_length_train_all = input_train_all.Sequence.astype(str).map(len) 
field_length_val_all = input_val_all.Sequence.astype(str).map(len)
field_length_test_all = input_test_all.Sequence.astype(str).map(len) 

print(max(field_length_train_all)) 
```

    ## 1959

``` python
print(max(field_length_val_all)) 
```

    ## 4034

``` python
print(max(field_length_test_all))
```

    ## 2574

``` python
# Change the data to list
x_train_all = input_train_all.Sequence.tolist()
x_val_all = input_val_all.Sequence.tolist()
x_test_all = input_test_all.Sequence.tolist()
```

``` python
# Define the maximum list
max_length = 4034

# Encoding by calling the function get_encoding()
one_hot_train_all = get_encoding(x_train_all, max_length)
one_hot_val_all = get_encoding(x_val_all, max_length)
one_hot_test_all = get_encoding(x_test_all, max_length)
```

``` python
# View the encoding results
input_train_all.info()
```

    ## <class 'pandas.core.frame.DataFrame'>
    ## RangeIndex: 448 entries, 0 to 447
    ## Data columns (total 1 columns):
    ## Sequence    448 non-null object
    ## dtypes: object(1)
    ## memory usage: 3.6+ KB

``` python
print(one_hot_test_all[1:2, :20, :20])
```

    ## [[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]]]

#### Change the label into list data format

``` python
# Change the data into 
y_train_all = label_train_all.label.tolist()
y_val_all = label_val_all.label.tolist()
y_test_all = label_test_all.label.tolist()
```

``` python
# View the label data

print(len(y_train_all))
```

    ## 448

#### Save all of the all all of the data encoded

``` python
# Save the input data
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_train_all.npy', one_hot_train_all)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_val_all.npy', one_hot_val_all)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/x_test_all.npy', one_hot_test_all)

# Save the label data 
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_train_all.npy', y_train_all)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_val_all.npy', y_val_all)
np.save('../../../data/secreted_data/ready_to_process/encoded_files/y_test_all.npy', y_test_all)
```
