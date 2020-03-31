Getting the effector data in Multiclass based on Domain
=======================================================

Introduction
------------

In the previous result of predicting effector and non-effector protein,
although it is obvious that the first 50 sequence are the most important
sequences on the prediction process, however, there is no certainly
exact position detected in all datasets.

Aims
----

### Questions

By limiting or group the effector data into several group of pathogens
and pests, will the pattern of feature importance will be obvious?

In order to get better understanding in how the models recognizing the
feature importance in effector prediction, I will group the effector and
non-effector data, into several groups, which are:

-   bacteria
-   fungi
-   oomycetes
-   nematodes
-   insects

Executions
----------

### Load all of the library needed

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(taxize)
library(caret)
```

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
reticulate::use_condaenv(condaenv = "tensorflow2", conda = "/anaconda3/bin/conda")
```

``` python
sys.executable
```

    ## '/anaconda3/envs/tensorflow2/bin/python'

### Load the effector and non-effector data

``` r
effector_data <- data.table::fread("../../../../data/getting-data-new/binary-class-data/effector_data.csv")
non_effector_data <- data.table::fread("../../../../data/getting-data-new/binary-class-data/non_effector_data.csv")
```

``` r
# Check the count of effector and non-effector data in total

eff_count <- effector_data %>%  nrow()
non_eff_count <- non_effector_data %>% nrow()

df_count <- data.frame("Effector" = eff_count, "Non-effector" = non_eff_count)

df_count %>% 
  knitr::kable()
```

|  Effector|  Non.effector|
|---------:|-------------:|
|       402|           398|

### Get the list of the PathogenID with the class

Using the function in `kingdom.R`, we can retrieve the set of list of
effector and non-effector data of ID and the class.

``` r
effector_data_multiclass <- readRDS("../../../../data/getting-data-new/multi-class-data/class_df_effectors.rds")
non_effector_data_multiclass <- readRDS("../../../../data/getting-data-new/multi-class-data/class_df_noneffectors.rds")
```

``` r
effector_data_multiclass %>% 
  group_by(class) %>% 
  summarise(count_id_uniq = n())
```

    ## # A tibble: 5 x 2
    ##   class    count_id_uniq
    ##   <fct>            <int>
    ## 1 bacteria            52
    ## 2 fungi               44
    ## 3 insecta              1
    ## 4 nematoda             2
    ## 5 oomycete            36

``` r
non_effector_data_multiclass %>% 
  group_by(class) %>% 
  summarise(count_id_uniq = n())
```

    ## # A tibble: 4 x 2
    ##   class    count_id_uniq
    ##   <fct>            <int>
    ## 1 bacteria            52
    ## 2 fungi               42
    ## 3 nematoda             2
    ## 4 oomycete            36

Getting the class of each datasets
----------------------------------

``` r
eff_all_with_class <- effector_data %>% 
  left_join(effector_data_multiclass %>% 
              ungroup() %>% 
              mutate("PathogenID" = as.integer(PathogenID)), by = "PathogenID")
```

``` r
non_eff_all_with_class <- non_effector_data %>% 
  left_join(non_effector_data_multiclass %>% 
              ungroup %>% 
              mutate("PathogenID" = as.integer(PathogenID)), by = c("txid" = "PathogenID"))
```

### Getting more information about each class

``` r
# Count the classes in effector data
eff_all_with_class %>% 
  group_by(class) %>% 
  dplyr::summarise(count_class = n())
```

    ## # A tibble: 5 x 2
    ##   class    count_class
    ##   <fct>          <int>
    ## 1 bacteria         190
    ## 2 fungi            113
    ## 3 insecta            2
    ## 4 nematoda           2
    ## 5 oomycete          95

``` r
# Count the classes in non-effector data
non_eff_all_with_class %>% 
  group_by(class) %>% 
  dplyr::summarise(count_class = n())
```

    ## # A tibble: 4 x 2
    ##   class    count_class
    ##   <fct>          <int>
    ## 1 bacteria         190
    ## 2 fungi            111
    ## 3 nematoda           2
    ## 4 oomycete          95

Considering how small teh number of `nematoda` and `insecta` data in
both datasets (will not be enough for training the model), then we will
only include `bacteria`, `fungi`, and `oomycete`.

``` r
eff_all_with_class <- eff_all_with_class %>% 
  ungroup() %>% 
  dplyr::filter(!(class %in% c("nematoda", "insecta")))

non_eff_all_with_class <- non_eff_all_with_class %>% 
  ungroup() %>% 
  dplyr::filter(!(class %in% c("nematoda", "insecta")))
```

### Split the data into training, development, and test datasets

As alternative of scikit-learn that we use in Python, we can use caret
package.

Before splittig data into three categories, we need to label the data

``` r
eff_all_with_class %>% 
  group_by(class) %>% 
  summarise(count = n())
```

    ## # A tibble: 3 x 2
    ##   class    count
    ##   <fct>    <int>
    ## 1 bacteria   190
    ## 2 fungi      113
    ## 3 oomycete    95

Multi-class processing data
---------------------------

Labeling the effector and non-effector data:

``` r
eff_all_with_class <- eff_all_with_class %>% 
  dplyr::mutate(data = 1) %>% 
  dplyr::select(Sequence, class, data) %>% 
  `colnames<-`(c("sequence", "class", "data"))
```

``` r
non_eff_all_with_class <- non_eff_all_with_class %>% 
  dplyr::mutate(data = 0) %>% 
  dplyr::select(sequence, class, data) %>% 
  `colnames<-`(c("sequence", "class", "data"))
```

``` r
# Combine all of the datasets

data_all <- eff_all_with_class %>% 
  rbind(non_eff_all_with_class)
```

``` r
# Labelling the data

data_with_label <- data_all %>%
  dplyr::mutate(label = case_when(
  (class == "bacteria" & data == 1) ~ 1,
  (class == "bacteria" & data == 0) ~ 2,
  (class == "fungi" & data == 1) ~ 3, 
  (class == "fungi" & data == 0) ~ 4, 
  (class == "oomycete" & data == 1) ~ 5, 
  TRUE ~ 6))
```

#### Combine both data into one, then we can split them

``` r
data_to_split <- data_with_label
```

``` r
# set.seed(100)  # For reproducibility

# Create index for testing and training data
training_id <- createDataPartition(y = data_to_split$label, 
                               p = 0.6, list = FALSE)

# subset power_plant data to training
training <- data_to_split[training_id,]


# subset the rest to test
rest <- data_to_split[-training_id,]

# Splitting the rest into two different class of datasets

val_id <- createDataPartition(y = rest$label, 
                               p = 0.5, list = FALSE)

validation <- rest[val_id,]

testing <- rest[-val_id,]
```

``` r
training %>% 
  group_by(label) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
```

    ## # A tibble: 6 x 3
    ##   label count  freq
    ##   <dbl> <int> <dbl>
    ## 1     1   119 0.249
    ## 2     2   109 0.229
    ## 3     3    68 0.143
    ## 4     4    67 0.140
    ## 5     5    60 0.126
    ## 6     6    54 0.113

``` r
validation %>% 
  group_by(label) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
```

    ## # A tibble: 6 x 3
    ##   label count   freq
    ##   <dbl> <int>  <dbl>
    ## 1     1    35 0.220 
    ## 2     2    41 0.258 
    ## 3     3    23 0.145 
    ## 4     4    22 0.138 
    ## 5     5    13 0.0818
    ## 6     6    25 0.157

``` r
testing %>% 
  group_by(label) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
```

    ## # A tibble: 6 x 3
    ##   label count  freq
    ##   <dbl> <int> <dbl>
    ## 1     1    36 0.228
    ## 2     2    40 0.253
    ## 3     3    22 0.139
    ## 4     4    22 0.139
    ## 5     5    22 0.139
    ## 6     6    16 0.101

``` r
# Save the dataframe to CSV
write_csv(training, "multi-class-data/data-sets/new_training_data.csv")
write_csv(validation, "multi-class-data/data-sets/new_validation_data.csv")
write_csv(testing, "multi-class-data/data-sets/new_testing_data.csv")
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

    ## /Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/r-scripts/getting-data-new/multi-class

### Load the data

``` python
training_data = pd.read_csv('../../../../data/getting-data-new/multi-class-data/data-sets/training.csv', index_col = False)
validation_data = pd.read_csv('../../../../data/getting-data-new/multi-class-data/data-sets/validation.csv', index_col = False)
testing_data = pd.read_csv('../../../../data/getting-data-new/multi-class-data/data-sets/testing.csv', index_col = False)
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
from collections import Counter
field_length_train = input_train.sequence.astype(str).map(len) 
field_length_val = input_val.sequence.astype(str).map(len)
field_length_test = input_test.sequence.astype(str).map(len) 

print(max(field_length_train)) 
```

    ## 2758

``` python
print(max(field_length_val)) 
```

    ## 1795

``` python
print(max(field_length_test))
```

    ## 4034

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
```

    ## <class 'pandas.core.frame.DataFrame'>
    ## RangeIndex: 477 entries, 0 to 476
    ## Data columns (total 1 columns):
    ## sequence    477 non-null object
    ## dtypes: object(1)
    ## memory usage: 3.9+ KB

``` python
print(one_hot_test[1:2, :20, :20])
```

    ## [[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]]]

#### Change the label into list data format

``` python
# Import library to change the data to encoded categorical data
from keras.utils.np_utils import to_categorical
```

``` python
# Change the data into list
y_train = to_categorical(label_train.label.tolist())
y_val = to_categorical(label_val.label.tolist())
y_test = to_categorical(label_test.label.tolist())
```

#### Save all of the data

``` python
# Save the input data
np.save('multi-class-data/data-sets/x_train.npy', one_hot_train)
np.save('multi-class-data/data-sets/x_val.npy', one_hot_val)
np.save('multi-class-data/data-sets/x_test.npy', one_hot_test)

# Save the label data 
np.save('multi-class-data/data-sets/y_train.npy', y_train)
np.save('multi-class-data/data-sets/y_val.npy', y_val)
np.save('multi-class-data/data-sets/y_test.npy', y_test)
```
