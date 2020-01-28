Encoding Sequence data after Blast
==================================

Background
----------

After we blast the data to make sure there is no effector inside, we
blast them again against each dataset to make sure there will not be.

Execution
---------

### Load the libraries

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

    ## /Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/r-scripts/getting-secreted-data

### Load the data

#### Oomycete data

In this process, the proocess of encoding the oomycete data will be
shown

``` python
# Oomycete
training_oomycete = pd.read_csv("../../../data/secreted_data/split-blast/csv_files/secreted_oomycete_training.csv", index_col = False)
validation_oomycete = pd.read_csv("../../../data/secreted_data/split-blast/csv_files/secreted_oomycete_validation.csv", index_col = False)
testing_oomycete = pd.read_csv("../../../data/secreted_data/split-blast/csv_files/secreted_oomycete_testing.csv", index_col = False)
```

``` python
training_oomycete.head(2)
```

    ##     ProteinID                                           Sequence  label
    ## 0  A0A0M5K865  MVKLYCAVVGVAGSAFSVRVDESDTVDDLKDAIKAKKPNDFKDIDA...      1
    ## 1      A5YTY8  MRLAQVVVVIAASFLVATDALSTTNANQAKIIKGTSPGGHSPRLLR...      1

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
# To get the information about the data

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

    ## 934

``` python
print(max(field_length_test_oomycete))
```

    ## 637

``` python
import matplotlib.pyplot as plt
plt.clf()
plt.hist(field_length_train_oomycete, bins=100)
```

    ## (array([2., 1., 3., 1., 0., 1., 3., 5., 4., 6., 7., 2., 3., 3., 3., 1., 2.,
    ##        6., 1., 2., 2., 3., 0., 3., 3., 1., 0., 1., 0., 1., 2., 2., 1., 1.,
    ##        0., 2., 2., 3., 0., 0., 1., 1., 0., 0., 1., 0., 0., 1., 1., 1., 1.,
    ##        0., 0., 0., 2., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0.,
    ##        0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0.,
    ##        0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1.]), array([ 54.  ,  61.66,  69.32,  76.98,  84.64,  92.3 ,  99.96, 107.62,
    ##        115.28, 122.94, 130.6 , 138.26, 145.92, 153.58, 161.24, 168.9 ,
    ##        176.56, 184.22, 191.88, 199.54, 207.2 , 214.86, 222.52, 230.18,
    ##        237.84, 245.5 , 253.16, 260.82, 268.48, 276.14, 283.8 , 291.46,
    ##        299.12, 306.78, 314.44, 322.1 , 329.76, 337.42, 345.08, 352.74,
    ##        360.4 , 368.06, 375.72, 383.38, 391.04, 398.7 , 406.36, 414.02,
    ##        421.68, 429.34, 437.  , 444.66, 452.32, 459.98, 467.64, 475.3 ,
    ##        482.96, 490.62, 498.28, 505.94, 513.6 , 521.26, 528.92, 536.58,
    ##        544.24, 551.9 , 559.56, 567.22, 574.88, 582.54, 590.2 , 597.86,
    ##        605.52, 613.18, 620.84, 628.5 , 636.16, 643.82, 651.48, 659.14,
    ##        666.8 , 674.46, 682.12, 689.78, 697.44, 705.1 , 712.76, 720.42,
    ##        728.08, 735.74, 743.4 , 751.06, 758.72, 766.38, 774.04, 781.7 ,
    ##        789.36, 797.02, 804.68, 812.34, 820.  ]), <a list of 100 Patch objects>)

``` python
plt.ylabel('Count')
plt.xlabel('Length')
plt.title('Histogram of Length of Sequences')
plt.show()
```

<img src="/Users/kristian/Documents/Workspace/ruth-effectors-prediction/reports/getting-data-secreted/0009-encode-data_files/figure-markdown_github/unnamed-chunk-8-1.png" width="672" />

### One hot encoding

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

``` python
# Change the data to list
x_train_oomycete = input_train_oomycete.Sequence.tolist()
x_val_oomycete = input_val_oomycete.Sequence.tolist()
x_test_oomycete = input_test_oomycete.Sequence.tolist()
```

``` python
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
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
    ##   [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
    ##   [0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
    ##   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
    ##   [0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
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

#### Save all of the data

``` python
# Save the input data
np.save('../../../data/secreted_data/split-blast/encoded_files/x_train_oomycete.npy', one_hot_train_oomycete)
np.save('../../../data/secreted_data/split-blast/encoded_files/x_val_oomycete.npy', one_hot_val_oomycete)
np.save('../../../data/secreted_data/split-blast/encoded_files/x_test_oomycete.npy', one_hot_test_oomycete)

# Save the label data 
np.save('../../../data/secreted_data/split-blast/encoded_files/y_train_oomycete.npy', y_train_oomycete)
np.save('../../../data/secreted_data/split-blast/encoded_files/y_val_oomycete.npy', y_val_oomycete)
np.save('../../../data/secreted_data/split-blast/encoded_files/y_test_oomycete.npy', y_test_oomycete)
```
