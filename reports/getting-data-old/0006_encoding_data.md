Encoding the amino acids protein
================================

In this report, both ways to encode amino acids protein for the effector
and non-effector data will be done. The first encoding method will be
one hot encoding, the second method will be the method where we assign
natural numbers to each amino acids.

``` r
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

    ## /Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/r-scripts/getting-data-old

Load the data
-------------

``` python
input_train = pd.read_csv('../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-sequence.csv', header = None, index_col = False, names = ["sequence"])
input_val = pd.read_csv('../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-sequence.csv', header = None, index_col = False, names = ["sequence"])
input_test = pd.read_csv('../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-sequence.csv', header = None, index_col = False, names = ["sequence"])
label_train = pd.read_csv('../../../data/getting-data-old/BLAST-data/0003-new-data-sets/training-label.csv', header = None, index_col = False, names = ["label"])
label_val = pd.read_csv('../../../data/getting-data-old/BLAST-data/0003-new-data-sets/validation-label.csv', header = None, index_col = False, names = ["label"])
label_test = pd.read_csv('../../../data/getting-data-old/BLAST-data/0003-new-data-sets/testing-label.csv', header = None, index_col = False, names = ["label"])
```

``` python
print(label_train.info())
```

    ## <class 'pandas.core.frame.DataFrame'>
    ## RangeIndex: 578 entries, 0 to 577
    ## Data columns (total 1 columns):
    ## label    578 non-null int64
    ## dtypes: int64(1)
    ## memory usage: 4.6 KB
    ## None

``` python
from collections import Counter
field_length_train = input_train.sequence.astype(str).map(len) 
field_length_val = input_val.sequence.astype(str).map(len)
field_length_test = input_test.sequence.astype(str).map(len) 
```

``` python
print(field_length_train)
# from collections import Counter
# print(Counter([field_length_train > 2500]))
```

    ## 0      2039
    ## 1      1048
    ## 2       350
    ## 3      1047
    ## 4       765
    ## 5       560
    ## 6        97
    ## 7       419
    ## 8       329
    ## 9       188
    ## 10     2349
    ## 11     2369
    ## 12      649
    ## 13      384
    ## 14     2108
    ## 15      432
    ## 16     2349
    ## 17      301
    ## 18      237
    ## 19      164
    ## 20      231
    ## 21     2185
    ## 22      117
    ## 23      251
    ## 24      820
    ## 25      216
    ## 26     2124
    ## 27      248
    ## 28     1785
    ## 29      961
    ##        ... 
    ## 548    2133
    ## 549    2058
    ## 550    1921
    ## 551    1798
    ## 552    1590
    ## 553    1575
    ## 554    1520
    ## 555    1679
    ## 556    2367
    ## 557    1388
    ## 558    1767
    ## 559    1910
    ## 560    1411
    ## 561    1411
    ## 562    1411
    ## 563    1411
    ## 564    1411
    ## 565    1411
    ## 566    1106
    ## 567     723
    ## 568    1593
    ## 569    1293
    ## 570    1107
    ## 571    1120
    ## 572    1115
    ## 573    1115
    ## 574    1778
    ## 575    1411
    ## 576    1334
    ## 577    1388
    ## Name: sequence, Length: 578, dtype: int64

``` python
import matplotlib.pyplot as plt
plt.clf()
plt.hist(field_length_train, bins=100)
```

    ## (array([18., 29., 28., 18., 26., 18., 19., 14., 10.,  7., 13.,  5.,  9.,
    ##         8., 10.,  4., 10.,  3.,  6.,  2.,  0.,  3.,  4.,  2., 11.,  4.,
    ##         7.,  5.,  2.,  1.,  2.,  2.,  1.,  3.,  7.,  1.,  8.,  5.,  6.,
    ##         2.,  3.,  1.,  2.,  4.,  4.,  7.,  4.,  2.,  6., 15., 11., 19.,
    ##        22., 14.,  8., 14.,  9., 25., 44., 13., 10.,  6.,  0.,  1.,  0.,
    ##         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
    ##         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
    ##         0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.]), array([  58.  ,   97.76,  137.52,  177.28,  217.04,  256.8 ,  296.56,
    ##         336.32,  376.08,  415.84,  455.6 ,  495.36,  535.12,  574.88,
    ##         614.64,  654.4 ,  694.16,  733.92,  773.68,  813.44,  853.2 ,
    ##         892.96,  932.72,  972.48, 1012.24, 1052.  , 1091.76, 1131.52,
    ##        1171.28, 1211.04, 1250.8 , 1290.56, 1330.32, 1370.08, 1409.84,
    ##        1449.6 , 1489.36, 1529.12, 1568.88, 1608.64, 1648.4 , 1688.16,
    ##        1727.92, 1767.68, 1807.44, 1847.2 , 1886.96, 1926.72, 1966.48,
    ##        2006.24, 2046.  , 2085.76, 2125.52, 2165.28, 2205.04, 2244.8 ,
    ##        2284.56, 2324.32, 2364.08, 2403.84, 2443.6 , 2483.36, 2523.12,
    ##        2562.88, 2602.64, 2642.4 , 2682.16, 2721.92, 2761.68, 2801.44,
    ##        2841.2 , 2880.96, 2920.72, 2960.48, 3000.24, 3040.  , 3079.76,
    ##        3119.52, 3159.28, 3199.04, 3238.8 , 3278.56, 3318.32, 3358.08,
    ##        3397.84, 3437.6 , 3477.36, 3517.12, 3556.88, 3596.64, 3636.4 ,
    ##        3676.16, 3715.92, 3755.68, 3795.44, 3835.2 , 3874.96, 3914.72,
    ##        3954.48, 3994.24, 4034.  ]), <a list of 100 Patch objects>)

``` python
plt.ylabel('Count')
plt.xlabel('Length')
plt.title('Histogram of Length of Sequences')
plt.show()
```

<img src="/Users/kristian/Documents/Workspace/ruth-effectors-prediction/reports/getting-data-old/0006_encoding_data_files/figure-markdown_github/unnamed-chunk-8-1.png" width="672" />

First encoding method - one hot encoding
----------------------------------------

``` python
def get_key(mydict, element):
    key = list(mydict.keys())[list(mydict.values()).index(element)]
    return(key)

amino = ['R', 'K', 'D', 'E', 'Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'W', 'A', 'I', 'L', 'M', 'F', 'V', 'P', 'G']
token_index = dict(zip(range(1, (len(amino)+1)), amino))

max_length = 2500
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
x_train = input_train[input_train['sequence'].map(len) < 2500].sequence.tolist()
x_val = input_val.sequence.tolist()
x_test = input_test.sequence.tolist()
one_hot_train = get_encoding(x_train, max_length)
one_hot_val = get_encoding(x_val, max_length)
one_hot_test = get_encoding(x_test, max_length)
```

``` python
print(input_train.loc[input_train.sequence.astype(str).map(len) > 2500])
```

    ##                                               sequence
    ## 105  MPSRMGYSRISSGLNASRGASPAPQPDTPPQTPPPDNRRRTRMDSP...
    ## 305  MRDEMWNTATEPIAIIGSGCKFPGGSTTPSKLWELLKDPKDIVSEI...

Since we only select the training data sequence with length less or
equal to 2500, then we also need to drop the data on that index from the
the label training.

``` python
drop_idx = [105, 305]
label_train = label_train.loc[~label_train.index.isin(drop_idx)]
print(label_train)
```

    ##      label
    ## 0        0
    ## 1        0
    ## 2        1
    ## 3        0
    ## 4        1
    ## 5        0
    ## 6        1
    ## 7        1
    ## 8        1
    ## 9        1
    ## 10       0
    ## 11       0
    ## 12       1
    ## 13       1
    ## 14       0
    ## 15       1
    ## 16       0
    ## 17       1
    ## 18       1
    ## 19       1
    ## 20       1
    ## 21       0
    ## 22       1
    ## 23       1
    ## 24       1
    ## 25       1
    ## 26       0
    ## 27       1
    ## 28       0
    ## 29       0
    ## ..     ...
    ## 548      0
    ## 549      0
    ## 550      0
    ## 551      0
    ## 552      0
    ## 553      0
    ## 554      0
    ## 555      0
    ## 556      0
    ## 557      1
    ## 558      1
    ## 559      1
    ## 560      1
    ## 561      1
    ## 562      1
    ## 563      1
    ## 564      1
    ## 565      1
    ## 566      1
    ## 567      1
    ## 568      1
    ## 569      1
    ## 570      1
    ## 571      1
    ## 572      1
    ## 573      1
    ## 574      1
    ## 575      1
    ## 576      1
    ## 577      1
    ## 
    ## [576 rows x 1 columns]

Second encoding method - assign natural numbers to each amino acids
-------------------------------------------------------------------

``` python
def get_value(mydict, element):
    key = mydict.get(element)
    return(key)
    
dic = {'A':1, 'B':22, 'U':23,'J':24,'Z':25,'O':26,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20,'X':21}
max_len = 2500

def encoding_each_sequence(input):
    results = np.zeros((len(input), max_len))
    for i, sample in enumerate(input):
        for j, character in enumerate(sample):
            results[i, j] = get_value(dic, character)
    return results
```

``` python
natural_number_assign_train = encoding_each_sequence(x_train)
natural_number_assign_val = encoding_each_sequence(x_val)
natural_number_assign_test = encoding_each_sequence(x_test)
```

``` python
print(natural_number_assign_val)
```

    ## [[11. 10. 10. ...  0.  0.  0.]
    ##  [11.  3. 10. ...  0.  0.  0.]
    ##  [11. 16.  3. ...  0.  0.  0.]
    ##  ...
    ##  [11.  6. 12. ...  0.  0.  0.]
    ##  [11.  1.  6. ...  0.  0.  0.]
    ##  [11.  1. 10. ...  0.  0.  0.]]

Saving all of the results in the data folder
--------------------------------------------

### One hot encoding

``` python
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/x_train.npy', one_hot_train)
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/x_val.npy', one_hot_val)
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/x_test.npy', one_hot_test)
```

### Encoding using Natural number

``` python
# np.save('../../../data/getting-data-old/0001-encoded-data/assign-natural-number/x_train.npy', natural_number_assign_train)
# np.save('../../../data/getting-data-old/0001-encoded-data/assign-natural-number/x_val.npy', natural_number_assign_val)
# np.save('../../../data/getting-data-old/0001-encoded-data/assign-natural-number/x_test.npy', natural_number_assign_test)
```

### Saving the data label

``` python
y_train = label_train.label.tolist()
y_val = label_val.label.tolist()
y_test = label_test.label.tolist()
# np.save('../../../data/getting-data-old/0001-encoded-data/y_train.npy', y_train)
# np.save('../../../data/getting-data-old/0001-encoded-data/y_val.npy', y_val)
# np.save('../../../data/getting-data-old/0001-encoded-data/y_test.npy', y_test)
```
