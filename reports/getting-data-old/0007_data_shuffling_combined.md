Shuffling data
==============

This report will present the part of data processing where the data will
be processed and shufllled. Before doing the data reprocessing, we need
to load all libraries that we need to process the data.

``` r
library(keras)
```

``` python
import pandas as pd
import numpy as np
import random
```

``` python
import os
cwd = os.getcwd()
print(cwd)
```

    ## /Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/r-scripts/getting-data-old

Data label
----------

Since two different encoded data will use the same data label, then we
will only need to load the label once.

``` python
y_train = np.load('../../../data/getting-data-old/0001-encoded-data/y_train.npy')
y_val = np.load('../../../data/getting-data-old/0001-encoded-data/y_val.npy')
y_test = np.load('../../../data/getting-data-old/0001-encoded-data/y_test.npy')
```

Data with Integer Encoding
--------------------------

We need to load all of the data

``` python
x_train_integer = np.load('../../../data/getting-data-old/0001-encoded-data/integer/x_train.npy')
x_val_integer = np.load('../../../data/getting-data-old/0001-encoded-data/integer/x_val.npy')
x_test_integer = np.load('../../../data/getting-data-old/0001-encoded-data/integer/x_test.npy')
```

``` python
# This Section is used to shuffle the data

# aggregates elements
data_training_integer = list(zip(x_train_integer, y_train))
data_development_integer = list(zip(x_val_integer, y_val))
data_testing_integer = list(zip(x_test_integer, y_test))

# shuffle the aggragated element on the list
random.shuffle(data_training_integer)
random.shuffle(data_development_integer)
random.shuffle(data_testing_integer)

# combine data training dan data development become one list of data train

data_train = data_training_integer + data_development_integer

# split the shuffled data
x_train_integer, y_train_integer = zip(*data_train)
x_test_integer, y_test_integer = zip(*data_testing_integer)

# unpack the tuples
x_train_integer = np.array(list(x_train_integer))
y_train_integer = np.array(list(y_train_integer))
x_test_integer = np.array(list(x_test_integer))
y_test_integer = np.array(list(y_test_integer))
```

Save the shuffled and combined (training and development) data that are ready to be used using Scikit-learn package
===================================================================================================================

``` python
# np.save('../../../data/getting-data-old/0001-encoded-data/integer/x_train_integer.npy', x_train_integer)
# np.save('../../../data/getting-data-old/0001-encoded-data/integer/y_train_integer.npy', y_train_integer)
# np.save('../../../data/getting-data-old/0001-encoded-data/integer/x_test_integer.npy', x_test_integer)
# np.save('../../../data/getting-data-old/0001-encoded-data/integer/y_test_integer.npy', y_test_integer)
```

We can check again if the size of the input and the label matches

``` python
print(x_test_integer.shape)
```

    ## (193, 2500)

``` python
print(y_test_integer.shape)
```

    ## (193,)

Data with Onehot Encoding
-------------------------

We need to load all of the data

``` python
x_train_onehot = np.load('../../../data/getting-data-old/0001-encoded-data/one-hot/x_train.npy')
x_val_onehot = np.load('../../../data/getting-data-old/0001-encoded-data/one-hot/x_val.npy')
x_test_onehot = np.load('../../../data/getting-data-old/0001-encoded-data/one-hot/x_test.npy')
```

``` python
# This Section is used to shuffle the data

# aggregates elements
data_training_onehot = list(zip(x_train_onehot, y_train))
data_development_onehot = list(zip(x_val_onehot, y_val))
data_testing_onehot = list(zip(x_test_onehot, y_test))

# shuffle the aggragated element on the list
random.shuffle(data_training_onehot)
random.shuffle(data_development_onehot)
random.shuffle(data_testing_onehot)

# combine data training dan data development become one list of data train

data_train = data_training_onehot + data_development_onehot

# split the shuffled data
x_train_onehot, y_train_onehot = zip(*data_train)
x_test_onehot, y_test_onehot = zip(*data_testing_onehot)

# unpack the tuples
x_train_onehot = np.array(list(x_train_onehot))
y_train_onehot = np.array(list(y_train_onehot))
x_test_onehot = np.array(list(x_test_onehot))
y_test_onehot = np.array(list(y_test_onehot))
```

Save the shuffled and combined (training and development) data that are ready to be used using Scikit-learn package
===================================================================================================================

``` python
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/x_train_onehot.npy', x_train_onehot)
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/y_train_onehot.npy', y_train_onehot)
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/x_test_onehot.npy', x_test_onehot)
# np.save('../../../data/getting-data-old/0001-encoded-data/one-hot/y_test_onehot.npy', y_test_onehot)
```

We can check again if the size of the input and the label matches

``` python
print(x_train_onehot.shape)
```

    ## (769, 2500, 20)

``` python
print(y_train_onehot.shape)
```

    ## (769,)
