Simple Fully Dense Model for Effector and Non-effector Protein
==============================================================

In this report, the fully conected model will be developed for the
encoded effector and non-effector data.

Introduction
------------

Since new encoding method was done, we will use two kind of data on the
same models.

Load the python library
-----------------------

``` r
library(keras)
```

``` python
import numpy as np
from keras import models
from keras import layers
from keras import losses
from keras import metrics
from keras import optimizers
from keras.layers import Dense, Dropout
```

Load the encoded data
---------------------

There are two kind of data that will be used for each model, which are
encoded by integer encoding and also onehot encoding.

``` python
x_train_integer = np.load('../../../data/0001-encoded-data/integer/x_train.npy')
x_dev_integer = np.load('../../../data/0001-encoded-data/integer/x_val.npy')
x_test_integer = np.load('../../../data/0001-encoded-data/integer/x_test.npy')

x_train_onehot = np.load('../../../data/0001-encoded-data/one-hot/x_train.npy')
x_dev_onehot = np.load('../../../data/0001-encoded-data/one-hot/x_val.npy')
x_test_onehot = np.load('../../../data/0001-encoded-data/one-hot/x_test.npy')

y_train = np.load('../../../data/0001-encoded-data/y_train.npy')
y_dev = np.load('../../../data/0001-encoded-data/y_val.npy')
y_test = np.load('../../../data/0001-encoded-data/y_test.npy')
```

We can print the shape of the one hot encoded data

``` python
print(x_train_onehot.shape)
print(x_dev_onehot.shape)
print(x_test_onehot.shape)
```

We need to reshape the data in order so it can be fed to fully connected
model.

``` python
x_train_onehot = x_train_onehot.reshape(576, 50000)
x_dev_onehot = x_dev_onehot.reshape(193, 50000)
x_test_onehot = x_test_onehot.reshape(193, 50000)
```

Fully connected model with integer encoded data
-----------------------------------------------

``` python
model_first = models.Sequential()
model_first.add(Dense(3, activation = 'relu', input_dim = 2500))
model_first.add(Dense(1, activation='sigmoid'))
model_first.compile(optimizer = 'sgd',
                    loss = 'binary_crossentropy',
                    metrics = ['accuracy'])
history = model_first.fit(x = x_train_integer, 
                          y = y_train, 
                          epochs = 30, 
                          batch_size = 16, 
                          validation_data = (x_dev_integer, y_dev), 
                          verbose = 0)
```

### Visualizing the result

``` python
import matplotlib.pyplot as plt

history_dict = history.history
history_dict.keys()
['val_loss', 'val_acc', 'loss', 'acc']
```

``` python
plt.clf()
loss_values = history_dict['loss']
val_loss_values = history_dict['val_loss']
epochs = range(1, len(loss_values) + 1)
print(epochs)

plt.plot(epochs, loss_values, 'bo', label="Training loss")
plt.plot(epochs, val_loss_values, 'b', label="Validation loss")
plt.title('Training and validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

plt.show()
```

``` python
plt.clf()
acc_values = history_dict['acc']
val_acc_values = history_dict['val_acc']
epochs = range(1, len(acc_values) + 1)
print(epochs)


plt.plot(epochs, acc_values, 'go', label="Training Acc")
plt.plot(epochs, val_acc_values, 'ro', label="Validation Acc")
plt.title('Training and validation accuracy')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

plt.show()
```

### Predict the test data

``` python
#  retrain the model
from keras.layers import Dense, Dropout
model_first = models.Sequential()
model_first.add(Dense(3, activation = 'relu', input_dim = 2500))
model_first.add(Dense(1, activation='sigmoid'))

model_first.compile(optimizer = 'sgd',
                    loss = 'binary_crossentropy',
                    metrics = ['accuracy'])
                    
history = model_first.fit(x = x_train_integer, 
                          y = y_train, 
                          epochs = 7, 
                          batch_size = 16,
                          verbose = 0)

#  evaluate using the test data set
results = model_first.evaluate(x_test_integer, y_test)
print(results)
```

Fully connected model with One hot encoded data
-----------------------------------------------

``` python
from keras.layers import Dense, Dropout

model_first = models.Sequential()
model_first.add(Dense(3, activation = 'relu', input_dim = 50000))
model_first.add(Dense(1, activation='sigmoid'))
model_first.compile(optimizer = 'sgd',
                    loss = 'binary_crossentropy',
                    metrics = ['accuracy'])
history_onehot = model_first.fit(x = x_train_onehot, 
                          y = y_train, 
                          epochs = 30, 
                          batch_size = 32, 
                          validation_data = (x_dev_onehot, y_dev), 
                          verbose = 0)
```

### Visualizing the results

``` python
import matplotlib.pyplot as plt

history_dict_onehot = history_onehot.history
history_dict_onehot.keys()
['val_loss', 'val_acc', 'loss', 'acc']
```

``` python
plt.clf()
loss_values = history_dict_onehot['loss']
val_loss_values = history_dict_onehot['val_loss']
epochs = range(1, len(loss_values) + 1)
print(epochs)

plt.plot(epochs, loss_values, 'bo', label="Training loss")
plt.plot(epochs, val_loss_values, 'b', label="Validation loss")
plt.title('Training and validation loss (One hot encoding data)')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

plt.show()
```

``` python
plt.clf()
acc_values = history_dict_onehot['acc']
val_acc_values = history_dict_onehot['val_acc']
epochs = range(1, len(acc_values) + 1)
print(epochs)


plt.plot(epochs, acc_values, 'go', label="Training Acc")
plt.plot(epochs, val_acc_values, 'g', label="Validation Acc")
plt.title('Training and validation accuracy (One hot encoding data)')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

plt.show()
```

### Predict the test data

``` python
from keras.layers import Dense, Dropout

model_first = models.Sequential()
model_first.add(Dense(3, activation = 'relu', input_dim = 50000))
model_first.add(Dense(1, activation='sigmoid'))

model_first.compile(optimizer = 'sgd',
                    loss = 'binary_crossentropy',
                    metrics = ['accuracy'])
                    
history_onehot = model_first.fit(x = x_train_onehot, 
                          y = y_train, 
                          epochs = 4, 
                          batch_size = 32,
                          verbose = 0)
                          
results = model_first.evaluate(x_test_onehot, y_test)
print(results)                  
```
