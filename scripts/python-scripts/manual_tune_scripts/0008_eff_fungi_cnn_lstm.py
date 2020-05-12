#!/usr/bin/env python3

# Import package to scan hyperparameter
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint as sp_randint
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold

# Import package to reprocess the data
import numpy as np
import random
import pandas as pd

# Import keras item
import keras
from keras.layers import *
from keras.optimizers import *
from keras.applications import *
from keras.models import *
from keras.models import Model
from keras.layers import Input, Dense
from keras import regularizers

# Get all of the data and reprocess them

# Get the reprocessed data from .npy file
x_train = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted_non_identical/x_train_fungi.npy')
y_train = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted_non_identical/y_train_fungi.npy')
x_dev = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted_non_identical/x_val_fungi.npy')
y_dev = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted_non_identical/y_val_fungi.npy')
x_test = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted_non_identical/x_test_fungi.npy')
y_test = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted_non_identical/y_test_fungi.npy')


# This Section is used to shuffle the data

# Aggregates elements
data_training = list(zip(x_train, y_train))
data_development = list(zip(x_dev, y_dev))
data_testing = list(zip(x_test, y_test))

# Shuffle the aggragated element on the list
random.shuffle(data_training)
random.shuffle(data_development)
random.shuffle(data_testing)

# Combine data training dan data development become one list of data train

data_train = data_training + data_development

# Split the shuffled data
x_train, y_train = zip(*data_train)
x_test, y_test = zip(*data_testing)

# Unpack the tuples
x_train = np.array(list(x_train))
y_train = np.array(list(y_train))
x_test = np.array(list(x_test))
y_test = np.array(list(y_test))


# Fix random seed for reproducibility
seed = 10
np.random.seed(seed)

# Define 5-fold cross validation test harness
kfold = StratifiedKFold(n_splits = 5, shuffle = True, random_state = seed)

# Define the convolutional layer with batch normalization
def conv1d_bn(x,
              filters,
              kernel_size,
              padding = 'same',
              strides = 1,
              activation_convolution = 'relu'):
    """Utility function to apply conv + BN.
    # Arguments
        x: input tensor.
        filters: filters in `Conv1D`.
        num_row: height of the convolution kernel.
        num_col: width of the convolution kernel.
        padding: padding mode in `Conv1D`.
        strides: strides in `Conv1D`.

    # Returns
        Output tensor after applying `Conv2D` and `BatchNormalization`.
    """

    x = layers.Conv1D(filters,
                    kernel_size,
                      strides = strides,
                      padding = padding,
                      use_bias = False,
                      kernel_regularizer = regularizers.l1(0.001),
                      bias_regularizer = regularizers.l1(0.001)
                      )(x)
    x = layers.BatchNormalization()(x)
    x = layers.Activation(activation_convolution)(x)
    return x

# Define the base model
def build_model_conv1D_lstm(filters,
              filters_LSTM,
              strides,
              padding,
              activation_convolution,
              activation_LSTM,
              optimizers,
              number_hidden_units):

  input = Input(shape = x_train.shape[1:])
  convolution_1 = conv1d_bn(x = input,
                filters = filters,
                kernel_size = 1,
                        strides = strides,
                        padding = padding,
                        activation_convolution = activation_convolution)
  convolution_2 = conv1d_bn(x = input,
                filters = filters,
                kernel_size = 3,
                        strides = strides,
                        padding = padding,
                        activation_convolution = activation_convolution)
  convolution_3 = conv1d_bn(x = input,
                        filters = filters,
                kernel_size = 5,
                        strides = strides,
                        padding = padding,
                        activation_convolution = activation_convolution)
  model = keras.layers.concatenate([convolution_1,
                                    convolution_2,
                                    convolution_3], axis=1)
  model = Conv1D(filters * 2,
               kernel_size = 3,
               strides = 1,
               padding = padding,
               activation = activation_convolution,
               use_bias = True,
               kernel_initializer = 'glorot_uniform',
               bias_initializer = 'zeros',
               kernel_regularizer = regularizers.l1(0.001),
                   bias_regularizer = regularizers.l1(0.001)
                   )(model)
  model_lstm_1 = LSTM(filters_LSTM,
                  activation = activation_LSTM,
                  recurrent_activation = 'hard_sigmoid',
                  use_bias = True,
                  kernel_initializer = 'glorot_uniform',
                  recurrent_initializer = 'orthogonal',
                  bias_initializer = 'zeros',
                  go_backwards = False)(model)

  model_lstm_2 = LSTM(filters_LSTM,
                      activation = activation_LSTM,
                      recurrent_activation = 'hard_sigmoid',
                      use_bias = True,
                      kernel_initializer = 'glorot_uniform',
                  recurrent_initializer = 'orthogonal',
                  bias_initializer = 'zeros',
                  go_backwards = True)(model)

  model_final = keras.layers.concatenate([model_lstm_1, model_lstm_2])

  output = Dense(number_hidden_units, activation = 'relu')(model_final)
  dropout = Dropout(0.5)(output)
  output = Dense(1, activation = 'sigmoid')(dropout)

  model = Model(inputs = input, outputs = output, name = 'model_1')

  model.compile(loss = 'binary_crossentropy',
                optimizer = optimizers,
                metrics = ['accuracy'])

  print(model.summary())

  return model

# Define the function to calculate sensitivity and specificity
def sensitivity_specificity(predictions, y_test, mode='binary'):
    if mode == 'binary':
        # Determine positive class predictions
        index = predictions > 0.5
        predictions = np.zeros(predictions.shape)
        predictions[index] = 1
        # No need to modify y_test since it consists of zeros and ones already
    else:
        y_test = y_test
        predictions = np.argmax(predictions, axis=-1)

    # In the binary classification case as we create, we can extract tn, fp, fn, tp as follows
    tn, fp, fn, tp = confusion_matrix(y_test, predictions, labels = [0, 1]).ravel()

    # Sensitivity = TP / (TP + FN)
    sensitivity = tp / (tp + fn)

    # Specificity = TN / (TN + FP)
    specificity = tn / (tn + fp)

    # Precision = TP / (TP + FP)
    precision = tp / (tp + fp)

    # Return sensitivity, specificity, precision
    return(sensitivity, specificity, precision)


# Define function to evaluate and predict
def evaluate_predict_fc_model():
  loss, acc = model.evaluate(x_test, y_test, verbose = 0)
  prediction = model.predict(x_test)
  sensitivity, specificity, precision = sensitivity_specificity(prediction, y_test, mode='binary')
  return acc, sensitivity, specificity, precision


# Run the model with CV = 5
# Make datalists to create dataframe of the results
result_list_train = []
result_list_val = []
result_list_test = []
col_names = ['CV1',
             'CV2',
             'CV3',
             'CV4',
             'CV5']

col_names_pred = ['acc',
                 'sensitivity',
                 'specifity']


for train, test in kfold.split(x_train, y_train):
    model = build_model_conv1D_lstm(filters = 8,
                                    filters_LSTM = 8,
                                    strides = 1,
                                    padding = "valid",
                                    activation_convolution = None,
                                    activation_LSTM = 'tanh',
                                    optimizers = 'Adadelta',
                                    number_hidden_units = 8)
    # Fit the model
    history_fit = model.fit(x_train[train], y_train[train], epochs = 30, batch_size = 16, verbose = 1, shuffle = 1, validation_data = (x_train[test], y_train[test]))

    # Get the results from the fitting
    acc_train = history_fit.history['acc']
    acc_val = history_fit.history['val_acc']

    result_line_train = np.array(acc_train)
    result_line_val = np.array(acc_val)

    # Create a list from the history of the fit
    result_list_train.append(result_line_train[:])
    result_list_val.append(result_line_val[:])

    # Evaluate the model
    acc, sensitivity, specifity, precision = evaluate_predict_fc_model()
    result_line_test = np.array((acc,
                                sensitivity,
                                specifity))

    result_list_test.append(result_line_test[:])

result_array_train = np.asarray(result_list_train)
result_array_val = np.asarray(result_list_val)
result_array_test = np.asarray(result_list_test)


# create dataframe of the result_array for train, val, and test dat
df_results_train = pd.DataFrame(np.transpose(result_array_train), columns = col_names)
df_results_val = pd.DataFrame(np.transpose(result_array_val), columns = col_names)
df_results_test = pd.DataFrame(result_array_test, columns = col_names_pred)

#  Save the dataframe into .csv
df_results_train.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-manual-train/results/fungi/df_results_train_cnn_lstm_best_fungi_filter8_filterlstm8_batch16.csv')
df_results_val.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-manual-train/results/fungi/df_results_val_cnn_lstm_best_fungi_filter8_filterlstm8_batch16.csv')
df_results_test.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-manual-train/results/fungi/df_results_test_cnn_lstm_best_fungi_filter8_filterlstm8_batch16.csv')


