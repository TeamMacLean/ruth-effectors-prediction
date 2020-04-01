#!/usr/bin/env python3

# Load of the library that needed

# Import package to scan hyperparameter
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint as sp_randint
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score

# Import package to reprocess the data
import numpy as np
import random
import pandas as pd
import glob
import os

# Import keras item
import keras

from keras.layers import *
from keras.optimizers import *
from keras.applications import *
from keras.models import *
from keras.models import Model
from keras.layers import Input, Dense, Dropout
from keras import regularizers

# Import keras properties to save models
from keras.callbacks import History
from keras.callbacks import ModelCheckpoint

# Get all of the data and reprocess them

# Get the reprocessed data from .npy file
x_train = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted/x_train_fungi.npy')
y_train = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/data-fungi/y_train_fungi.npy')

x_dev = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted/x_val_fungi.npy')
y_dev = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/data-fungi/y_val_fungi.npy')

x_test = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-fungi/secreted/x_val_fungi.npy')
y_test = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/data-fungi/y_test.npy')

# This Section is used to shuffle the data

# Aggregates elements
data_training = list(zip(x_train, y_train))
data_development = list(zip(x_dev, y_dev))
data_testing = list(zip(x_test, y_test))

# Shuffle the aggragated element on the list
random.shuffle(data_training)
random.shuffle(data_development)
random.shuffle(data_testing)


# Split the shuffled data
x_train, y_train = zip(*data_training)
x_val, y_val = zip(*data_development)
x_test, y_test = zip(*data_testing)

# Unpack the tuples
x_train = np.array(list(x_train))
y_train = np.array(list(y_train))
x_val = np.array(list(x_val))
y_val = np.array(list(y_val))
x_test = np.array(list(x_test))
y_test = np.array(list(y_test))


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
                      # activity_regularizer = regularizers.l1(0.01)
                      )(x)
    x = layers.BatchNormalization()(x)
    x = layers.Activation(activation_convolution)(x)
    return x

# Define the base model

def build_model_conv1D_lstm(filters = 48,
                            filters_LSTM = 48,
                            strides = 1,
                            padding = "valid",
                            activation_convolution = None,
                            activation_LSTM = 'tanh',
                            number_hidden_units = 64,
                            x_train = x_train):

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

    model = Model(inputs = input, outputs = output, name = 'cnn_lstm_fungi_secreted_data')

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
def evaluate_predict_fc_model(model):
  loss, acc = model.evaluate(x_test, y_test, verbose = 0)
  prediction = model.predict(x_test)
  sensitivity, specificity, precision = sensitivity_specificity(prediction, y_test, mode='binary')
  return acc, sensitivity, specificity, precision


# Define the column header for the dataframes

col_names = ['acc',
             'loss']

col_names_pred = ['acc',
                 'sensitivity',
                 'specifity']


model = build_model_conv1D_lstm(filters = 4,
                                filters_LSTM = 8,
                                strides = 1,
                                padding = 'valid',
                                activation_convolution = None,
                                activation_LSTM = 'tanh',
                                number_hidden_units = 4)

def compile_and_train(model, num_batchs, num_epochs):

    model.compile(loss = 'binary_crossentropy', optimizer = 'Adadelta', metrics=['acc'])
    filepath = '/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/weights/' + model.name + '.{epoch:02d}-{loss:.2f}.hdf5'
    checkpoint = ModelCheckpoint(filepath, monitor='loss', verbose=0, save_weights_only=False, save_best_only = False, mode='auto', period=1)
    history = model.fit(x = x_train,
                        y = y_train,
                        batch_size = num_batchs,
                        epochs = num_epochs,
                        verbose = 1,
                        callbacks = [checkpoint],
                        validation_data = (x_val, y_val))
    weight_files = glob.glob(os.path.join(os.getcwd(), '/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/weights/secreted/*'))
    weight_file = max(weight_files, key = os.path.getctime) # most recent file
    return history, weight_file


# Compile, fit, and evaluate the model
history_fit, model_1_weight_file = compile_and_train(model, num_batchs = 8, num_epochs = 30)

# Load the weight
model.load_weights(model_1_weight_file)
acc, sensitivity, specifity, precision = evaluate_predict_fc_model(model)

# Saving the results and history into dataframe
result_train = []
result_val = []

acc_train = np.array(history_fit.history['acc'])
acc_val = np.array(history_fit.history['val_acc'])
loss_train = np.array(history_fit.history['loss'])
loss_val = np.array(history_fit.history['val_loss'])


result_train = np.concatenate((acc_train.reshape(30, 1), loss_train.reshape(30, 1)), axis = 1)
result_val = np.concatenate((acc_val.reshape(30, 1), loss_val.reshape(30, 1)), axis = 1)
result_test = np.array((acc, sensitivity, specifity)).reshape(1, 3)

# create dataframe of the result_array for train, val, and test dat
df_results_train = pd.DataFrame(result_train, columns = col_names)
df_results_val = pd.DataFrame(result_val, columns = col_names)
df_results_test = pd.DataFrame(result_test, columns = col_names_pred)

#  Save the dataframe into .csv
df_results_train.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/results/secreted/df_results_train_saved_cnn_lstm_fungi_secreted_data.csv')
df_results_val.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/results/secreted/df_results_val_saved_cnn_lstm_fungi_secreted_data.csv')
df_results_test.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/fungi/results/secreted/df_results_test_saved_cnn_lstm_fungi_secreted_data.csv')
