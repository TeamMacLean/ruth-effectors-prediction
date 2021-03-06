#!/usr/bin/env python3

# Load of the library that needed

# Import package to scan hyperparameter
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint as sp_randint
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from keras.models import load_model
from sklearn.externals import joblib

# Import package to reprocess the data
import numpy as np
import random
import pandas as pd
import ast

# Import keras item
import keras

from keras.layers import *
from keras.optimizers import *
from keras.applications import *
from keras.models import *
from keras.models import Model
from keras.layers import Input, Dense
from keras import regularizers
from keras.wrappers.scikit_learn import KerasClassifier
from keras.layers.embeddings import Embedding



# Get the reprocessed data from .npy file
x_train = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-bacteria/secreted_non_identical/x_train_int_bacteria.npy')
y_train = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-bacteria/secreted_non_identical/y_train_bacteria.npy')
x_dev = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-bacteria/secreted_non_identical/x_val_int_bacteria.npy')
y_dev = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-bacteria/secreted_non_identical/y_val_bacteria.npy')
x_test = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-bacteria/secreted_non_identical/x_test_int_bacteria.npy')
y_test = np.load('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/data-bacteria/secreted_non_identical/y_test_bacteria.npy')


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

def simple_LSTM(inputdim = 23,
               outputdim = 32,
               lstm_hidden_units = 16,
               opt_dropout = 0.5,
               opt_dropout_recurrent = 0.5,
               opt_go_backwards = 'FALSE',
               reg_rate = 0.01,
               optimizers = 'Adam'
              ):
  # create the model
    model = Sequential()
    model.add(Embedding(input_dim = inputdim,
                        output_dim = outputdim,
                        input_length = 2574))
    model.add(Bidirectional(LSTM(units = lstm_hidden_units,
                                activation = 'tanh',
                                recurrent_activation = 'hard_sigmoid',
                                kernel_regularizer = regularizers.l2(reg_rate),
                                recurrent_regularizer = regularizers.l2(reg_rate),
                                dropout = opt_dropout,
                                recurrent_dropout = opt_dropout_recurrent,
                                go_backwards = opt_go_backwards)))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers,
                  metrics = ['accuracy'])
    print(model.summary())

    return model

# Pass the model design to KerasClassifier() wrapper
model = KerasClassifier(build_fn =  simple_LSTM, verbose = 1)

# Define the parameters that will be tuned randomly
keras_param_options = {'outputdim' : [16, 32, 64],
                       'lstm_hidden_units' : [16, 32],
                       'opt_dropout' : [0, 0.25],
                       'opt_dropout_recurrent' : [0],
                       'opt_go_backwards' : ['TRUE', 'FALSE'],
                       'reg_rate' : [0.01, 0.001],
                       'optimizers' : ['Adam', 'sgd'],
                       'batch_size' : [4, 8],
                       'epochs' : [30]}

# using RandomizedSearchCV to find the best model randomly
random_search = RandomizedSearchCV(model,
                                   param_distributions = keras_param_options,
                                   n_iter = 50,
                                   cv = 5,
                                   verbose = 10)


 # Fit to the training data
random_search.fit(x_train, y_train)
import numpy as np
df_result_hyper_tuned = pd.DataFrame.from_dict(random_search.cv_results_)
df_result_hyper_tuned.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-scan-multiclass/fungi/results/all_scan_results_lstm_emb_scan_bacteria.csv.csv')


# Save all of the params to be used to predict
df_result_hyper_tuned['mean_test_score']= pd.to_numeric(df_result_hyper_tuned['mean_test_score'])
param_best_model_dict = df_result_hyper_tuned.nlargest(10, 'mean_test_score')['params']
print(param_best_model_dict)


# Define function to fit the model
def train_fc_model(batch_sizes = None, num_epochs = None):
      model.fit(x = x_train,
                y = y_train,
                batch_size = batch_sizes,
                epochs = num_epochs,
                verbose = 1,
                shuffle = 1)

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

# Make prediction
result_list = []
columns_names = ['Parameters',
                'Accuracy',
                'Sensitivity',
                'Specifity']

    model = simple_LSTM(inputdim = 23,
                       outputdim = list_par[1],
                       lstm_hidden_units = list_par[6],
                       opt_dropout = list_par[5],
                       opt_dropout_recurrent = list_par[4],
                       opt_go_backwards = list_par[3],
                       reg_rate = list_par[0],
                       optimizers = list_par[2]
                       )
    train_fc_model(batch_sizes = list_par[8], num_epochs = list_par[7])
    acc, sensitivity, specifity, precision = evaluate_predict_fc_model()
    result_line = np.array((param_best_model_dict[i],
                            acc,
                            sensitivity,
                            specifity))
    result_list.append(result_line[:])
    result_array = np.asarray(result_list)

    df_results = pd.DataFrame(result_array,
                              columns = columns_names)

df_results.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-scan-multiclass/bacteria/results/df_pred_results_lstm_emb_scan_bacteria.csv')


