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

# Import keras item
import keras

from keras.layers import *
from keras.optimizers import *
from keras.applications import *
from keras.models import *
from keras.models import Model
from keras.layers import Input, Dense
from keras.wrappers.scikit_learn import KerasClassifier

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

	model = Model(inputs = input, outputs = output)

	model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers,
                  metrics = ['accuracy'])

	print(model.summary())
	return model


# Pass the model design to KerasClassifier() wrapper

model = KerasClassifier(build_fn =  build_model_conv1D_lstm, verbose = 1)

# Define the parameters that will be tuned randomly
keras_param_options = {'filters' : [4, 8, 16],
					   'filters_LSTM' : [4, 8, 16],
					   'strides' : [1],
					   'padding' : ['valid'],
					   'activation_convolution' : [None],
					   'activation_LSTM' : ['tanh'],
					   'optimizers' : ['Adam', 'Adadelta'],
					   'number_hidden_units' : [4, 8],
					   'epochs' : [30],
					   'batch_size' : [4, 8, 16]}

# Using RandomizedSearchCV to find the best model randomly
random_search = RandomizedSearchCV(model,
                                   param_distributions = keras_param_options,
                                   n_iter = 50,
                                   cv = 5,
                                   verbose = 10)

# Fit to the training data
random_search.fit(x_train, y_train)
df_result_hyper_tuned = pd.DataFrame.from_dict(random_search.cv_results_)
df_result_hyper_tuned.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-scan-multiclass/fungi/results/all_scan_results_cnn_lstm_scan_fungi_secreted.csv')

# Save all of the params to be used to predict on the test data
df_result_hyper_tuned['mean_test_score']= pd.to_numeric(df_result_hyper_tuned['mean_test_score'])
param_best_model_dict = dict(df_result_hyper_tuned.nlargest(30, 'mean_test_score')['params'])
params = list(param_best_model_dict.values())
print(params)

# Get info ahead about the best model obtained
print('Best score obtained: {0}'.format(random_search.best_score_))
print('Parameters:')
for param, value in random_search.best_params_.items():
    print('\t{}: {}'.format(param, value))


# Predict the prediction of the best model
print('Predict using test data using random_search:')
y_pred_random_search = random_search.predict(x_test)
acc_pred_random_search = accuracy_score(y_test, y_pred_random_search)
print('acc y_pred_random_search:', acc_pred_random_search)

# Predict the results of hyperparamaters tuning for all parameters

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


# Make prediction of test data

result_list = []
columns_names = ['Parameters',
                'Accuracy',
                'Sensitivity',
                'Specifity']


# For loop to train model and get prediction for each combination of dataset
for i in range(len(params)):
    list_par = list(params[i].values())
    # Define the models
    model = build_model_conv1D_lstm(filters = list_par[5],
									filters_LSTM = list_par[4],
									strides = list_par[0],
									padding = list_par[1],
									activation_convolution = list_par[8],
									activation_LSTM = list_par[9],
									optimizers = list_par[2],
									number_hidden_units = list_par[3])

    # Train the model one by one based on the parameters combination
    train_fc_model(batch_sizes = list_par[7], num_epochs = list_par[6])
    acc, sensitivity, specifity, precision = evaluate_predict_fc_model()
    result_line = np.array((params[i],
                            acc,
                            sensitivity,
                            specifity))
    result_list.append(result_line[:])
    result_array = np.asarray(result_list)

    df_results = pd.DataFrame(result_array,
                              columns = columns_names)

df_results.to_csv('/hpc-home/kristian/effector-non-effector/scripts-cnn-lstm-separate-group/scripts-scan-multiclass/fungi/results/df_pred_results_cnn_lstm_scan_fungi_secreted.csv')


