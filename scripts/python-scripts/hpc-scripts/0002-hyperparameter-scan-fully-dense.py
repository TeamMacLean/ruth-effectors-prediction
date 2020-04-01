# Import all packages and library

# Import package to scan hyperparameter
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import confusion_matrix

# Import package to reprocess the data
import numpy as np
import pandas as pd
import random

# Import properties from keras
from keras import models
from keras.layers import Dense, Dropout, Activation
from keras import regularizers

# Import keras items
from keras.optimizers import Adam, Adadelta, SGD
from keras.activations import relu, sigmoid
from keras.layers.advanced_activations import PReLU
from keras.losses import binary_crossentropy
from keras.layers.normalization import BatchNormalization
from keras.wrappers.scikit_learn import KerasClassifier

# Get all of the data and reprocess them

# Get the reprocessed data from .npy file
x_train = np.load('../r-scripts/getting-data-current/data-sets/x_train.npy')
y_train = np.load('../r-scripts/getting-data-current/data-sets/y_train.npy')

x_dev = np.load('../r-scripts/getting-data-current/data-sets/x_val.npy')
y_dev = np.load('../r-scripts/getting-data-current/data-sets/y_val.npy')

x_test = np.load('../r-scripts/getting-data-current/data-sets/x_test.npy')
y_test = np.load('../r-scripts/getting-data-current/data-sets/y_test.npy')

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

# Reshape the datasets
x_train = x_train.reshape(615, 4034 * 20)
x_test = x_test.reshape(150, 4034 * 20)


# Define the model and function

# Define the model base
def build_fc_model(# Hyperparameters as parts of model designs (building blocks)
	               input_num_hidden_units = 2,
                   num_hidden_layers = [0],
                   activation_function = 'relu',

                   # Hyperparameters as part of optimization and regularization of the models
                   l2_rate = 0.001,
                   input_dropout_rates = 0.5,
                   dropout_rates = 0.5,
                   optim_methods = 'Adam',
                   batch_norm = "yes"
                   ):

	# Add the input layer
    model = models.Sequential()
    model.add(Dense(input_num_hidden_units,
              activation = 'relu',
              kernel_regularizer = regularizers.l2(l2_rate),
              input_dim = x_train.shape[1]))
    model.add(Dropout(input_dropout_rates))

    # Add the hidden layers
    for num in range(len(num_hidden_layers)):
        if num_hidden_layers[num] == 0:
            continue
        else:
            model.add(Dense(num_hidden_layers[num]))

            # Add batch normaization before adding the activation layers
            if batch_norm == "yes":
            	model.add(BatchNormalization())
            else:
            	continue

            model.add(Activation(activation_function))
            model.add(Dropout(dropout_rates))

    # Add the output layer
    model.add(Dense(1,
              activation = 'sigmoid'))

    # Compile the model defined
    model.compile(optimizer = optim_methods,
                  loss = 'binary_crossentropy',
                  metrics = ['acc'])

    # Print the summary of the model
    print(model.summary())

    return model

# Pass the model design to KerasClassifier() wrapper -------------------------------

model = KerasClassifier(build_fn = build_fc_model,
	                         verbose = 1)

# Define the parameters that will be tuned randomly
keras_param_options = {
                       # Hyperparameters as parts of model designs (building blocks)
                       'input_num_hidden_units': [1, 2, 3, 4, 8, 16],
                       'num_hidden_layers': [[0], [1], [2], [3], [1, 2], [1, 3], [2, 2], [2, 3]],
                       'activation_function': ['relu', 'PReLU'],
                       # Hyperparameters as part of optimization and regularization of the models
                       'optim_methods' : ['Adam', 'Adadelta', 'SGD'],
                       'l2_rate':[0.001, 0.01],
                       'input_dropout_rates': [0.25, 0.5],
                       'dropout_rates': [0.25, 0.5],
                       'batch_norm' : ['yes', 'no'],
                       # Fitting parameters
                       'batch_size': [8, 16, 32],
                       'epochs': [20, 30],
                       'shuffle': [True]
                      }

# Using RandomizedSearchCV to find the best model randomly
random_search = RandomizedSearchCV(model,
                                   param_distributions = keras_param_options,
                                   return_train_score=True,
                                   n_iter = 100,
                                   cv = 5,
                                   verbose = 10)

# Fit to the training data
random_search.fit(x_train, y_train)
df_result_hyper_tuned = pd.DataFrame.from_dict(random_search.cv_results_)
df_result_hyper_tuned.to_csv('result_hyper_tuned.csv')

# save all of the params to be used to predict
df_result_hyper_tuned['mean_test_score']= pd.to_numeric(df_result_hyper_tuned['mean_test_score'])
param_best_model_dict = dict(df_result_hyper_tuned.nlargest(30, 'mean_test_score')['params'])
params = list(param_best_model_dict.values())
print(params)

# Predict the results of hyperparamaters tuning

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

# make prediction
result_list = []
columns_names = ['Parameters',
                'Accuracy',
                'Sensitivity',
                'Specifity']

for i in range(len(params)):
    list_par = list(params[i].values())
    model = build_fc_model(input_num_hidden_units = list_par[4],
                           num_hidden_layers = list_par[2],
                           activation_function = list_par[10],
                           l2_rate = list_par[3],
                           input_dropout_rates = list_par[5],
                           dropout_rates = list_par[7],
                           optim_methods = list_par[1],
                           batch_norm = list_par[9])
    train_fc_model(batch_sizes = list_par[8], num_epochs = list_par[6])
    acc, sensitivity, specifity, precision = evaluate_predict_fc_model()
    result_line = np.array((params[i],
                            acc,
                            sensitivity,
                            specifity))
    result_list.append(result_line[:])
    result_array = np.asarray(result_list)

    df_results = pd.DataFrame(result_array,
                         columns = columns_names)

df_results.to_csv('df_result_prediction.csv')
