# import all packages and library --------------------------------------------------

# import package to scan hyperparameter
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint as sp_randint
from sklearn.metrics import confusion_matrix

# import package to reprocess the data
import numpy as np
import pandas as pd
import random

# import properties from keras
from keras import models
from keras.layers import Dense, Dropout
from keras import regularizers

# import keras items
from keras.optimizers import Adam, Adadelta, SGD
from keras.activations import relu, elu, sigmoid
from keras.losses import binary_crossentropy
from keras.layers import Dense, Dropout
from keras.layers.normalization import BatchNormalization
from keras.wrappers.scikit_learn import KerasClassifier
