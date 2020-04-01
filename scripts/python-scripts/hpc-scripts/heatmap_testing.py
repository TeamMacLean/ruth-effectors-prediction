import tensorflow as tf
from tensorflow.keras.models import load_model
tf.disable_v2_behavior()
import cv2
import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

import tensorflow.keras.backend as K

# Load the model
model = load_model()



def get_heatmap_matrix(dataset, nth_data, layer, verbose = False):

    #Get the data
    data = dataset[nth_data:nth_data+1, :, :]
    if verbose: print(data.shape)

    # Get the prediction from the model
    preds = model.predict(data)

    # Get the position that maximally activated
    position = np.argmax(preds[0])
    if verbose: print(position)

    # Get the output of the layer and the data
    data_output = model.output[:, 0]

    # Get layer output
    get_layer = model.get_layer(layer)
    layer_output = get_layer.output

    # Calculate the gradients
    grads = K.gradients(data_output, layer_output)[0]
    if verbose: print(grads.shape)

    pooled_grads = K.mean(grads, axis = (0, 1))

    iterate = K.function([model.input], [pooled_grads, layer_output[0]])

    if verbose: print("Starting K iteration...")
    t = time.time()
    pooled_grads_value, layer_output_value = iterate([data])
    if verbose: print("Finished K iteration " + "(" + str(time.time() - t) + "seconds)")
    if verbose: print(layer_output_value.shape)

    if verbose: print("Starting loop...")
    t = time.time()
    # Looping
    for i in range(32):
        layer_output_value[:, i] *= pooled_grads_value[i]

    if verbose: print("Loop finished " + " (" + str(time.time() - t) + "seconds)" )
    # Get the heatmap matrix
    heatmap = np.mean(layer_output_value, axis = -1)
    heatmap = np.maximum(heatmap, 0)
    heatmap /= np.max(heatmap)

    # Expand the dimensionality of heatmap so that it can be plot
    heatmap = np.expand_dims(heatmap, axis=0)

    return heatmap


# Create list of data name that will be loaded:

list_data = ["x_train", "x_test", "x_val"]

for nth_list in range(len(list_data)):

    # Define pattern of the models
    dataset = list_data[nth_list]

    data_loading_pattern = dataset + ".npy"
    data_name_pattern = "data_" + dataset

    data_loading_path = glob.glob(data_loading_pattern)
    data_name = glob.glob(data_name_pattern)

    # Load the model
    print("Loading data:", dataset)
    data_name = np.load(data_loading_path)


# Getting all of the heatmap for all of the data together in one matrices

# Create list of data name that will be loaded:
result = []

def getting_sum_heatmap(data, layer):

    # For loop to get the heatmap for each matrix
    for nth_sample in data.shape[0]:
        print("Getting the heatmap for the data:", nth_sample)
        heatmap = get_heatmap_matrix(data, nth_sample, layer, verbose = False)

        # Making sure that the shape is correct which is (1, 4034)
        print("Heatmap size:", heatmap.shape)

        # Put all of the reults together
        result.append(heatmap)

        # Change the list to numpy array
        all_matrices = np.array(result)

        # Get sum of all matrices
        sum_all_matrices = np.sum(all_matrices, axis = 0)

    # Save the sum of all matrices
    saving_pattern = "sum_all_matrices_" + data + ".npy"
    saving_path = glob.glob(saving_pattern)
    np.save(sum_all_matrices, saving_path)







