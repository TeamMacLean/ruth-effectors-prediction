Smaller network on deep learning models
=======================================

``` r
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(taxize)
library(caret)
```

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
reticulate::use_condaenv(condaenv = "tensorflow2", conda = "/anaconda3/bin/conda")
```

``` python
sys.executable
```

    ## '/anaconda3/envs/tensorflow2/bin/python'

``` python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from tensorflow.compat.v1.keras.models import load_model
```

CNN-GRU models
--------------

``` python
# Get the pretarined model and show the summary of the model
model_cnn_gru = load_model("../../../../data/getting-data-new/multi-class-data/data-sets/models/bacteria/cnn-gru/sequential_1.30-0.27.hdf5")
model_cnn_gru.summary()
```

    ## Model: "sequential_1"
    ## _________________________________________________________________
    ## Layer (type)                 Output Shape              Param #   
    ## =================================================================
    ## conv1d_1 (Conv1D)            (None, 4033, 16)          656       
    ## _________________________________________________________________
    ## max_pooling1d_1 (MaxPooling1 (None, 1344, 16)          0         
    ## _________________________________________________________________
    ## bidirectional_1 (Bidirection (None, 16)                1200      
    ## _________________________________________________________________
    ## dense_1 (Dense)              (None, 1)                 17        
    ## =================================================================
    ## Total params: 1,873
    ## Trainable params: 1,873
    ## Non-trainable params: 0
    ## _________________________________________________________________

Now, we can get the results of the testing data of each class: Bacteria,
Oomycete, and Fungi

``` r
cnn_gru_bacteria <- data.table::fread("../../../../results/cnn-gru-separate-class/bacteria/df_results_test_cnn_gru_saved_model1.csv", drop = "V1") %>% select("acc")

cnn_gru_fungi <- data.table::fread("../../../../results/cnn-gru-separate-class/fungi/df_results_test_cnn_gru_saved_model1.csv", drop = "V1") %>% select("acc")

cnn_gru_oomycete <- data.table::fread("../../../../results/cnn-gru-separate-class/oomycete/df_results_test_cnn_gru_saved_model1.csv", drop = "V1") %>% select("acc")
```

``` r
cnn_gru_oomycete
```

    ##          acc
    ## 1: 0.6842105

CNN-LSTM models
---------------

``` python
# Get the pretarined model and show the summary of the model
model_cnn_lstm = load_model("../../../../data/getting-data-new/multi-class-data/data-sets/models/bacteria/cnn-lstm/model_1.30-0.41.hdf5")
model_cnn_lstm.summary()
```

    ## Model: "model_1"
    ## __________________________________________________________________________________________________
    ## Layer (type)                    Output Shape         Param #     Connected to                     
    ## ==================================================================================================
    ## input_1 (InputLayer)            [(None, 4034, 20)]   0                                            
    ## __________________________________________________________________________________________________
    ## conv1d_1 (Conv1D)               (None, 4034, 4)      80          input_1[0][0]                    
    ## __________________________________________________________________________________________________
    ## conv1d_2 (Conv1D)               (None, 4032, 4)      240         input_1[0][0]                    
    ## __________________________________________________________________________________________________
    ## conv1d_3 (Conv1D)               (None, 4030, 4)      400         input_1[0][0]                    
    ## __________________________________________________________________________________________________
    ## batch_normalization_1 (BatchNor (None, 4034, 4)      16          conv1d_1[0][0]                   
    ## __________________________________________________________________________________________________
    ## batch_normalization_2 (BatchNor (None, 4032, 4)      16          conv1d_2[0][0]                   
    ## __________________________________________________________________________________________________
    ## batch_normalization_3 (BatchNor (None, 4030, 4)      16          conv1d_3[0][0]                   
    ## __________________________________________________________________________________________________
    ## activation_1 (Activation)       (None, 4034, 4)      0           batch_normalization_1[0][0]      
    ## __________________________________________________________________________________________________
    ## activation_2 (Activation)       (None, 4032, 4)      0           batch_normalization_2[0][0]      
    ## __________________________________________________________________________________________________
    ## activation_3 (Activation)       (None, 4030, 4)      0           batch_normalization_3[0][0]      
    ## __________________________________________________________________________________________________
    ## concatenate_1 (Concatenate)     (None, 12096, 4)     0           activation_1[0][0]               
    ##                                                                  activation_2[0][0]               
    ##                                                                  activation_3[0][0]               
    ## __________________________________________________________________________________________________
    ## conv1d_4 (Conv1D)               (None, 12094, 8)     104         concatenate_1[0][0]              
    ## __________________________________________________________________________________________________
    ## lstm_1 (LSTM)                   (None, 8)            544         conv1d_4[0][0]                   
    ## __________________________________________________________________________________________________
    ## lstm_2 (LSTM)                   (None, 8)            544         conv1d_4[0][0]                   
    ## __________________________________________________________________________________________________
    ## concatenate_2 (Concatenate)     (None, 16)           0           lstm_1[0][0]                     
    ##                                                                  lstm_2[0][0]                     
    ## __________________________________________________________________________________________________
    ## dense_1 (Dense)                 (None, 4)            68          concatenate_2[0][0]              
    ## __________________________________________________________________________________________________
    ## dropout_1 (Dropout)             (None, 4)            0           dense_1[0][0]                    
    ## __________________________________________________________________________________________________
    ## dense_2 (Dense)                 (None, 1)            5           dropout_1[0][0]                  
    ## ==================================================================================================
    ## Total params: 2,033
    ## Trainable params: 2,009
    ## Non-trainable params: 24
    ## __________________________________________________________________________________________________

``` r
cnn_lstm_bacteria <- data.table::fread("../../../../results/cnn-lstm-separate-class/bacteria/df_results_test_saved_model.csv", drop = "V1") %>% select("acc")

cnn_lstm_fungi <- data.table::fread("../../../../results/cnn-lstm-separate-class/fungi/df_results_test_saved_model.csv", drop = "V1") %>% select("acc")

cnn_lstm_oomycete <- data.table::fread("../../../../results/cnn-lstm-separate-class/oomycete/df_results_test_saved_model.csv", drop = "V1") %>% select("acc")
```

### Summary of results

``` r
data.frame(
  name = c("Bacteria", "Oomycete", "Fungi"),
  acc_cnn_gru = c(cnn_gru_bacteria[["acc"]], cnn_gru_oomycete[["acc"]], cnn_gru_fungi[["acc"]]),
  acc_cnn_lstm = c(cnn_lstm_bacteria[["acc"]], cnn_lstm_fungi[["acc"]], cnn_lstm_oomycete[["acc"]])
)
```

    ##       name acc_cnn_gru acc_cnn_lstm
    ## 1 Bacteria   0.8552632    0.8026316
    ## 2 Oomycete   0.6842105    0.5227273
    ## 3    Fungi   0.6363636    0.7368421
