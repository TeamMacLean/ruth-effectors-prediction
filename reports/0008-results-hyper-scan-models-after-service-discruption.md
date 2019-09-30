Current Results of Hyperparameters Scanning after the jobs got killed on 20 – 22 September 2019
===============================================================================================

Introduction
------------

According to the last results after there was an HPC service
discruption, there are severak items need to be done according to Dan’s
suggestions:

1.  Do prediction to the test data for the hyperparamaters space that
    gave the best accuracy, for the model:

    -   CNN-LSTM
    -   CNN-GRU

2.  Do grid search on the model LSTM

3.  Do again RandomizedSearchCV on the model GRU

Execution and Results
---------------------

### Evaluating on the test data for the model with the combination of hyperparamaters that gives the best accuracy training dan validation data

#### CNN - LSTM with the best parameters

    Predict using test data using random_search:

      8/150 [>.............................] - ETA: 23s
     16/150 [==>...........................] - ETA: 15s
     24/150 [===>..........................] - ETA: 13s
     32/150 [=====>........................] - ETA: 11s
     40/150 [=======>......................] - ETA: 10s
     48/150 [========>.....................] - ETA: 9s 
     56/150 [==========>...................] - ETA: 8s
     64/150 [===========>..................] - ETA: 7s
     72/150 [=============>................] - ETA: 6s
     80/150 [===============>..............] - ETA: 5s
     88/150 [================>.............] - ETA: 5s
     96/150 [==================>...........] - ETA: 4s
    104/150 [===================>..........] - ETA: 3s
    112/150 [=====================>........] - ETA: 3s
    120/150 [=======================>......] - ETA: 2s
    128/150 [========================>.....] - ETA: 1s
    136/150 [==========================>...] - ETA: 1s
    144/150 [===========================>..] - ETA: 0s
    150/150 [==============================] - 12s 81ms/step
    acc y_pred_random_search: 0.7133333333333334

#### CNN - GRU with the best paramaters

    Parameters:
        optimizers: Adam
        opt_recurrent_regs: <keras.regularizers.L1L2 object at 0x7f7d0257ed68>
        opt_kernel_regs: <keras.regularizers.L1L2 object at 0x7f7d0257ed68>
        opt_go_backwards: TRUE
        opt_dropout_recurrent: 0
        opt_dropout: 0
        maxpool_size: 3
        kernel_size: 2
        gru_hidden_units: 64
        filter_conv: 48
        epochs: 30
        batch_size: 8
        activation_conv: relu

      8/150 [>.............................] - ETA: 38s
     16/150 [==>...........................] - ETA: 24s
     24/150 [===>..........................] - ETA: 19s
     32/150 [=====>........................] - ETA: 16s
     40/150 [=======>......................] - ETA: 13s
     48/150 [========>.....................] - ETA: 12s
     56/150 [==========>...................] - ETA: 10s
     64/150 [===========>..................] - ETA: 9s 
     72/150 [=============>................] - ETA: 8s
     80/150 [===============>..............] - ETA: 7s
     88/150 [================>.............] - ETA: 6s
     96/150 [==================>...........] - ETA: 5s
    104/150 [===================>..........] - ETA: 4s
    112/150 [=====================>........] - ETA: 3s
    120/150 [=======================>......] - ETA: 3s
    128/150 [========================>.....] - ETA: 2s
    136/150 [==========================>...] - ETA: 1s
    144/150 [===========================>..] - ETA: 0s
    150/150 [==============================] - 15s 102ms/step
    acc: 0.7133333333333334

### CNN - LSTM parallel

According to those results above, using the models with the combination
of hyperparameters with the best results did not really give the good or
massive improvements on the accuracy (comparing with the best results
from fully connected dense layers). Therefore, for these models, a grid
search will be done, and the scripts will be divided into 4 to minimize
the time for the jobs to run:

Here are the hyperparamaters and parameters division for each script:

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 7%" />
<col style="width: 6%" />
<col style="width: 10%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 14%" />
<col style="width: 10%" />
<col style="width: 13%" />
<col style="width: 13%" />
<col style="width: 4%" />
</colgroup>
<thead>
<tr class="header">
<th>Script</th>
<th>Batch</th>
<th>Hidden Unit</th>
<th>Filters</th>
<th>Filters LSTM</th>
<th>Strides</th>
<th>Padding</th>
<th>Activation Convolution</th>
<th>Activation LSTM</th>
<th>Optimizers</th>
<th>Batch Normalizations</th>
<th>Epochs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>16</td>
<td>16</td>
<td>[16, 32]</td>
<td>[8, 16, 32, 48]</td>
<td>1</td>
<td>Valid</td>
<td>None</td>
<td>tanh</td>
<td>[‘Adam’, ‘Adadelta’]</td>
<td>yes</td>
<td>30</td>
</tr>
<tr class="even">
<td>2</td>
<td>16</td>
<td>32</td>
<td>[16, 32]</td>
<td>[8, 16, 32, 48]</td>
<td>1</td>
<td>Valid</td>
<td>None</td>
<td>tanh</td>
<td>[‘Adam’, ‘Adadelta’]</td>
<td>yes</td>
<td>30</td>
</tr>
<tr class="odd">
<td>3</td>
<td>32</td>
<td>16</td>
<td>[16, 32]</td>
<td>[8, 16, 32, 48]</td>
<td>1</td>
<td>Valid</td>
<td>None</td>
<td>tanh</td>
<td>[‘Adam’, ‘Adadelta’]</td>
<td>yes</td>
<td>30</td>
</tr>
<tr class="even">
<td>4</td>
<td>32</td>
<td>32</td>
<td>[16, 32]</td>
<td>[8, 16, 32, 48]</td>
<td>1</td>
<td>Valid</td>
<td>None</td>
<td>tanh</td>
<td>[‘Adam’, ‘Adadelta’]</td>
<td>yes</td>
<td>30</td>
</tr>
</tbody>
</table>

#### Script 1

``` r
src_file_cnn_lstm_1 <- "../../../results/results/0001-cnn-lstm-scan-batch16-unit16.log"

src_params_list_cnn_lstm <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw", "bn",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params_cnn_lstm <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)
```

``` r
data_cnn_lstm_1 <- src_file_cnn_lstm_1 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm, src_numeric_params_cnn_lstm) %>%
  clean_log_df_with_params()
```

    ## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
    ## This warning is displayed once per session.

    ## Warning: Expected 7 pieces. Additional pieces discarded in 118 rows [385,
    ## 465, 466, 467, 468, 469, 470, 2353, 2369, 2370, 2371, 2959, 3145, 3151,
    ## 3152, 3153, 3176, 3177, 3178, 3179, ...].

``` r
new_data_cnn_lstm_1 <- data_cnn_lstm_1 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm)
```

``` r
new_data_cnn_lstm_1 %>% 
  dplyr::select(-c(loss_sd, acc_sd)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw| bn  |  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  mean\_train\_score|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|:----|------------:|:------------------------|:-----------------|-----------:|-------------------:|
| 1     | 30     |        1| valid   | Adadelta   |                     16|             48|       32|           30| yes |           32| None                    | tanh             |     0.01180|             0.99510|
| 2     | 30     |        1| valid   | Adadelta   |                     32|             16|       48|           30| yes |           16| None                    | tanh             |     0.00894|             0.99674|
| 3     | 30     |        1| valid   | Adadelta   |                     16|              8|       16|           30| yes |           32| None                    | tanh             |     0.02870|             0.99595|

#### Script 2

``` r
src_file_cnn_lstm_2 <- "../../../results/results/0002-cnn-lstm-scan-batch16-unit32.log"
```

``` r
data_cnn_lstm_2 <- src_file_cnn_lstm_2 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm, src_numeric_params_cnn_lstm) %>%
  clean_log_df_with_params()
```

    ## Warning: Expected 7 pieces. Additional pieces discarded in 25 rows [6418,
    ## 6449, 6450, 6451, 6452, 6453, 6454, 6455, 6456, 6457, 6480, 6481, 6482,
    ## 8185, 8309, 8310, 9088, 9089, 9115, 9240, ...].

``` r
new_data_cnn_lstm_2 <- data_cnn_lstm_2 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm)
```

``` r
new_data_cnn_lstm_2 %>% 
  dplyr::select(-c(loss_sd, acc_sd)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw| bn  |  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  mean\_train\_score|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|:----|------------:|:------------------------|:-----------------|-----------:|-------------------:|
| 1     | 30     |        1| valid   | Adadelta   |                     32|              8|       16|           30| yes |           16| None                    | tanh             |     0.02146|             0.99390|
| 2     | 30     |        1| valid   | Adadelta   |                     32|             16|       32|           30| yes |           16| None                    | tanh             |     0.01046|             0.99594|
| 3     | 30     |        1| valid   | Adadelta   |                     32|             16|       16|           30| yes |           16| None                    | tanh             |     0.01075|             0.99695|

#### Script 3

``` r
src_file_cnn_lstm_3 <- "../../../results/results/0003-cnn-lstm-scan-batch32-unit16.log"
```

``` r
data_cnn_lstm_3 <- src_file_cnn_lstm_3 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm, src_numeric_params_cnn_lstm) %>%
  clean_log_df_with_params()

new_data_cnn_lstm_3 <- data_cnn_lstm_3 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm)
```

``` r
new_data_cnn_lstm_3 %>% 
  dplyr::select(-c(loss_sd, acc_sd)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw| bn  |  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  mean\_train\_score|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|:----|------------:|:------------------------|:-----------------|-----------:|-------------------:|
| 1     | 30     |        1| valid   | Adadelta   |                     16|             16|       32|           30| yes |           32| None                    | tanh             |   0.0140800|           0.9959200|
| 2     | 30     |        1| valid   | Adadelta   |                     16|              8|       16|           30| yes |           32| None                    | tanh             |   0.0482000|           0.9890400|
| 3     | 30     |        1| valid   | Adam       |                     16|              8|       16|           30| yes |           32| None                    | tanh             |   0.0857800|           0.9821000|
| 4     | 30     |        1| valid   | Adam       |                     16|             32|       16|           30| yes |           32| None                    | tanh             |   0.0136600|           0.9959200|
| 5     | 30     |        1| valid   | Adam       |                     16|             16|       16|           30| yes |           32| None                    | tanh             |   0.0275600|           0.9951200|
| 6     | 30     |        1| valid   | Adadelta   |                     16|             16|       16|           30| yes |           32| None                    | tanh             |   0.0168800|           0.9963400|
| 7     | 30     |        1| valid   | Adam       |                     16|              8|       32|           30| yes |           32| None                    | tanh             |   0.0321200|           0.9947200|
| 8     | 30     |        1| valid   | Adam       |                     16|             32|       32|           30| yes |           32| None                    | tanh             |   0.0120667|           0.9945667|

#### Script 4 (Finished)

``` r
results_script4_train_val <- data.table::fread("../../../results/results/all_scan_results_0004-cnn-lstm-scan-batch32-unit32.csv")

results_script4_train_val %>% 
  select(params, mean_train_score, mean_test_score) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                                                   |  mean\_train\_score|  mean\_test\_score|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------:|------------------:|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |           0.9967480|          0.7008130|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |           0.9967480|          0.7186992|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |           0.9943089|          0.6943089|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 48, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |           0.9967480|          0.7024390|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |           0.9922764|          0.6861789|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |           0.9963415|          0.7186992|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |           0.9967480|          0.6764228|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |           0.9967480|          0.6959350|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |           0.9967480|          0.6829268|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |           0.9955285|          0.7121951|

``` r
results_script4_test <- data.table::fread("../../../results/results/df_pred_results_0004-cnn-lstm-scan-batch32-unit32.csv",  drop = 'V1')

results_script4_test %>% 
  knitr::kable()
```

| Parameters                                                                                                                                                                                                                               |   Accuracy|  Sensitivity|  Specifity|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.7333333|    0.7631579|  0.7027027|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6866667|    0.7763158|  0.5945946|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6866667|    0.8026316|  0.5675676|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 48, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7066667|    0.7500000|  0.6621622|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7066667|    0.6842105|  0.7297297|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6733333|    0.7368421|  0.6081081|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6933333|    0.7105263|  0.6756757|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6600000|    0.5921053|  0.7297297|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6666667|    0.6578947|  0.6756757|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7200000|    0.7763158|  0.6621622|

### CNN - GRU

<table style="width:100%;">
<colgroup>
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 9%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 8%" />
<col style="width: 5%" />
<col style="width: 10%" />
<col style="width: 13%" />
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 10%" />
<col style="width: 4%" />
</colgroup>
<thead>
<tr class="header">
<th>Script</th>
<th>Batch</th>
<th>GRU Hidden Unit</th>
<th>Filters</th>
<th>Kernel Size</th>
<th>Maxpool Size</th>
<th>Dropout</th>
<th>Dropout Recurrent</th>
<th>Activation Convolution</th>
<th>Go Backwards</th>
<th>Optimizers</th>
<th>Regularizer rates</th>
<th>Epochs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>8</td>
<td>16</td>
<td>[32, 48]</td>
<td>2</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>True</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="even">
<td>2</td>
<td>16</td>
<td>32</td>
<td>[32, 48]</td>
<td>2</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>True</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="odd">
<td>3</td>
<td>8</td>
<td>64</td>
<td>[32, 48]</td>
<td>2</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>True</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="even">
<td>4</td>
<td>16</td>
<td>16</td>
<td>[32, 48]</td>
<td>2</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>True</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="odd">
<td>5</td>
<td>8</td>
<td>32</td>
<td>[32, 48]</td>
<td>2</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>True</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td></td>
</tr>
<tr class="even">
<td>6</td>
<td>16</td>
<td>64</td>
<td>[32, 48]</td>
<td>2</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>True</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td></td>
</tr>
</tbody>
</table>

``` r
src_file_cnn_gru_3 <- "../../../results/results/0007-cnn-gru-grid-batch8_gru_unit64.log"
src_file_cnn_gru_4 <- "../../../results/results/0008-cnn-gru-grid-batch16_gru_unit16.log"
src_file_cnn_gru_5 <- "../../../results/results/0017-cnn-gru-grid-batch8_gru_unit32.log"
src_file_cnn_gru_6 <- "../../../results/results/0018-cnn-gru-grid-batch16_gru_unit64.log"

src_params_list_cnn_gru <- c(
  "reg_rate", "optimizers", "opt_go_backwards",
  "opt_dropout_recurrent", "opt_dropout", "maxpool_size", "kernel_size",
  "gru_hidden_units", "filter_conv", "epochs_raw", "batch_size", "activation_conv"
)
src_numeric_params_cnn_gru <- c(
  "opt_dropout_recurrent", "opt_dropout", "maxpool_size", "kernel_size",
  "gru_hidden_units", "filter_conv", "epochs_raw", "batch_size"
)
```

#### Script 1 (Pending)

This job is still pending to be run in GPU.

#### Script 2

``` r
results_cnn_gru_script6_train_val <- data.table::fread("../../../results/results/all_scan_results_0006-cnn-gru-grid-batch16_gru_unit32.csv")

results_cnn_gru_script6_train_val %>% 
  select(params, mean_train_score, mean_test_score) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                                                                                      |  mean\_train\_score|  mean\_test\_score|
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------:|------------------:|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |           0.8987805|          0.6894309|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |           0.9821138|          0.7024390|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |           0.9426829|          0.7089431|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |           0.9918699|          0.6894309|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |           0.8426829|          0.6520325|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |           0.9670732|          0.7154472|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |           0.9227642|          0.6796748|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |           0.9955285|          0.7024390|

``` r
results_cnn_gru_script6_test <- data.table::fread("../../../results/results/df_pred_results_0006-cnn-gru-grid-batch16_gru_unit32.csv",  drop = 'V1')

results_cnn_gru_script6_test %>% 
  knitr::kable()
```

| Parameters                                                                                                                                                                                                                                                                  |   Accuracy|  Sensitivity|  Specifity|
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6666667|    0.7105263|  0.6216216|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.6733333|    0.8026316|  0.5405405|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6800000|    0.6842105|  0.6756757|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.7066667|    0.6052632|  0.8108108|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.7000000|    0.7105263|  0.6891892|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.6866667|    0.7105263|  0.6621622|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7666667|    0.7236842|  0.8108108|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7200000|    0.8157895|  0.6216216|

#### Script 3

``` r
data_cnn_gru_3 <- src_file_cnn_gru_3 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_gru, src_numeric_params_cnn_gru) %>%
  clean_log_df_with_params()

new_data_cnn_gru_3 <- data_cnn_gru_3 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_gru)
```

``` r
new_data_cnn_gru_3 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate | optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  maxpool\_size|  kernel\_size|  gru\_hidden\_units|  filter\_conv|  epochs\_raw|  batch\_size| activation\_conv |  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|:-----------|:-------------------|------------------------:|-------------:|--------------:|-------------:|-------------------:|-------------:|------------:|------------:|:-----------------|-----------:|-------------------:|
| 2     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            32|           30|            8| None             |   0.0883400|           0.9939200|
| 4     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            48|           30|            8| None             |   0.0851333|           0.9925667|
| 3     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            48|           30|            8| None             |   0.2020400|           0.9711200|
| 1     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            32|           30|            8| None             |   0.2410200|           0.9540600|

#### Script 4

``` r
data_cnn_gru_4 <- src_file_cnn_gru_4 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_gru, src_numeric_params_cnn_gru) %>%
  clean_log_df_with_params()

new_data_cnn_gru_4 <- data_cnn_gru_4 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_gru)
```

``` r
new_data_cnn_gru_4 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate | optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  maxpool\_size|  kernel\_size|  gru\_hidden\_units|  filter\_conv|  epochs\_raw|  batch\_size| activation\_conv |  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|:-----------|:-------------------|------------------------:|-------------:|--------------:|-------------:|-------------------:|-------------:|------------:|------------:|:-----------------|-----------:|-------------------:|
| 4     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  16|            48|           30|           16| None             |     0.10826|             0.99308|
| 2     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  16|            32|           30|           16| None             |     0.16070|             0.96788|
| 3     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  16|            48|           30|           16| None             |     0.35550|             0.90078|
| 1     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  16|            32|           30|           16| None             |     0.40014|             0.86624|
| 5     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  16|            32|           30|           16| relu             |     0.43598|             0.83904|

#### Script 5

``` r
data_cnn_gru_5 <- src_file_cnn_gru_5 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_gru, src_numeric_params_cnn_gru) %>%
  clean_log_df_with_params()

new_data_cnn_gru_5 <- data_cnn_gru_5 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_gru)
```

``` r
new_data_cnn_gru_5 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate | optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  maxpool\_size|  kernel\_size|  gru\_hidden\_units|  filter\_conv|  epochs\_raw|  batch\_size| activation\_conv |  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|:-----------|:-------------------|------------------------:|-------------:|--------------:|-------------:|-------------------:|-------------:|------------:|------------:|:-----------------|-----------:|-------------------:|
| 4     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  32|            48|           30|            8| None             |     0.08364|             0.99510|
| 2     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  32|            32|           30|            8| None             |     0.10282|             0.98902|
| 3     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  32|            48|           30|            8| None             |     0.20220|             0.97968|
| 1     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  32|            32|           30|            8| None             |     0.28688|             0.93050|
| 5     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  32|            32|           30|            8| relu             |     0.33836|             0.89878|

#### Script 6

``` r
data_cnn_gru_6 <- src_file_cnn_gru_6 %>%
  read_log_into_df_with_params_list(src_params_list_cnn_gru, src_numeric_params_cnn_gru) %>%
  clean_log_df_with_params()

new_data_cnn_gru_6 <- data_cnn_gru_6 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_gru)
```

``` r
new_data_cnn_gru_6 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate | optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  maxpool\_size|  kernel\_size|  gru\_hidden\_units|  filter\_conv|  epochs\_raw|  batch\_size| activation\_conv |  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|:-----------|:-------------------|------------------------:|-------------:|--------------:|-------------:|-------------------:|-------------:|------------:|------------:|:-----------------|-----------:|-------------------:|
| 4     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            48|           30|           16| None             |     0.11006|             0.98942|
| 2     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            32|           30|           16| None             |     0.13074|             0.98252|
| 6     | 30     | 0.001     | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            32|           30|           16| relu             |     0.15040|             0.97760|
| 3     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            48|           30|           16| None             |     0.27600|             0.94552|
| 1     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            32|           30|           16| None             |     0.32754|             0.90976|
| 5     | 30     | 0.01      | Adam       | TRUE               |                        0|             0|              3|             2|                  64|            32|           30|           16| relu             |     0.40380|             0.85896|

### LSTM

<table style="width:100%;">
<colgroup>
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 11%" />
<col style="width: 14%" />
<col style="width: 11%" />
<col style="width: 7%" />
<col style="width: 11%" />
<col style="width: 4%" />
</colgroup>
<thead>
<tr class="header">
<th>Script</th>
<th>Batch</th>
<th>LSTM Units</th>
<th>Output Dim</th>
<th>Maxpool Size</th>
<th>Dropout</th>
<th>Dropout Recurrent</th>
<th>Activation Convolution</th>
<th>Go Backwards</th>
<th>Optimizers</th>
<th>Regularizer rates</th>
<th>Epochs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>8</td>
<td>16</td>
<td>64</td>
<td>3</td>
<td>[0, 0.25]</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>[‘True’, ‘False’]</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="even">
<td>2</td>
<td>16</td>
<td>32</td>
<td>64</td>
<td>3</td>
<td>[0, 0.25]</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>[‘True’, ‘False’]</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="odd">
<td>3</td>
<td>8</td>
<td>16</td>
<td>64</td>
<td>3</td>
<td>[0, 0.25]</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>[‘True’, ‘False’]</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
<tr class="even">
<td>4</td>
<td>16</td>
<td>32</td>
<td>64</td>
<td>3</td>
<td>[0, 0.25]</td>
<td>0</td>
<td>[‘None’, ‘relu’]</td>
<td>[‘True’, ‘False’]</td>
<td>Adam</td>
<td>[0.01, 0.001]</td>
<td>30</td>
</tr>
</tbody>
</table>

``` r
src_file_lstm_1 <- "../../../results/results/0013-lstm-embedding-grid-batch8_lstm_unit16.log"
src_file_lstm_2 <- "../../../results/results/0014-lstm-embedding-grid-batch16_lstm_unit16_old.log"
src_file_lstm_3 <- "../../../results/results/0015-lstm-embedding-grid-batch8_lstm_unit32.log"
src_file_lstm_4 <- "../../../results/results/0016-lstm-embedding-grid-batch16_lstm_unit32.log"

src_params_list_lstm <- c(
  "reg_rate",
  "outputdim", "optimizers", 
  "opt_go_backwards", "opt_dropout_recurrent", "opt_dropout",
  "lstm_hidden_units", "epochs_raw", "batch_size"
)
src_numeric_params_lstm <- c(
  "outputdim", "opt_dropout_recurrent", "opt_dropout",
  "lstm_hidden_units", "epochs_raw", "batch_size"
)
```

#### Script 1

``` r
data_lstm_1 <- src_file_lstm_1 %>%
  read_log_into_df_with_params_list(src_params_list_lstm, src_numeric_params_lstm) %>%
  clean_log_df_with_params()

new_data_lstm_1 <- data_lstm_1 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_lstm)
```

``` r
new_data_lstm_1 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate |  outputdim| optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  lstm\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|----------:|:-----------|:-------------------|------------------------:|-------------:|--------------------:|------------:|------------:|-----------:|-------------------:|
| 2     | 30     | 0.001     |         64| Adam       | TRUE               |                        0|             0|                   16|           30|            8|   0.4427400|             0.81258|
| 1     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|             0|                   16|           30|            8|   0.4963000|             0.79146|
| 3     | 30     | 0.01      |         64| Adam       | FALSE              |                        0|             0|                   16|           30|            8|   0.4887667|             0.78650|

#### Script 2

``` r
data_lstm_2 <- src_file_lstm_2 %>%
  read_log_into_df_with_params_list(src_params_list_lstm, src_numeric_params_lstm) %>%
  clean_log_df_with_params()

new_data_lstm_2 <- data_lstm_2 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_lstm)
```

``` r
new_data_lstm_2 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate |  outputdim| optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  lstm\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|----------:|:-----------|:-------------------|------------------------:|-------------:|--------------------:|------------:|------------:|-----------:|-------------------:|
| 4     | 30     | 0.001     |         64| Adam       | FALSE              |                        0|          0.00|                   16|           30|           16|     0.44626|             0.81544|
| 2     | 30     | 0.001     |         64| Adam       | TRUE               |                        0|          0.00|                   16|           30|           16|     0.46294|             0.80812|
| 8     | 30     | 0.001     |         64| Adam       | FALSE              |                        0|          0.25|                   16|           30|           16|     0.46902|             0.79226|
| 6     | 30     | 0.001     |         64| Adam       | TRUE               |                        0|          0.25|                   16|           30|           16|     0.48498|             0.78902|
| 1     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|          0.00|                   16|           30|           16|     0.50012|             0.78820|
| 7     | 30     | 0.01      |         64| Adam       | FALSE              |                        0|          0.25|                   16|           30|           16|     0.49880|             0.78782|
| 3     | 30     | 0.01      |         64| Adam       | FALSE              |                        0|          0.00|                   16|           30|           16|     0.49894|             0.78456|
| 5     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|          0.25|                   16|           30|           16|     0.50266|             0.78252|
| 9     | 30     | 0.001     |         64| Adam       | FALSE              |                        0|          0.25|                   16|           30|           16|     0.50050|             0.78210|

#### Script 3

``` r
data_lstm_3 <- src_file_lstm_3 %>%
  read_log_into_df_with_params_list(src_params_list_lstm, src_numeric_params_lstm) %>%
  clean_log_df_with_params()

new_data_lstm_3 <- data_lstm_3 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_lstm)
```

``` r
new_data_lstm_3 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate |  outputdim| optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  lstm\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|----------:|:-----------|:-------------------|------------------------:|-------------:|--------------------:|------------:|------------:|-----------:|-------------------:|
| 2     | 30     | 0.001     |         64| Adam       | TRUE               |                        0|          0.00|                   32|           30|            8|     0.46420|             0.80488|
| 4     | 30     | 0.001     |         64| Adam       | FALSE              |                        0|          0.00|                   32|           30|            8|     0.46580|             0.80448|
| 6     | 30     | 0.001     |         64| Adam       | TRUE               |                        0|          0.25|                   32|           30|            8|     0.47426|             0.79754|
| 8     | 30     | 0.001     |         64| Adam       | FALSE              |                        0|          0.25|                   32|           30|            8|     0.48798|             0.79754|
| 3     | 30     | 0.01      |         64| Adam       | FALSE              |                        0|          0.00|                   32|           30|            8|     0.49644|             0.79350|
| 7     | 30     | 0.01      |         64| Adam       | FALSE              |                        0|          0.25|                   32|           30|            8|     0.50080|             0.78498|
| 1     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|          0.00|                   32|           30|            8|     0.49830|             0.78294|
| 5     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|          0.25|                   32|           30|            8|     0.50446|             0.78010|

#### Script 4

``` r
data_lstm_4 <- src_file_lstm_4 %>%
  read_log_into_df_with_params_list(src_params_list_lstm, src_numeric_params_lstm) %>%
  clean_log_df_with_params()

new_data_lstm_4 <- data_lstm_4 %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_lstm)
```

``` r
new_data_lstm_4 %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate |  outputdim| optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  lstm\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|----------:|:-----------|:-------------------|------------------------:|-------------:|--------------------:|------------:|------------:|-----------:|-------------------:|
| 4     | 30     | 0.001     |         64| Adam       | FALSE              |                        0|          0.00|                   32|           30|           16|    0.469700|             0.79798|
| 1     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|          0.00|                   32|           30|           16|    0.500520|             0.78660|
| 3     | 30     | 0.01      |         64| Adam       | FALSE              |                        0|          0.00|                   32|           30|           16|    0.496940|             0.78052|
| 5     | 30     | 0.01      |         64| Adam       | TRUE               |                        0|          0.25|                   32|           30|           16|    0.506925|             0.77490|
| 2     | 30     | 0.001     |         64| Adam       | TRUE               |                        0|          0.00|                   32|           30|           16|    0.499360|             0.77396|

### GRU (RandomGridSearchCV)

``` r
src_file_gru <- "../../../results/results/0010_gru-embedding_scan.log"

src_params_list_gru <- c(
  "reg_rate", "outputdim", "optimizers", 
  "opt_go_backwards", "opt_dropout_recurrent", "opt_dropout", "gru_hidden_units",
  "epochs_raw", "batch_size"
)
src_numeric_params_gru <- c(
  "outputdim", "opt_dropout_recurrent", "opt_dropout",
  "gru_hidden_units", "epochs_raw", "batch_size"
)
```

``` r
data_gru <- src_file_gru %>%
  read_log_into_df_with_params_list(src_params_list_gru, src_numeric_params_gru) %>%
  clean_log_df_with_params()

new_data_gru <- data_gru %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_gru)
```

``` r
new_data_gru %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  knitr::kable()
```

| model | epochs | reg\_rate |  outputdim| optimizers | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  gru\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  mean\_train\_score|
|:------|:-------|:----------|----------:|:-----------|:-------------------|------------------------:|-------------:|-------------------:|------------:|------------:|-----------:|-------------------:|
| 13    | 30     | 0.001     |         32| Adam       | FALSE              |                     0.25|          0.50|                  32|           30|           16|     0.54448|             0.75204|
| 7     | 30     | 0.001     |         64| Adadelta   | TRUE               |                     0.50|          0.25|                  32|           30|           16|     0.56978|             0.72316|
| 14    | 30     | 0.001     |         48| Adadelta   | TRUE               |                     0.50|          0.25|                  64|           30|           32|     0.61692|             0.70892|
| 4     | 30     | 0.001     |         64| Adam       | TRUE               |                     0.50|          0.00|                   8|           30|           16|     0.58350|             0.70244|
| 9     | 30     | 0.01      |         64| Adam       | TRUE               |                     0.50|          0.00|                   8|           30|           32|     0.61362|             0.68740|
| 1     | 30     | 0.01      |         48| Adam       | FALSE              |                     0.50|          0.00|                   8|           30|           16|     0.59968|             0.68496|
| 12    | 30     | 0.01      |         48| Adadelta   | TRUE               |                     0.25|          0.25|                   8|           30|           32|     0.60800|             0.68454|
| 8     | 30     | 0.01      |         64| Adadelta   | FALSE              |                     0.50|          0.25|                  16|           30|           32|     0.60732|             0.68172|
| 5     | 30     | 0.01      |         48| Adadelta   | TRUE               |                     0.50|          0.25|                  32|           30|           32|     0.60434|             0.68092|
| 15    | 30     | 0.01      |         64| Adam       | FALSE              |                     0.50|          0.50|                  32|           30|           16|     0.60325|             0.67785|
| 3     | 30     | 0.01      |         32| sgd        | FALSE              |                     0.00|          0.00|                   8|           30|           16|     1.18326|             0.57966|
| 10    | 30     | 0.001     |         32| sgd        | FALSE              |                     0.00|          0.00|                  16|           30|           16|     0.79496|             0.57072|
| 11    | 30     | 0.001     |         64| sgd        | TRUE               |                     0.00|          0.25|                  16|           30|           16|     0.82488|             0.56828|
| 6     | 30     | 0.001     |         48| sgd        | FALSE              |                     0.00|          0.00|                  64|           30|           16|     0.96220|             0.55934|
| 2     | 30     | 0.001     |         48| sgd        | FALSE              |                     0.00|          0.00|                  32|           30|           16|     0.87576|             0.54636|
