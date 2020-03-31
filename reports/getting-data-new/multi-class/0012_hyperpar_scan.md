Getting the results of hyperparameter scans
===========================================

CNN-GRU
-------

``` r
# Load library
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
# Get the results of Bacteria
df_result_bacteria <- data.table::fread("../../../../results/scan_separate_multi_class/bacteria/df_pred_results_0002-cnn-gru-grid-bacteria.csv")

df_result_bacteria %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable() 
```

|   V1| Parameters                                                                                                                                                                                                                                                                |   Accuracy|  Sensitivity|  Specifity|
|----:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}   |  0.9078947|    0.8947368|  0.9210526|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}   |  0.8947368|    0.8421053|  0.9473684|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.8947368|    0.8421053|  0.9473684|
|    8| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’} |  0.8947368|    0.8421053|  0.9473684|
|    9| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’} |  0.8947368|    0.8157895|  0.9736842|
|   12| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}     |  0.8947368|    0.8421053|  0.9473684|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}   |  0.8815790|    0.8684211|  0.8947368|
|   13| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.8815790|    0.8684211|  0.8947368|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.8815789|    0.8684211|  0.8947368|
|   14| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.8684211|    0.8684211|  0.8684211|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.8684211|    0.8157895|  0.9210526|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}   |  0.8684211|    0.8684211|  0.8684211|
|   10| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.8684211|    0.8421053|  0.8947368|
|   15| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.8421053|    0.7368421|  0.9473684|
|   11| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.7763158|    0.6578947|  0.8947368|
|    6| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.7631579|    0.7368421|  0.7894737|

``` r
# Get the results of Bacteria
df_result_fungi <- data.table::fread("../../../../results/scan_separate_multi_class/fungi/df_pred_results_0002-cnn-gru-grid-fungi.csv")

df_result_fungi %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable() 
```

|   V1| Parameters                                                                                                                                                                                                                                                                |   Accuracy|  Sensitivity|  Specifity|
|----:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    7| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7045455|    0.7727273|  0.6363636|
|   14| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}   |  0.7045455|    0.9090909|  0.5000000|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}   |  0.6818182|    0.8181818|  0.5454545|
|   10| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.6818182|    0.7727273|  0.5909091|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’} |  0.6818182|    0.7727273|  0.5909091|
|    8| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.6818182|    0.7727273|  0.5909091|
|    4| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.6590909|    0.7727273|  0.5454545|
|    0| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.6363636|    0.7272727|  0.5454545|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.6363636|    0.7272727|  0.5454545|
|    9| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.6363636|    0.7727273|  0.5000000|
|   11| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}     |  0.6363636|    0.7272727|  0.5454545|
|   12| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.6363636|    0.8181818|  0.4545455|
|   13| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}   |  0.6363636|    0.7727273|  0.5000000|
|    3| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}   |  0.6136364|    0.6818182|  0.5454545|
|   15| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}   |  0.6136364|    0.6818182|  0.5454545|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.5909091|    0.7272727|  0.4545455|

``` r
# Get the results of Bacteria
df_result_oomycete <- data.table::fread("../../../../results/scan_separate_multi_class/oomycete/df_pred_results_0002-cnn-gru-grid-oomycete.csv")

df_result_oomycete %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable() 
```

|   V1| Parameters                                                                                                                                                                                                                                                                 |   Accuracy|  Sensitivity|  Specifity|
|----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.8157895|    0.7894737|  0.8421053|
|    8| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}    |  0.7894737|    0.7894737|  0.7894737|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.7894737|    0.7894737|  0.7894737|
|   14| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}   |  0.7368421|    0.7368421|  0.7368421|
|    0| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.7368421|    0.6842105|  0.7894737|
|    3| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}      |  0.7368421|    0.7368421|  0.7368421|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.7368421|    0.7368421|  0.7368421|
|   10| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.7368421|    0.7368421|  0.7368421|
|   12| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7368421|    0.7368421|  0.7368421|
|    1| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7105263|    0.6842105|  0.7368421|
|   13| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.7105263|    0.6315789|  0.7894737|
|   15| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}      |  0.7105263|    0.7368421|  0.6842105|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}     |  0.7105263|    0.7894737|  0.6315789|
|    9| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}     |  0.7105263|    0.6842105|  0.7368421|
|    6| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.6842105|    0.8947368|  0.4736842|
|   11| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.6315790|    0.7368421|  0.5263158|

CNN-LSTM
--------

### CNN-LSTM Bacteria

``` r
src_file_cnn_lstm_bacteria <- "../../../../results/scan_separate_multi_class/bacteria/0001-cnn-lstm-scan-bacteria.log"

src_params_list_cnn_lstm_bacteria <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params_cnn_lstm_bacteria <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)
```

``` r
data_cnn_lstm_bacteria <- src_file_cnn_lstm_bacteria %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm_bacteria, src_numeric_params_cnn_lstm_bacteria) %>%
  clean_log_df_with_params()
```

    ## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
    ## This warning is displayed once per session.

``` r
new_data_cnn_lstm_bacteria <- data_cnn_lstm_bacteria %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm_bacteria)
```

``` r
new_data_cnn_lstm_bacteria %>% 
  dplyr::select(-c(loss_sd, acc_sd)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  arrange(desc(mean_train_score)) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw|  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  mean\_train\_score|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|------------:|:------------------------|:-----------------|-----------:|-------------------:|
| 1     | 30     |        1| valid   | Adadelta   |                      8|             16|        8|           30|           16| None                    | tanh             |   0.3513000|             0.90624|
| 3     | 30     |        1| valid   | Adam       |                      8|             16|        8|           30|           16| None                    | tanh             |   0.3332600|             0.88898|
| 4     | 30     |        1| valid   | Adadelta   |                      4|              8|        8|           30|            8| None                    | tanh             |   0.4660600|             0.83388|
| 8     | 30     |        1| valid   | Adam       |                      4|              8|        8|           30|            8| None                    | tanh             |   0.4552667|             0.81480|
| 7     | 30     |        1| valid   | Adam       |                      4|              4|        8|           30|            8| None                    | tanh             |   0.5069600|             0.77628|
| 2     | 30     |        1| valid   | Adam       |                      4|              8|        4|           30|           16| None                    | tanh             |   0.5773600|             0.71868|
| 6     | 30     |        1| valid   | Adam       |                      8|              4|        4|           30|           16| None                    | tanh             |   0.6080800|             0.70564|
| 5     | 30     |        1| valid   | Adadelta   |                      4|              4|        4|           30|            8| None                    | tanh             |   0.6330600|             0.68594|

### CNN-LSTM Fungi

``` r
src_file_cnn_lstm_fungi <- "../../../../results/scan_separate_multi_class/fungi/0001-cnn-lstm-scan-fungi.log"

src_params_list_cnn_lstm_fungi <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params_cnn_lstm_fungi <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)
```

``` r
data_cnn_lstm_fungi <- src_file_cnn_lstm_fungi %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm_fungi, src_numeric_params_cnn_lstm_fungi) %>%
  clean_log_df_with_params()

new_data_cnn_lstm_fungi <- data_cnn_lstm_fungi %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm_fungi)
```

``` r
new_data_cnn_lstm_fungi %>% 
  dplyr::select(-c(loss_sd, acc_sd)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  arrange(desc(mean_train_score)) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw|  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  mean\_train\_score|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|------------:|:------------------------|:-----------------|-----------:|-------------------:|
| 7     | 30     |        1| valid   | Adam       |                      8|              8|        8|           30|            8| None                    | tanh             |     0.40808|             0.89304|
| 3     | 30     |        1| valid   | Adadelta   |                      8|              8|        8|           30|            8| None                    | tanh             |     0.36230|             0.88612|
| 9     | 30     |        1| valid   | Adam       |                      8|              8|        4|           30|            8| None                    | tanh             |     0.43034|             0.84306|
| 8     | 30     |        1| valid   | Adadelta   |                      8|              8|        4|           30|            8| None                    | tanh             |     0.45976|             0.82084|
| 5     | 30     |        1| valid   | Adadelta   |                      8|              8|        4|           30|           16| None                    | tanh             |     0.52928|             0.78056|
| 4     | 30     |        1| valid   | Adam       |                      4|              8|        4|           30|            8| None                    | tanh             |     0.50104|             0.77084|
| 6     | 30     |        1| valid   | Adadelta   |                      4|              8|        4|           30|            8| None                    | tanh             |     0.53002|             0.73472|
| 12    | 30     |        1| valid   | Adam       |                      4|              4|        4|           30|            8| None                    | tanh             |     0.58328|             0.71032|
| 11    | 30     |        1| valid   | Adam       |                      4|              8|        4|           30|           16| None                    | tanh             |     0.58912|             0.70972|
| 2     | 30     |        1| valid   | Adadelta   |                      4|              4|        4|           30|           16| None                    | tanh             |     0.62322|             0.70694|
| 10    | 30     |        1| valid   | Adam       |                      4|              4|        4|           30|           16| None                    | tanh             |     0.63926|             0.69028|
| 1     | 30     |        1| valid   | Adadelta   |                      4|              4|        4|           30|            8| None                    | tanh             |     0.60214|             0.65830|

### CNN-LSTM Oomycete

``` r
src_file_cnn_lstm_oomycete <- "../../../../results/scan_separate_multi_class/oomycete/0001-cnn-lstm-scan-oomycete.log"

src_params_list_cnn_lstm_oomycete <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params_cnn_lstm_oomycete <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)
```

``` r
data_cnn_lstm_oomycete <- src_file_cnn_lstm_oomycete %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm_oomycete, src_numeric_params_cnn_lstm_oomycete) %>%
  clean_log_df_with_params()

new_data_cnn_lstm_oomycete <- data_cnn_lstm_oomycete %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm_oomycete)
```

``` r
new_data_cnn_lstm_oomycete %>% 
  dplyr::select(-c(loss_sd, acc_sd)) %>% 
  dplyr::rename(mean_train_score = acc_mean) %>% 
  arrange(desc(mean_train_score)) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw|  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  mean\_train\_score|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|------------:|:------------------------|:-----------------|-----------:|-------------------:|
| 7     | 30     |        1| valid   | Adam       |                      8|              8|        8|           30|           16| None                    | tanh             |     0.44610|             0.88430|
| 6     | 30     |        1| valid   | Adam       |                      8|              8|        8|           30|            8| None                    | tanh             |     0.38216|             0.88312|
| 1     | 30     |        1| valid   | Adam       |                      4|              8|        8|           30|            8| None                    | tanh             |     0.51726|             0.79616|
| 4     | 30     |        1| valid   | Adam       |                      4|              4|        8|           30|            8| None                    | tanh             |     0.56448|             0.78130|
| 3     | 30     |        1| valid   | Adadelta   |                      8|              8|        4|           30|           16| None                    | tanh             |     0.56136|             0.76156|
| 2     | 30     |        1| valid   | Adadelta   |                      4|              8|        4|           30|            8| None                    | tanh             |     0.57622|             0.70550|
| 5     | 30     |        1| valid   | Adam       |                      4|              8|        4|           30|           16| None                    | tanh             |     0.63182|             0.68250|
