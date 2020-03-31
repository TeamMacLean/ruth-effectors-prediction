Results of hyperoar-scan
========================

CNN - LSTM
----------

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
result1 <- data.table::fread("../../../../results/results/df_pred_results_0002-cnn-lstm-scan-batch16-unit32.csv")
result2 <- data.table::fread("../../../../results/results/df_pred_results_0003-cnn-lstm-scan-batch32-unit16.csv")
result3 <- data.table::fread("../../../../results/results/df_pred_results_0004-cnn-lstm-scan-batch32-unit32.csv")
```

``` r
result1 %>% 
  rbind(., result2) %>% 
  rbind(., result3) %>%
  arrange(desc(Accuracy)) %>%   
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                               |   Accuracy|  Sensitivity|  Specifity|
|----:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    0| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.7333333|    0.7631579|  0.7027027|
|    3| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.7200000|    0.7500000|  0.6891892|
|    9| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7200000|    0.7763158|  0.6621622|
|    1| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 48, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7133333|    0.7763158|  0.6486486|
|    3| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 48, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.7133333|    0.7500000|  0.6756757|
|    0| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 32, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.7133333|    0.7763158|  0.6486486|
|    9| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.7133333|    0.7500000|  0.6756757|
|    6| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7133333|    0.8157895|  0.6081081|
|    3| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 48, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7066667|    0.7500000|  0.6621622|
|    4| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7066667|    0.6842105|  0.7297297|
|    6| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6933333|    0.7105263|  0.6756757|
|    7| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6933333|    0.7105263|  0.6756757|
|    0| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 48, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6933333|    0.7500000|  0.6351351|
|    5| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6866667|    0.6710526|  0.7027027|
|    1| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 48, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6866667|    0.7631579|  0.6081081|
|    2| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6866667|    0.8026316|  0.5675676|
|    1| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6866667|    0.7763158|  0.5945946|
|    4| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6800000|    0.7368421|  0.6216216|
|    6| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6800000|    0.7236842|  0.6351351|
|    5| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6733333|    0.7368421|  0.6081081|
|    4| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6733333|    0.7105263|  0.6351351|
|    7| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6733333|    0.6447368|  0.7027027|
|    9| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6666667|    0.7105263|  0.6216216|
|    2| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 32, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6666667|    0.7894737|  0.5405405|
|    8| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 16, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6666667|    0.6578947|  0.6756757|
|    7| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6600000|    0.5921053|  0.7297297|
|    8| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6600000|    0.6973684|  0.6216216|
|    5| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 16, ‘filters\_LSTM’: 8, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 32, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6600000|    0.7236842|  0.5945946|
|    2| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 32, ‘filters’: 32, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6466667|    0.6578947|  0.6351351|
|    8| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 32, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘bn’: ‘yes’, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6200000|    0.6052632|  0.6351351|

CNN-GRU
-------

``` r
res_cnn_gru_1 <- data.table::fread("../../../../results/results/df_pred_results_0005-cnn-gru-grid-batch8_gru_unit16.csv")
res_cnn_gru_2 <- data.table::fread("../../../../results/results/df_pred_results_0006-cnn-gru-grid-batch16_gru_unit32.csv")
res_cnn_gru_3 <- data.table::fread("../../../../results/results/df_pred_results_0007-cnn-gru-grid-batch8_gru_unit64.csv")
res_cnn_gru_4 <- data.table::fread("../../../../results/results/df_pred_results_0008-cnn-gru-grid-batch16_gru_unit16.csv")
res_cnn_gru_5 <- data.table::fread("../../../../results/results/df_pred_results_0017-cnn-gru-grid-batch8_gru_unit32.csv")
res_cnn_gru_6 <- data.table::fread("../../../../results/results/df_pred_results_0018-cnn-gru-grid-batch16_gru_unit64.csv")
```

``` r
res_cnn_gru_1 %>% 
  rbind(., res_cnn_gru_2) %>% 
  rbind(., res_cnn_gru_3) %>% 
  rbind(., res_cnn_gru_4) %>% 
  rbind(., res_cnn_gru_5) %>% 
  rbind(., res_cnn_gru_6) %>% 
  arrange(desc(Accuracy)) %>%   
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                  |   Accuracy|  Sensitivity|  Specifity|
|----:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7666667|    0.7236842|  0.8108108|
|    2| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7400000|    0.7500000|  0.7297297|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7400000|    0.7105263|  0.7702703|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7400000|    0.7105263|  0.7702703|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7333333|    0.7368421|  0.7297297|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.7333333|    0.6578947|  0.8108108|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7333333|    0.6710526|  0.7972973|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7266667|    0.6973684|  0.7567568|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.7200000|    0.7105263|  0.7297297|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7200000|    0.7631579|  0.6756757|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.7200000|    0.6973684|  0.7432432|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7200000|    0.8157895|  0.6216216|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7200000|    0.6052632|  0.8378378|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7200000|    0.8815789|  0.5540541|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7133333|    0.5657895|  0.8648649|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7133333|    0.7105263|  0.7162162|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.7066667|    0.7631579|  0.6486486|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7066667|    0.7631579|  0.6486486|
|    3| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7066667|    0.6315789|  0.7837838|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.7066667|    0.6052632|  0.8108108|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7066667|    0.6578947|  0.7567568|
|    4| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.7000000|    0.7105263|  0.6891892|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7000000|    0.6447368|  0.7567568|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.7000000|    0.8552632|  0.5405405|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7000000|    0.6973684|  0.7027027|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.6933333|    0.6973684|  0.6891892|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.6933333|    0.7368421|  0.6486486|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.6933333|    0.7236842|  0.6621622|
|    4| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.6933333|    0.7631579|  0.6216216|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.6866667|    0.5921053|  0.7837838|
|    5| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.6866667|    0.7105263|  0.6621622|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.6866667|    0.6973684|  0.6756757|
|    3| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.6866667|    0.7500000|  0.6216216|
|    4| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.6866667|    0.6578947|  0.7162162|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.6866667|    0.7105263|  0.6621622|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6800000|    0.6842105|  0.6756757|
|    1| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.6733333|    0.6710526|  0.6756757|
|    1| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.6733333|    0.8026316|  0.5405405|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6733333|    0.5526316|  0.7972973|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.6666667|    0.7236842|  0.6081081|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6666667|    0.7105263|  0.6216216|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.6600000|    0.6710526|  0.6486486|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6600000|    0.6315789|  0.6891892|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6600000|    0.6578947|  0.6621622|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6466667|    0.6578947|  0.6351351|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.6333333|    0.6315789|  0.6351351|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 64, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.6266667|    0.6578947|  0.5945946|
|    5| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}  |  0.6000000|    0.6052632|  0.5945946|

``` r
res_cnn_gru_2 %>% 
  arrange(desc(Accuracy)) %>%   
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                  |   Accuracy|  Sensitivity|  Specifity|
|----:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7666667|    0.7236842|  0.8108108|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7200000|    0.8157895|  0.6216216|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.7066667|    0.6052632|  0.8108108|
|    4| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.7000000|    0.7105263|  0.6891892|
|    5| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.6866667|    0.7105263|  0.6621622|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6800000|    0.6842105|  0.6756757|
|    1| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 48, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.6733333|    0.8026316|  0.5405405|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.6666667|    0.7105263|  0.6216216|

LSTM-Embedding
--------------

``` r
res_lstm_emd_1 <- data.table::fread("../../../../results/results/lstm-embedding/df_pred_results_0013-lstm-embedding-grid-batch8_lstm_unit16.csv")
res_lstm_emd_2 <- data.table::fread("../../../../results/results/lstm-embedding/df_pred_results_0014-lstm-embedding-grid-batch16_lstm_unit16.csv")
```

``` r
res_lstm_emd_1 %>% 
  rbind(., res_lstm_emd_2) %>% 
  arrange(desc(Accuracy)) %>%   
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                 |   Accuracy|  Sensitivity|  Specifity|
|----:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    6| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.7466667|    0.7894737|  0.7027027|
|    1| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.7266667|    0.6842105|  0.7702703|
|    5| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}       |  0.7200000|    0.7500000|  0.6891892|
|    1| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.7133333|    0.7763158|  0.6486486|
|    0| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |  0.7133333|    0.8289474|  0.5945946|
|    6| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |  0.7133333|    0.7763158|  0.6486486|
|    3| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}    |  0.7133333|    0.6842105|  0.7432432|
|    4| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |  0.7066667|    0.7631579|  0.6486486|
|    4| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.7000000|    0.7500000|  0.6486486|
|    5| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}   |  0.7000000|    0.7631579|  0.6351351|
|    7| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16} |  0.7000000|    0.7368421|  0.6621622|
|    2| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |  0.6866667|    0.7500000|  0.6216216|
|    2| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.6866667|    0.8157895|  0.5540541|
|    3| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}  |  0.6866667|    0.7105263|  0.6621622|
|    0| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.6733333|    0.7236842|  0.6216216|
|    7| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |  0.6533333|    0.7631579|  0.5405405|

GRU-Embedding
-------------

``` r
gru_emb_res <- data.table::fread("../../../../results/results/gru-embedding/df_pred_results_0010-gru-embedding_scan.csv")
```

``` r
gru_emb_res %>% 
  arrange(desc(Accuracy)) %>%   
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                      |   Accuracy|  Sensitivity|  Specifity|
|----:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    0| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 32}       |  0.7066667|    0.7631579|  0.6486486|
|    1| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}   |  0.7000000|    0.7894737|  0.6081081|
|    8| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}       |  0.7000000|    0.7105263|  0.6891892|
|   26| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}   |  0.6866667|    0.7236842|  0.6486486|
|   20| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |  0.6866667|    0.7105263|  0.6621622|
|   24| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}      |  0.6800000|    0.7105263|  0.6486486|
|   18| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 16}        |  0.6733333|    0.7105263|  0.6351351|
|   27| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 32}        |  0.6733333|    0.6973684|  0.6486486|
|   29| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}        |  0.6666667|    0.6710526|  0.6621622|
|   31| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}          |  0.6666667|    0.7105263|  0.6216216|
|    5| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 32}    |  0.6666667|    0.6052632|  0.7297297|
|   33| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}          |  0.6666667|    0.7105263|  0.6216216|
|   10| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}    |  0.6666667|    0.7236842|  0.6081081|
|    6| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.6600000|    0.7500000|  0.5675676|
|   14| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16} |  0.6600000|    0.7368421|  0.5810811|
|   19| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}         |  0.6600000|    0.6710526|  0.6486486|
|    3| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.6533333|    0.8157895|  0.4864865|
|   12| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}   |  0.6533333|    0.7631579|  0.5405405|
|   13| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 32}  |  0.6533333|    0.7368421|  0.5675676|
|   30| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}         |  0.6533333|    0.7105263|  0.5945946|
|    7| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.6533333|    0.6842105|  0.6216216|
|   42| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}         |  0.6466667|    0.4210526|  0.8783784|
|   17| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}       |  0.6466667|    0.6052632|  0.6891892|
|    9| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.6466667|    0.6447368|  0.6486486|
|   22| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 32}      |  0.6400000|    0.7236842|  0.5540541|
|   11| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}   |  0.6400000|    0.7236842|  0.5540541|
|   23| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.6400000|    0.7763158|  0.5000000|
|   21| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}          |  0.6400000|    0.6052632|  0.6756757|
|   32| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 32}           |  0.6333333|    0.6447368|  0.6216216|
|    2| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.6333333|    0.7236842|  0.5405405|
|    4| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.6333333|    0.7105263|  0.5540541|
|   15| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 32}      |  0.6333333|    0.6315789|  0.6351351|
|   28| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 32}   |  0.6266667|    0.5657895|  0.6891892|
|   47| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}        |  0.6266667|    0.8421053|  0.4054054|
|   16| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 16} |  0.6200000|    0.8026316|  0.4324324|
|   34| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 32}         |  0.6200000|    0.5921053|  0.6486486|
|   38| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}        |  0.6133333|    0.5394737|  0.6891892|
|   45| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}           |  0.6066667|    0.6973684|  0.5135135|
|   44| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 16}           |  0.6000000|    0.6973684|  0.5000000|
|   40| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}           |  0.5933333|    0.4605263|  0.7297297|
|   41| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}       |  0.5933333|    0.7105263|  0.4729730|
|   25| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}    |  0.5933333|    0.8026316|  0.3783784|
|   37| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}          |  0.5866667|    0.7500000|  0.4189189|
|   43| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}             |  0.5733333|    0.4210526|  0.7297297|
|   46| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}       |  0.5733333|    0.4210526|  0.7297297|
|   48| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}         |  0.5733333|    0.8026316|  0.3378378|
|   35| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 16}           |  0.5666667|    0.3552632|  0.7837838|
|   36| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 32}      |  0.5333333|    0.6842105|  0.3783784|
|   39| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 16}          |  0.5266667|    0.3289474|  0.7297297|
|   49| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}        |  0.5200000|    0.9473684|  0.0810811|
