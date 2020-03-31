Results of Hyperscan of the Effector and Non-effector (All secreted data) on different architecture
===================================================================================================

Background
----------

In this report, all differenct results in hyperparameter scan for
several different models for different datasets will be shown.

``` r
# Load libabry needed
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

Oomycete
--------

### CNN-GRU

The memory used in cluster:

    Job ID: 25614469
    Cluster: ciscluster
    User/Group: kristian/TSL_2t
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 13:18:32
    CPU Efficiency: 99.91% of 13:19:17 core-walltime
    Job Wall-clock time: 13:19:17
    Memory Utilized: 8.71 GB
    Memory Efficiency: 21.79% of 40.00 GB

``` r
oomycete_cnn_gru_secreted <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/df_pred_results_cnn_gru_scan_oomycete_secreted.csv")

oomycete_cnn_gru_all_scan_result <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/all_scan_results_cnn_gru_scan_oomycete_secreted.csv") %>% 
  dplyr::select(c(params, mean_test_score))

oomycete_cnn_gru_secreted %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                      |   Accuracy|  Sensitivity|  Specifity|
|----:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}  |  0.7647059|    0.8823529|  0.6470588|
|   28| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}        |  0.7647059|    0.9411765|  0.5882353|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}  |  0.7352941|    0.8235294|  0.6470588|
|    6| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}  |  0.7352941|    0.8823529|  0.5882353|
|    8| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |  0.7352941|    0.7647059|  0.7058824|
|   18| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7352941|    0.9411765|  0.5294118|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7058824|    0.8235294|  0.5882353|
|   10| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}      |  0.7058824|    0.8235294|  0.5882353|
|   16| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.7058824|    0.8235294|  0.5882353|
|   17| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}       |  0.7058824|    0.7058824|  0.7058824|
|   21| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}       |  0.7058824|    0.7058824|  0.7058824|
|   22| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}        |  0.7058824|    0.8235294|  0.5882353|
|   27| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.7058824|    0.7058824|  0.7058824|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |  0.6764706|    0.7647059|  0.5882353|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}   |  0.6764706|    0.6470588|  0.7058824|
|   12| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}        |  0.6764706|    0.8823529|  0.4705882|
|   15| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}        |  0.6764706|    0.6470588|  0.7058824|
|   19| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}       |  0.6764706|    0.8823529|  0.4705882|
|   20| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.6764706|    0.5882353|  0.7647059|
|   23| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}     |  0.6764706|    0.7647059|  0.5882353|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |  0.6470588|    0.6470588|  0.6470588|
|   14| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  0.6470588|    0.6470588|  0.6470588|
|   29| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.6470588|    0.6470588|  0.6470588|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |  0.6176471|    0.6470588|  0.5882353|
|    9| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  0.6176471|    0.6470588|  0.5882353|
|   13| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}  |  0.6176471|    0.6470588|  0.5882353|
|   24| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  0.5882353|    0.5294118|  0.6470588|
|   25| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.5882353|    0.5882353|  0.5882353|
|   11| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}     |  0.5588235|    0.5882353|  0.5294118|
|   26| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.5588235|    0.5882353|  0.5294118|

``` r
oomycete_cnn_gru_all_scan_result %>% 
  arrange(desc(mean_test_score)) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                                                                                          |  mean\_test\_score|
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |          0.6911765|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}  |          0.6838235|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |          0.6838235|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |          0.6764706|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}  |          0.6764706|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}  |          0.6764706|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}  |          0.6617647|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}   |          0.6617647|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |          0.6617647|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |          0.6544118|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}      |          0.6470588|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}     |          0.6470588|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}        |          0.6397059|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}  |          0.6323529|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |          0.6323529|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}        |          0.6323529|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |          0.6323529|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}       |          0.6250000|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |          0.6176471|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}       |          0.6102941|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |          0.6029412|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}       |          0.5955882|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}        |          0.5955882|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}     |          0.5882353|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |          0.5882353|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |          0.5882353|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |          0.5808824|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |          0.5808824|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}        |          0.5735294|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |          0.5661765|

### GRU-Embedding

``` r
oomycete_gru_emb_results_all <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/oomycete_scan_results_gru_embedding_scan_oomycete.csv") %>% 
  dplyr::select(c(params, mean_test_score))


oomycete_gru_emb_results <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/df_pred_results_gru_embedding_scan_oomycete.csv")
```

``` r
oomycete_gru_emb_results %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                    |   Accuracy|  Sensitivity|  Specifity|
|----:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|   32| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}      |  0.7058824|    0.7058824|  0.7058824|
|   19| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.6764706|    0.7647059|  0.5882353|
|   21| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 4}        |  0.6764706|    0.7647059|  0.5882353|
|   27| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}       |  0.6764706|    0.7058824|  0.6470588|
|   42| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}        |  0.6764706|    0.5882353|  0.7647059|
|   43| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}           |  0.6764706|    0.5294118|  0.8235294|
|    0| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}  |  0.6470588|    0.7647059|  0.5294118|
|    1| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 8}    |  0.6470588|    0.7058824|  0.5882353|
|    2| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 4}  |  0.6470588|    0.7058824|  0.5882353|
|    4| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 4}      |  0.6470588|    0.7058824|  0.5882353|
|    5| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}      |  0.6470588|    0.7058824|  0.5882353|
|    6| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.6470588|    0.7058824|  0.5882353|
|    7| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}         |  0.6470588|    0.7058824|  0.5882353|
|    8| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 4} |  0.6470588|    0.7058824|  0.5882353|
|    9| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4} |  0.6470588|    0.7058824|  0.5882353|
|   11| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 4}       |  0.6470588|    0.7058824|  0.5882353|
|   12| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 8}       |  0.6470588|    0.7058824|  0.5882353|
|   15| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.6470588|    0.7058824|  0.5882353|
|   17| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.6470588|    0.7058824|  0.5882353|
|   18| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}      |  0.6470588|    0.7058824|  0.5882353|
|   20| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |  0.6470588|    0.7058824|  0.5882353|
|   26| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 4}           |  0.6470588|    0.6470588|  0.6470588|
|   29| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 4}         |  0.6470588|    0.2941176|  1.0000000|
|   30| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}        |  0.6470588|    0.3529412|  0.9411765|
|   31| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}       |  0.6470588|    0.2941176|  1.0000000|
|   33| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 4}     |  0.6470588|    0.6470588|  0.6470588|
|   34| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |  0.6470588|    0.2941176|  1.0000000|
|   44| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 8}         |  0.6470588|    0.3529412|  0.9411765|
|   45| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 8}        |  0.6470588|    0.2941176|  1.0000000|
|   47| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}           |  0.6470588|    0.2941176|  1.0000000|
|   49| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}        |  0.6470588|    0.2941176|  1.0000000|
|    3| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 4}        |  0.6176471|    0.6470588|  0.5882353|
|   10| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}      |  0.6176471|    0.6470588|  0.5882353|
|   14| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 4}    |  0.6176471|    0.7058824|  0.5294118|
|   16| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |  0.6176471|    0.6470588|  0.5882353|
|   24| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}        |  0.6176471|    0.3529412|  0.8823529|
|   25| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 4}    |  0.6176471|    0.6470588|  0.5882353|
|   28| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}             |  0.6176471|    0.2941176|  0.9411765|
|   35| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}         |  0.6176471|    0.2941176|  0.9411765|
|   36| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |  0.6176471|    0.4705882|  0.7647059|
|   38| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}        |  0.6176471|    0.4117647|  0.8235294|
|   40| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 8, ‘epochs’: 30, ‘batch\_size’: 8}           |  0.6176471|    0.2352941|  1.0000000|
|   41| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 8}           |  0.6176471|    0.3529412|  0.8823529|
|   46| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}       |  0.6176471|    0.2941176|  0.9411765|
|   13| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |  0.5882353|    0.4705882|  0.7058824|
|   23| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.5, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}         |  0.5882353|    0.5882353|  0.5882353|
|   37| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}         |  0.5882353|    0.2352941|  0.9411765|
|   39| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}       |  0.5882353|    0.3529412|  0.8235294|
|   48| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 64, ‘epochs’: 30, ‘batch\_size’: 8}            |  0.5882353|    0.4705882|  0.7058824|
|   22| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.5, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}        |  0.5294118|    0.4705882|  0.5882353|

### LSTM-Embedding

``` r
oomycete_lstm_embd_all_results <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/all_scan_results_lstm_emb_scan_oomycete.csv") %>% 
  dplyr::select(c(params, mean_test_score))

oomycete_lstm_embd_all_results %>% 
  arrange(desc(mean_test_score)) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                    |  mean\_test\_score|
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}    |          0.6691176|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.6544118|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.6397059|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.6250000|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}    |          0.6250000|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.6250000|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.6176471|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.6102941|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.6102941|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.6102941|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}      |          0.6029412|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.6029412|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.6029412|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4} |          0.6029412|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.5955882|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}  |          0.5882353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.5882353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.5808824|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.5808824|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.5661765|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}    |          0.5220588|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}      |          0.5073529|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}   |          0.4926471|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}      |          0.4852941|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}  |          0.4779412|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}  |          0.4779412|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}      |          0.4705882|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.4705882|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.4705882|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}       |          0.4705882|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}   |          0.4632353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}       |          0.4632353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.4632353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}       |          0.4632353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.4632353|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}  |          0.4632353|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.4632353|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.4558824|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.4558824|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.4558824|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.4485294|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.4485294|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.4485294|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.4485294|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.4485294|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.4485294|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.4485294|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}   |          0.4411765|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.4411765|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.4411765|

### CNN-LSTM

``` r
oomycete_cnn_lstm_pred <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/df_pred_results_cnn_lstm_scan_oomycete_secreted_data.csv")
oomycete_cnn_lstm_all <- data.table::fread("../../../../../data/secreted_data/training-results/oomycete/results/secreted_hyper_scan/all_scan_results_cnn_lstm_scan_oomycete_secreted_data.csv")

oomycete_cnn_lstm_pred %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                 |   Accuracy|  Sensitivity|  Specifity|
|----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|   11| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.7352941|    0.6470588|  0.8235294|
|   12| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.7352941|    0.8823529|  0.5882353|
|   16| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.7352941|    0.7647059|  0.7058824|
|   20| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}        |  0.7352941|    0.8823529|  0.5882353|
|    4| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.7058824|    0.7058824|  0.7058824|
|    7| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.7058824|    0.8235294|  0.5882353|
|   13| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.6764706|    0.7058824|  0.6470588|
|   21| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.6764706|    0.7058824|  0.6470588|
|   27| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.6764706|    0.8823529|  0.4705882|
|   28| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6764706|    0.7647059|  0.5882353|
|    1| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6470588|    0.6470588|  0.6470588|
|    5| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.6470588|    0.7058824|  0.5882353|
|   25| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.6470588|    0.8235294|  0.4705882|
|    0| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.6176471|    0.5294118|  0.7058824|
|    2| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}        |  0.6176471|    0.6470588|  0.5882353|
|    8| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.6176471|    0.7058824|  0.5294118|
|   17| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.6176471|    0.5882353|  0.6470588|
|   22| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.6176471|    0.8235294|  0.4117647|
|   26| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.6176471|    0.6470588|  0.5882353|
|    9| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.5882353|    0.7058824|  0.4705882|
|   19| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}    |  0.5882353|    0.6470588|  0.5294118|
|   29| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.5882353|    0.4705882|  0.7058824|
|    3| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.5588235|    0.5294118|  0.5882353|
|   10| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.5588235|    0.7058824|  0.4117647|
|   18| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.5294118|    0.6470588|  0.4117647|
|   23| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.5294118|    0.8823529|  0.1764706|
|   14| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.5000000|    0.5294118|  0.4705882|
|   15| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}        |  0.5000000|    0.5294118|  0.4705882|
|   24| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.5000000|    0.8235294|  0.1764706|
|    6| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.4411765|    0.4705882|  0.4117647|

Fungi
-----

### CNN-GRU

``` r
fungi_cnn_gru_secreted <- data.table::fread("../../../../../data/secreted_data/training-results/fungi/df_pred_results_cnn_gru_fungi_secreted.csv")

fungi_cnn_gru_secreted %>% 
  arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                       |   Accuracy|  Sensitivity|  Specifity|
|----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|   22| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}           |  0.8157895|    0.7368421|  0.8947368|
|   21| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}         |  0.7894737|    0.7894737|  0.7894737|
|    2| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}          |  0.7631579|    0.6315789|  0.8947368|
|    5| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}       |  0.7631579|    0.8421053|  0.6842105|
|   11| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7631579|    0.7894737|  0.7368421|
|    8| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}       |  0.7631579|    0.7894737|  0.7368421|
|    6| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}       |  0.7368421|    0.7894737|  0.6842105|
|   10| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.7368421|    0.7368421|  0.7368421|
|   19| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}         |  0.7105263|    0.5789474|  0.8421053|
|   14| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |  0.7105263|    0.6315789|  0.7894737|
|   15| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  0.7105263|    0.6842105|  0.7368421|
|   17| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}   |  0.7105263|    0.7368421|  0.6842105|
|   24| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7105263|    0.6842105|  0.7368421|
|    4| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}       |  0.7105263|    0.6315789|  0.7894737|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}        |  0.6842105|    0.5263158|  0.8421053|
|    7| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}       |  0.6842105|    0.8421053|  0.5263158|
|   13| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}      |  0.6578947|    0.6842105|  0.6315789|
|   25| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.6578947|    0.5789474|  0.7368421|
|   27| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}      |  0.6578947|    0.4736842|  0.8421053|
|   18| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}      |  0.6315790|    0.5789474|  0.6842105|
|    9| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}     |  0.6315789|    0.6842105|  0.5789474|
|   16| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’} |  0.6052632|    0.6315789|  0.5789474|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}          |  0.5789474|    0.5789474|  0.5789474|
|   20| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.5789474|    0.3684211|  0.7894737|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: ‘relu’}         |  0.5789474|    0.3157895|  0.8421053|
|   12| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}            |  0.5789474|    0.4210526|  0.7368421|
|   26| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}         |  0.5526316|    0.4736842|  0.6315789|
|   29| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}            |  0.5263158|    0.4736842|  0.5789474|
|   28| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 4, ‘activation\_conv’: None}    |  0.4736842|    0.3157895|  0.6315789|
|   23| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.4210526|    0.4736842|  0.3684211|

### LSTM Emb

``` r
fungi_lstm_emb_all <- data.table::fread("../../../../../data/secreted_data/training-results/fungi/all_scan_results_lstm_emb_scan_fungi.csv") %>% 
   dplyr::select(c(params, mean_test_score))
 

fungi_lstm_emb_all %>% 
  dplyr::arrange(desc(mean_test_score)) %>% 
  head(10) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                    |  mean\_test\_score|
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.7179487|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.7179487|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.7115385|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}  |          0.7115385|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.6987180|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.6987180|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4}  |          0.6987179|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 4} |          0.6987179|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 4}     |          0.6987179|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.6923077|

Bacteria
--------

CNN-GRU
-------

``` r
bacteria_cnn_gru_secreted <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/df_pred_results_cnn_gru_bacteria_secreted.csv") 


bacteria_cnn_gru_pred_all <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/all_scan_results_cnn_gru_bacteria_secreted.csv") %>% 
  dplyr::select(c(params, mean_test_score))


bacteria_cnn_gru_secreted %>% 
  dplyr::arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                        |   Accuracy|  Sensitivity|  Specifity|
|----:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    6| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}         |  1.0000000|    1.0000000|  1.0000000|
|    7| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  1.0000000|    1.0000000|  1.0000000|
|   12| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}        |  1.0000000|    1.0000000|  1.0000000|
|   19| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}        |  0.9868421|    1.0000000|  0.9736842|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}    |  0.9868421|    0.9736842|  1.0000000|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}     |  0.9868421|    0.9736842|  1.0000000|
|   10| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}         |  0.9868421|    1.0000000|  0.9736842|
|   16| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}          |  0.9868421|    1.0000000|  0.9736842|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  0.9736842|    0.9736842|  0.9736842|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}      |  0.9736842|    1.0000000|  0.9473684|
|    8| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.9736842|    0.9736842|  0.9736842|
|    9| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}      |  0.9736842|    0.9473684|  1.0000000|
|   14| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.9736842|    0.9736842|  0.9736842|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.9736842|    1.0000000|  0.9473684|
|    5| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.9605263|    0.9736842|  0.9473684|
|   13| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}      |  0.9605263|    0.9473684|  0.9736842|
|   21| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}        |  0.9605263|    0.9736842|  0.9473684|
|   11| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}            |  0.9605263|    0.9473684|  0.9736842|
|   15| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.9605263|    0.9473684|  0.9736842|
|   17| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}           |  0.9605263|    0.9473684|  0.9736842|
|   18| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.9605263|    0.9473684|  0.9736842|
|   20| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.9605263|    0.9736842|  0.9473684|
|   22| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}      |  0.9210526|    0.9210526|  0.9210526|
|   25| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.8289474|    0.8684211|  0.7894737|
|   29| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |  0.8289474|    0.8421053|  0.8157895|
|   23| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}            |  0.8026316|    0.8157895|  0.7894737|
|   24| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}      |  0.7894737|    0.7631579|  0.8157895|
|   26| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |  0.7631579|    0.7368421|  0.7894737|
|   27| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}          |  0.7631579|    0.7894737|  0.7368421|
|   28| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}      |  0.7368421|    0.7368421|  0.7368421|

``` r
bacteria_cnn_gru_pred_all %>% 
  dplyr::arrange(desc(mean_test_score)) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                                                                                            |  mean\_test\_score|
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |          0.9407895|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}    |          0.9407895|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |          0.9309211|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}     |          0.9276316|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}      |          0.9276316|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |          0.9243421|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}         |          0.9243421|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |          0.9243421|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |          0.9210526|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}      |          0.9210526|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}         |          0.9177632|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}            |          0.9144737|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}        |          0.9111842|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}      |          0.9078947|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |          0.9078947|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |          0.9078947|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}          |          0.9046053|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}           |          0.9046053|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |          0.8980263|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}        |          0.8980263|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |          0.8947368|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}        |          0.8815789|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}      |          0.8815789|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}            |          0.7828947|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}      |          0.7730263|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |          0.7697368|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}        |          0.7565789|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}          |          0.7532895|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}      |          0.7500000|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |          0.7401316|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}          |          0.7236842|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |          0.7203947|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}          |          0.7072368|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}            |          0.7039474|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}        |          0.6907895|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}         |          0.6447368|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}          |          0.6414474|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |          0.6381579|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}           |          0.6085526|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}       |          0.6085526|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}     |          0.5822368|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}      |          0.5789474|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}     |          0.5789474|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}      |          0.5690789|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}        |          0.5657895|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}      |          0.5526316|
| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}          |          0.5328947|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 32, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}          |          0.5328947|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 1, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}        |          0.5263158|
| {‘reg\_rate’: 0.01, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0.25, ‘opt\_dropout’: 0.25, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}   |          0.5197368|

CNN-LSTM
--------

``` r
bacteria_cnn_lstm_secreted <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/df_pred_results_cnn_lstm_scan_bacteria_secreted.csv") 


bacteria_cnn_lstm_pred_all <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/all_scan_results_cnn_lstm_scan_bacteria_secreted.csv") %>% 
  dplyr::select(c(params, mean_test_score))


bacteria_cnn_lstm_secreted %>% 
  dplyr::arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                |   Accuracy|  Sensitivity|  Specifity|
|----:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|   13| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  1.0000000|    1.0000000|  1.0000000|
|    0| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.9868421|    0.9736842|  1.0000000|
|   11| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.9868421|    0.9736842|  1.0000000|
|    3| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.9736842|    0.9736842|  0.9736842|
|    4| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.9736842|    0.9736842|  0.9736842|
|   22| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.9736842|    1.0000000|  0.9473684|
|    1| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |  0.9605263|    1.0000000|  0.9210526|
|    7| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.9605263|    1.0000000|  0.9210526|
|   15| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.9605263|    0.9736842|  0.9473684|
|   16| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.9605263|    1.0000000|  0.9210526|
|    5| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.9473684|    0.9736842|  0.9210526|
|   20| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.9473684|    0.9473684|  0.9473684|
|    6| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.9473684|    0.9736842|  0.9210526|
|   26| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.9473684|    0.9736842|  0.9210526|
|   27| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.9342105|    0.8947368|  0.9736842|
|    2| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.9342105|    0.9473684|  0.9210526|
|    9| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.9342105|    0.9473684|  0.9210526|
|   12| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.9342105|    0.9736842|  0.8947368|
|    8| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |  0.9210526|    0.8421053|  1.0000000|
|   17| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.9210526|    0.8947368|  0.9473684|
|   25| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.9210526|    0.9736842|  0.8684211|
|   23| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.9210526|    0.9736842|  0.8684211|
|   28| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.9210526|    0.9736842|  0.8684211|
|   10| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.9078947|    0.9210526|  0.8947368|
|   24| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.9078947|    0.9736842|  0.8421053|
|   18| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |  0.9078947|    0.9736842|  0.8421053|
|   14| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.8815789|    0.9473684|  0.8157895|
|   19| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |  0.8157895|    0.8684211|  0.7631579|
|   29| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |  0.8157895|    0.7894737|  0.8421053|
|   21| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |  0.7763158|    0.8421053|  0.7105263|

``` r
bacteria_cnn_lstm_pred_all %>% 
  dplyr::arrange(desc(mean_test_score)) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                                    |  mean\_test\_score|
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.9276316|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |          0.9177632|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |          0.9111842|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |          0.9078947|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |          0.9046053|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 16, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |          0.9046053|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.8980263|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.8980263|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}     |          0.8947368|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |          0.8914474|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.8881579|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 16, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’} |          0.8881579|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |          0.8881579|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.8848684|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.8815789|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.8782895|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |          0.8684211|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |          0.8684211|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |          0.8651316|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 8, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |          0.8585526|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}       |          0.8453947|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.8421053|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.8355263|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.8322368|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.8256579|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |          0.8223684|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.8157895|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adam’, ‘number\_hidden\_units’: 8, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}      |          0.8059211|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}   |          0.7763158|
| {‘strides’: 1, ‘padding’: ‘valid’, ‘optimizers’: ‘Adadelta’, ‘number\_hidden\_units’: 4, ‘filters\_LSTM’: 4, ‘filters’: 4, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_convolution’: None, ‘activation\_LSTM’: ‘tanh’}  |          0.7434211|

GRU-EMb
-------

``` r
bacteria_gru_emb_pred <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/df_pred_results_gru_embedding_scan_all.csv") 


bacteria_gru_emb_all <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/bacteria_scan_results_gru_embedding_scan_all.csv") %>% 
  dplyr::select(c(params, mean_test_score))

bacteria_gru_emb_all %>% 
  dplyr::arrange(desc(mean_test_score)) %>% 
  head(10) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                     |  mean\_test\_score|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.9407895|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16} |          0.9407895|
| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32} |          0.9375000|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.9375000|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16} |          0.9342105|
| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.9309211|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32} |          0.9276316|
| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32} |          0.9276316|
| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.9276316|
| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.9243421|

``` r
bacteria_gru_emb_pred %>% 
  dplyr::arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                 |   Accuracy|  Sensitivity|  Specifity|
|----:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    5| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |  0.9868421|    0.9736842|          1|
|    6| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32} |  0.9868421|    0.9736842|          1|
|    7| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32} |  0.9868421|    0.9736842|          1|
|    8| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |  0.9868421|    0.9736842|          1|
|    9| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.9868421|    0.9736842|          1|
|   10| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}  |  0.9868421|    0.9736842|          1|
|   11| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |  0.9868421|    0.9736842|          1|
|   16| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}       |  0.9868421|    0.9736842|          1|
|    0| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |  0.9736842|    0.9473684|          1|
|    2| {‘reg\_rate’: 0.001, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32} |  0.9736842|    0.9473684|          1|
|    3| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.9736842|    0.9473684|          1|
|    4| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16} |  0.9736842|    0.9473684|          1|
|   12| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}   |  0.9736842|    0.9473684|          1|
|   13| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |  0.9736842|    0.9473684|          1|
|   14| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |  0.9736842|    0.9473684|          1|
|   15| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}   |  0.9736842|    0.9473684|          1|
|   17| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |  0.9736842|    0.9473684|          1|
|   18| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}      |  0.9736842|    0.9473684|          1|
|   19| {‘reg\_rate’: 0.01, ‘outputdim’: 48, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}      |  0.9736842|    0.9473684|          1|
|    1| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adadelta’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘gru\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16} |  0.9605263|    0.9210526|          1|

LSTM-Emb
--------

``` r
bacteria_latm_emb_scan_all <- data.table::fread("../../../../../data/secreted_data/training-results/bacteria/all_scan_results_lstm_emb_scan_bacteria.csv") %>% 
   dplyr::select(c(params, mean_test_score))
 

bacteria_latm_emb_scan_all %>% 
  dplyr::arrange(desc(mean_test_score)) %>% 
  head(10) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                     |  mean\_test\_score|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.9506579|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}       |          0.9473684|
| {‘reg\_rate’: 0.01, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.9473684|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.9440789|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16} |          0.9440789|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}     |          0.9440789|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}    |          0.9407895|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.9407895|
| {‘reg\_rate’: 0.001, ‘outputdim’: 16, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.9407895|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}       |          0.9407895|

All data together
-----------------

``` r
all_data_cnn_gru_pred_res <- data.table::fread("../../../../../data/secreted_data/training-results/all/df_pred_results_cnn_gru_all.csv")

all_data_cnn_gru_pred_res %>% 
  dplyr::arrange(desc(Accuracy)) %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                  |   Accuracy|  Sensitivity|  Specifity|
|----:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|   25| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: ‘relu’}  |  0.8445946|    0.8108108|  0.8783784|
|    0| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.8310811|    0.8108108|  0.8513514|
|   15| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.8243243|    0.7837838|  0.8648649|
|   21| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: ‘relu’}  |  0.8175676|    0.7702703|  0.8648649|
|    8| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: None}     |  0.8108108|    0.7297297|  0.8918919|
|    2| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’} |  0.8040541|    0.7702703|  0.8378378|
|   26| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |  0.8040541|    0.7702703|  0.8378378|
|   18| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}     |  0.7972973|    0.6756757|  0.9189189|
|    9| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |  0.7972973|    0.7297297|  0.8648649|
|   23| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: None}     |  0.7972973|    0.8378378|  0.7567568|
|   22| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}    |  0.7972973|    0.8648649|  0.7297297|
|   19| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}     |  0.7905405|    0.7837838|  0.7972973|
|    7| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7905405|    0.8243243|  0.7567568|
|   12| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7905405|    0.7297297|  0.8513514|
|   14| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}      |  0.7905405|    0.8108108|  0.7702703|
|   17| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: None}     |  0.7837838|    0.7027027|  0.8648649|
|   13| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}    |  0.7837838|    0.7972973|  0.7702703|
|   24| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: None}    |  0.7837838|    0.7702703|  0.7972973|
|    5| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.7770270|    0.7702703|  0.7837838|
|   11| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}    |  0.7770270|    0.7567568|  0.7972973|
|    6| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7770270|    0.7567568|  0.7972973|
|   10| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 32, ‘activation\_conv’: ‘relu’} |  0.7770270|    0.7702703|  0.7837838|
|   16| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 8, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}   |  0.7702703|    0.7432432|  0.7972973|
|    1| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}   |  0.7702703|    0.7432432|  0.7972973|
|   20| {‘reg\_rate’: 0.01, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: None}     |  0.7702703|    0.7027027|  0.8378378|
|    4| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}     |  0.7500000|    0.8378378|  0.6621622|
|   27| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 3, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.7500000|    0.7432432|  0.7567568|
|    3| {‘reg\_rate’: 0.001, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 32, ‘epochs’: 30, ‘batch\_size’: 16, ‘activation\_conv’: ‘relu’}  |  0.7500000|    0.8378378|  0.6621622|
|   29| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 8, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: None}      |  0.6756757|    0.6081081|  0.7432432|
|   28| {‘reg\_rate’: 0.001, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘maxpool\_size’: 2, ‘kernel\_size’: 2, ‘gru\_hidden\_units’: 16, ‘filter\_conv’: 16, ‘epochs’: 30, ‘batch\_size’: 8, ‘activation\_conv’: ‘relu’}   |  0.6081081|    0.5810811|  0.6351351|

``` r
all_data_lstm_emb_all <- data.table::fread("../../../../../data/secreted_data/training-results/all/all_scan_results_lstm_emb_scan_all_data.csv") %>% select(params, mean_test_score)

all_data_lstm_emb_all %>% 
  dplyr::arrange(desc(mean_test_score)) %>% 
  knitr::kable()
```

| params                                                                                                                                                                                                     |  mean\_test\_score|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.8120805|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.8087248|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.8070470|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.8070470|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.8053691|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}   |          0.8020134|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}    |          0.8020134|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}     |          0.8003356|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.8003356|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32} |          0.7986577|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}    |          0.7986577|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7969799|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7969799|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.7953020|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16} |          0.7953020|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.7953020|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.7919463|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.7902685|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7902685|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.7902685|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.7885906|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}      |          0.7869128|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}      |          0.7852349|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.5318792|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.5151007|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}    |          0.5117450|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}   |          0.5117450|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}   |          0.5100671|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}   |          0.5100671|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.5100671|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}    |          0.5016779|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.5016779|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}       |          0.5000000|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.4983221|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}      |          0.4983221|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}   |          0.4983221|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}      |          0.4983221|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.4983221|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.4916107|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}     |          0.4882550|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.4865772|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 32}   |          0.4832215|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}       |          0.4798658|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}     |          0.4798658|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}      |          0.4798658|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}    |          0.4765101|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.4748322|
| {‘reg\_rate’: 0.01, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}    |          0.4731544|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}       |          0.4664430|
| {‘reg\_rate’: 0.001, ‘outputdim’: 32, ‘optimizers’: ‘sgd’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 32}  |          0.4614094|

Summary Comparing results of training data with the effector secreted and non-secreted
--------------------------------------------------------------------------------------

``` r
data.frame(model = c("cnn_gru_bacteria", 
                         "cnn_gru_bacteria", 
                         "cnn_gru_fungi", 
                         "cnn_gru_fungi", 
                         "cnn_gru_oomycete", 
                         "cnn_gru_oomycete"), 
           data = c("non-secreted", 
                    "secreted", 
                    "non-secreted", 
                    "secreted", 
                    "non-secreted", 
                    "secreted"),
           accuracy = c(bacteria_cnn_gru %>% dplyr::select(Accuracy) %>% max(), 
                        bacteria_cnn_gru_secreted %>% dplyr::select(Accuracy) %>% max(),
                        fungi_cnn_gru %>% dplyr::select(Accuracy) %>% max(),
                        fungi_cnn_gru_secreted %>% dplyr::select(Accuracy) %>% max(),
                        oomycete_cnn_gru %>% dplyr::select(Accuracy) %>% max(),
                        oomycete_cnn_gru_secreted %>% dplyr::select(Accuracy) %>% max())
)
```
