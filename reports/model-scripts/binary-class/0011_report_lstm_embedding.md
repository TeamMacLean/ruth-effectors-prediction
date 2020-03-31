Report on GridSearchCV the model LSTM-Embedding
===============================================

Load library
------------

``` r
# Read the all of results 

result1 <- data.table::fread("../../../../results/results/lstm-embedding/all_scan_results_0013-lstm-embedding-grid-batch8_lstm_unit16.csv")
result2 <- data.table::fread("../../../../results/results/lstm-embedding/all_scan_results_0014-lstm-embedding-grid-batch16_lstm_unit16.csv")
result3 <- data.table::fread("../../../../results/results/lstm-embedding/all_scan_results_0015-lstm-embedding-grid-batch8_lstm_unit32.csv")
result4 <- data.table::fread("../../../../results/results/lstm-embedding/all_scan_results_0016-lstm-embedding-grid-batch16_lstm_unit32.csv")
```

``` r
result1 %>% 
  rbind(., result2) %>% 
  rbind(., result3) %>% 
  rbind(., result4) %>% 
  select("params", "mean_test_score") %>% 
  arrange(desc(mean_test_score)) %>% 
  knitr::kable() 
```

| params                                                                                                                                                                                                     |  mean\_test\_score|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------:|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.7495935|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7479675|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.7463415|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.7398374|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7398374|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.7398374|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.7398374|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}   |          0.7365854|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}    |          0.7349594|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.7349593|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.7349593|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.7349593|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.7333333|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}   |          0.7333333|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}  |          0.7333333|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}    |          0.7317073|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.7300813|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}       |          0.7284553|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}   |          0.7284553|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7284553|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.7268293|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16}      |          0.7252033|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}       |          0.7203252|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}  |          0.7186992|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 16} |          0.7186992|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.7170732|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}    |          0.7138211|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 16, ‘epochs’: 30, ‘batch\_size’: 8}     |          0.7138211|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0.25, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16} |          0.7121951|
| {‘reg\_rate’: 0.01, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘FALSE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.7121951|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 8}      |          0.6991870|
| {‘reg\_rate’: 0.001, ‘outputdim’: 64, ‘optimizers’: ‘Adam’, ‘opt\_go\_backwards’: ‘TRUE’, ‘opt\_dropout\_recurrent’: 0, ‘opt\_dropout’: 0, ‘lstm\_hidden\_units’: 32, ‘epochs’: 30, ‘batch\_size’: 16}     |          0.6878049|

``` r
# Read the result of the prediction

pred_1 <- data.table::fread("../../../../results/results/lstm-embedding/df_pred_results_0013-lstm-embedding-grid-batch8_lstm_unit16.csv")
pred_2 <- data.table::fread("../../../../results/results/lstm-embedding/df_pred_results_0014-lstm-embedding-grid-batch16_lstm_unit16.csv")
```

``` r
pred_1 %>% 
  rbind(., pred_2) %>% 
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
