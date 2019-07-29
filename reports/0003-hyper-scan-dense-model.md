Hyperparameters scan of fully connected dense network to predict effector protein
=================================================================================

Aim
---

### Question

Can the simplest deep learning model (fully connected model) with
certain combination of hyperparamaters give us a good accuracy to
predict effector and non-effector protein? Which hyperparamater
combination will give the best accuracy?

### Purpose

My purpose is to find the hyperparamater combination that can give us
best accuracy model.

Method
------

### Procedure

In order to achieve what I aim for, I will conduct a hyperparameter scan
using simple fully connected dense network as a base model. I will use
RandomizedSearchCV() from Scikit-learn that will scan the combination of
hyperparameters randomly.

Execution
---------

I sent a job in GPUs, the scripts is available on github.

Results
-------

``` r
result_hyper_scan <- data.table::fread("data/result_hyper_tuned_old.csv", drop = 'V1') %>% 
  dplyr::select(params, mean_train_score, mean_test_score)

result_hyper_scan %>% 
  top_n(25) %>% 
  knitr::kable() 
```

    ## Selecting by mean_test_score

| params                                                                                                                                                                                                                                                                           |  mean\_train\_score|  mean\_test\_score|
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------:|------------------:|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}            |           0.9979675|          0.6975610|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 3\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |           0.9378049|          0.7024390|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9959350|          0.7154472|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}        |           0.9894309|          0.6991870|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |           0.9979675|          0.6975610|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |           0.9857724|          0.7089431|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2, 3\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9983740|          0.7040650|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |           0.9878049|          0.7349593|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |           0.9971545|          0.6975610|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9979675|          0.7024390|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}    |           0.9975610|          0.6991870|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[2, 3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9983740|          0.7089431|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1, 2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}    |           0.9983740|          0.7154472|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘yes’, ‘activation\_function’: ‘relu’}          |           0.9861789|          0.7008130|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9983740|          0.7073171|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}  |           0.9983740|          0.7056911|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |           0.9955285|          0.7170732|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}     |           0.9983740|          0.7154472|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9979675|          0.7138211|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |           0.9979675|          0.7252033|
| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}           |           0.9967480|          0.7219512|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |           0.9979675|          0.7008130|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}        |           0.9983740|          0.7073171|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |           0.9979675|          0.7073171|
| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’} |           0.9975610|          0.7154472|

``` r
predict_results <- data.table::fread("data/df_result_prediction_old.csv")

predict_results %>% 
  knitr::kable()
```

|   V1| Parameters                                                                                                                                                                                                                                                                       |   Accuracy|  Sensitivity|  Specifity|
|----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|------------:|----------:|
|    0| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |  0.6800000|    0.6710526|  0.6891892|
|    1| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.7200000|    0.7105263|  0.7297297|
|    2| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}           |  0.6733333|    0.6973684|  0.6486486|
|    3| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |  0.6933333|    0.7500000|  0.6351351|
|    4| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.7000000|    0.6184211|  0.7837838|
|    5| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1, 2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}    |  0.6666667|    0.5789474|  0.7567568|
|    6| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’} |  0.7066667|    0.6052632|  0.8108108|
|    7| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}     |  0.7333333|    0.7500000|  0.7162162|
|    8| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.7266667|    0.8421053|  0.6081081|
|    9| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |  0.7066667|    0.7368421|  0.6756757|
|   10| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[2, 3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.6933333|    0.5921053|  0.7972973|
|   11| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.7333333|    0.7105263|  0.7567568|
|   12| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |  0.7000000|    0.5657895|  0.8378378|
|   13| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}        |  0.7000000|    0.6973684|  0.7027027|
|   14| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1, 2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}  |  0.7666667|    0.8157895|  0.7162162|
|   15| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2, 3\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 16, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.7400000|    0.8289474|  0.6486486|
|   16| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 3\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}       |  0.7266667|    0.7368421|  0.7162162|
|   17| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 4, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}      |  0.7066667|    0.5526316|  0.8648649|
|   18| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |  0.6866667|    0.6842105|  0.6891892|
|   19| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘yes’, ‘activation\_function’: ‘relu’}          |  0.6400000|    0.6447368|  0.6351351|
|   20| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}        |  0.7266667|    0.8026316|  0.6486486|
|   21| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[2, 2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}    |  0.7333333|    0.7763158|  0.6891892|
|   22| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}            |  0.7400000|    0.7631579|  0.7162162|
|   23| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1, 3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |  0.7733333|    0.8026316|  0.7432432|
|   24| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |  0.6400000|    0.4473684|  0.8378378|
|   25| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[3\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 2, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 20, ‘dropout\_rates’: 0.25, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |  0.7533333|    0.8684211|  0.6351351|
|   26| {‘shuffle’: True, ‘optim\_methods’: ‘Adadelta’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 8, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}     |  0.7400000|    0.7894737|  0.6891892|
|   27| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[2\], ‘l2\_rate’: 0.01, ‘input\_num\_hidden\_units’: 1, ‘input\_dropout\_rates’: 0.5, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 32, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}           |  0.6533333|    0.6315789|  0.6756757|
|   28| {‘shuffle’: True, ‘optim\_methods’: ‘SGD’, ‘num\_hidden\_layers’: \[1\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 16, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 30, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘no’, ‘activation\_function’: ‘relu’}          |  0.7400000|    0.8026316|  0.6756757|
|   29| {‘shuffle’: True, ‘optim\_methods’: ‘Adam’, ‘num\_hidden\_layers’: \[0\], ‘l2\_rate’: 0.001, ‘input\_num\_hidden\_units’: 3, ‘input\_dropout\_rates’: 0.25, ‘epochs’: 20, ‘dropout\_rates’: 0.5, ‘batch\_size’: 8, ‘batch\_norm’: ‘yes’, ‘activation\_function’: ‘relu’}         |  0.6933333|    0.7500000|  0.6351351|
