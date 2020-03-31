
# October 2019 Lab Report <img src="figures/tsl-logo.png" align="right" width="120" />

### [Ruth Kristianingsih](https://github.com/ruthkr)

## 01 October 2019 (Tuesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Read about having K-fold for manually training.
  - Made a script using K-fold to train models manually (running on
    cluster).
  - Continued the reports.

## 02–04 October 2019– Sick Days

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

-----

## 07–08 October 2019– Sick Days

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

## 09 October 2019 (Wednesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Made a test script for subcellular multicompartment to see how long
    approximately one model finish to run:
    `0001-getting-started-cnn-lstm.py`.
  - Made a script using the base CNN-LSTM with CVs with all parameters
    that has a very good perfomance in nuclear/non-nuclear subcellular
    localisation (failed since the kfold() from Scikit learn did not
    support multiclass classification), manual scripts has not worked
    yet either:
      - `best_CV5.py`.
      - `test.py`.

## 10 October 2019 (Thursday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Made scripts for Grid-search for LSTM-Embedding.
  - Continued project-report.
  - Made progress-report.

## 11 October 2019 (Friday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Continued project-report.
  - Made progress-report.
  - Created several scripts for subcellular multicompartment prediction,
    involving:
      - Use the base model: CNN-LSTM.
      - The scripts are:
          - `filter48_LSTM48.py`.
          - `filter48_LSTM64.py`.
          - `filter64_LSTM64.py`.
          - `filter64_LSTM48.py`.

-----

## 14 October 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Weekly catch up minutes

##### Sub-cellular protein project

  - Progress:
      - From the RandomSearch hyperparameter scan with total (16 models
        with CV=5, so in total is 80 models, only 16 models have been
        trained ==\> 20%). Hyperparamater scan will be finished in
        around two weeks.
  - Points:
      - While waiting for hyperparameter scan finished, then take one
        model and tune manually.

##### Effector and non-effector prediction

  - Progress:
      - CNN-LSTM: got 6 models.
      - CNN-GRU: got 3 models.
      - LSTM-Embedding: still in progress of Grid search.
      - GRU-Embedding: still in progress of RandomSearch.
  - Points:
      - Dan suggested to train manually for the CNN-LSTM and GRU and see
        if it;s overfitting or not. Also for LSTM-Embdedding and
        GRU-Embdedding if the hyperparameters scan finish. To reduce
        overfitting, can try some regularization methods, such as:
        dropout, pre-norm step, optimizer:sgd, can change the batch
        size.
      - Check also overfitting section on the Deep Learning Books.
      - Can also try to shrink the network and Earlystopping.
      - New challenge: Do the ensemble models using weighted average
        method, where the weight would be the accuracy of each model.

## 15–18 October 2019

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Wrote scripts for model ensemble (effector and non-effector
    projects).
  - Contined the report for both projects.
  - Worked on the regularization method for the model.
  - Read some papers about ensembling models.

-----

## 21 October 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Weekly catch up minutes

##### Sub-cellular protein project

  - Points:
      - Based on parallel jobs running in GPU of protein localisation on
        multicompartment projects (log files), the best accuracy on the
        training data so far is already 99% (potentially overfitting).
        Then Dan suggested to stop it and run it manually (tuning), and
        see if the model is overfitting or not.
  - Plans:
      - Stop all jobs and tune one by one the best models so far.

##### Effector and non-effector prediction

  - Points:
      - Try another regularization method.
      - Making ensemble models ==\> weighted average ensemble, and use
        accuracy for each model for the weights.
  - Plans:
      - For each models, pick the best accuracy so far, and tune it to
        get least overfitting models.
      - Construct the ensemble model, if the best model has not found
        yet, then use the current models.

## 22 October 2019 (Tuesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Stopped all of the hyperparameters scan of Subcellular compartments
    jobs: 13 models so far with accuracy in training data
  - Continued the scripts of ensemble, previously using the
    layers.Average() layer in keras to ensemble, but since it needs to
    be wieghted therefore, I can not use layers.Average() layer in keras
  - Manually doing regularization tuning for the efector-non effector
    models (run in GPUs):
      - CNN-LSTM:
          - `best1_reg0001_withoutCV_dropout.py`.
          - `best1_reg001_withoutCV.py`.
          - `best1_reg0001_withoutCV.py`.
          - `best1_with_reg.py`.
          - `best1_without_actv_reg_sgd.py`.

## 23 October 2019 (Wednesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Tried to make a script of subcellular multicompartment with
    multi-GPU job (Sam Gallop gave me the scripts).
  - Continued regularization tuning for the efector-non effector models
    (run in GPUs):
      - CNN-GRU:
          - `best1_with_reg.py`.
          - `best1_reg001_withoutCV.py`.
          - `best1_without_actv_bias.py`.

## 24 October 2019 (Thursday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Wrote scripts of manually tuning jobs for the hyperparameters scan
    results of CNN-LSTM models for subcellular multicompartment (without
    any regularization, and with regularization (kernel + bias), also
    using sgd), the file name are:
      - `sub_cnn_lstm_best1.py`.
      - `sub_cnn_lstm_best1_with_reg.py`.
      - `sub_cnn_lstm_best1_with_reg_with_sgd.py`.
  - Continued then regularization tuning for LSTM and GRU Embedding
    (effector and non-effector models):
      - LSTM-Embedding:
          - `best1_with_reg.py`.
      - GRU-Embedding:
          - `best1_with_reg.py`.

## 25 October 2019 (Friday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Finishing the scripts of ensemble models.
  - Run the subcellular multicompartment models again to GPUs.

-----

## 28 October 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Subcellular multicompartment jobs that were sent on 25rd of October
    failed due to the confusion matrix that could not be used to
    evaluate the prediction results of multiclass classification.
  - Fixed the function for ensemble models so that it has binary
    classification (1 or 0) instead of only probabilities.

#### Weekly catch up minutes

##### Sub-cellular protein project

  - Points:
      - Dan suggested if one model does not take so much to run, then
        may not need to run in parallel as what Sam did (not to waste
        time to figure out how to use the parallelism).
      - Also try to reduce the size of filter for both LSTM and CNN. And
        also for the unit of dense layers.

##### Effector and non-effector prediction

  - Points:
      - Dan asked, why did not I just test the ensemble data on the tes
        nting data, instead of test it one by one.
      - Instead of probability, make it 0 and 1 (binary results).

## 29 October 2019 (Tuesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Created scripts for subcellular multicompartment with filter 64,
    filter LTSM 64, hidden unit in dense layer 64 \[Model 4\]:
      - `0001_model4.py` original model without any regularization.
      - `0002_model4_reg_dropout.py` with l2 regularization and dropout.
      - `0003_model4_batchsize.py` with reduced batch size.
      - `filter16_LSTM16_unit16.py` with all reduced filters.
  - Tried to visualize the ensemble results together with all prediction
    of each model.
  - Updated the October lab book.

## 30 October 2019 (Wednesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Fixed the function of the ensemble scripts (generalized the pattern
    to read the `.hdf5` files).
  - In `.Rmd` made function to calculate the accuracy of each model
    using `caret` library and also a function to create venn diagram to
    represent each the prediction results of each models.
  - Created a script with higher number of filter (CNN, LSTM) and also
    hidden units `filter32_LSTM32_unit32.py` since a script for
    subcellular localisation model with `filter16_LSTM16_unit16.py` did
    not give a good accuracy (`acc = 0.1` for the test data).
  - Created scripts for subcellular multicompartment with filter 48,
    filter LTSM 48, hidden unit in dense layer 48 \[Model 1\]:
      - `0001_model1.py` original model without any regularization.
      - `0002_model1_reg_dropout.py` with l2 regularization and dropout.
      - `0003_model1_batchsize.py` with reduced batch size.

## 31 October 2019 (Thursday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Fixed the mistake on subcellular multi-compartment on evaluating the
    test model.
  - Created the diagram Venn for visializing the prediction results.
  - After showing Dan the results, found that for ensemble model, I
    should not change the prediction to binary result (0 and 1), instead
    putting them first to the ensemble formula (the ensemble improved
    the accuracy).
  - Made the confusion matrix.
