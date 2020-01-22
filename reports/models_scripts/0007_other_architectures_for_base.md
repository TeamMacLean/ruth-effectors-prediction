Getting the architecture for effector and non-effector prediction
=================================================================

Aim
---

### Background

Previously, hyperparamaters scan on the fully connected dense layer
models has been done, the highest accuracy can be achieved was around
77.33 and it could not go higher than that.

### Question

What is the deep learning model which can give better accuracy than the
accuracy obtained from fully connected dense layer?

### Purpose

The purpose of this experiment is to get model architecture that will
give better accuracy than the best accuracy that we obtained from fully
connected dense layer.

Method
------

I will construct some different deep learning model architectures using
random hyperparamaters choice and investigate how the accuracy of the
model will behave. If the model can give accuracy better than 60%, then
we can use the related model as the base models.

### Procedure

#### GRU model

``` python
def simple_lstm(inputdim = 23,
                  outputdim = 32,
                gru_hidden_units = 16,
                ):
      # Create the model
    emb_vecor_length = 32
    model = Sequential()
    model.add(Embedding(input_dim = inputdim,
                        output_dim = outputdim,
                        input_length = 4034))
    model.add(GRU(gru_hidden_units))
    model.add(Bidirectional(GRU(gru_hidden_units * 2)))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    print(model.summary())
    return model
```

##### Results

``` bash
Parameters:
    outputdim: 32
    inputdim: 23
    gru_hidden_units: 8
    epochs: 30
    batch_size: 16

 16/150 [==>...........................] - ETA: 24s
 32/150 [=====>........................] - ETA: 16s
 48/150 [========>.....................] - ETA: 12s
 64/150 [===========>..................] - ETA: 10s
 80/150 [===============>..............] - ETA: 7s 
 96/150 [==================>...........] - ETA: 5s
112/150 [=====================>........] - ETA: 4s
128/150 [========================>.....] - ETA: 2s
144/150 [===========================>..] - ETA: 0s
150/150 [==============================] - 17s 111ms/step
acc y_pred: 0.68
```

#### LSTM model

``` python
def simple_lstm(inputdim = 23,
                outputdim = 32,
                ):
    # create the model
    emb_vecor_length = 32
    model = Sequential()
    model.add(Embedding(input_dim = inputdim,
                        output_dim = outputdim,
                        input_length = 4034))
    model.add(Bidirectional(LSTM(20)))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    print(model.summary())

    return model
```

##### Results

``` bash
Parameters:
    outputdim: 32
    inputdim: 23
    epochs: 30
    batch_size: 16

 16/150 [==>...........................] - ETA: 13s
 32/150 [=====>........................] - ETA: 9s 
 48/150 [========>.....................] - ETA: 7s
 64/150 [===========>..................] - ETA: 5s
 80/150 [===============>..............] - ETA: 4s
 96/150 [==================>...........] - ETA: 3s
112/150 [=====================>........] - ETA: 2s
128/150 [========================>.....] - ETA: 1s
144/150 [===========================>..] - ETA: 0s
150/150 [==============================] - 10s 65ms/step
acc y_pred: 0.6933333333333334
```

#### GRU – CNN model

``` python
# Define the model architecture
def simple_CNN_GRU(filter_conv = 32,
                   kernel_size = 1,
                   maxpool_size = 2,
                   filter_gru = 32,
                   dropout_gru = 0.1,
                   dropout_reccurent = 0.5,
                   optimizers = 'sgd'
                    ):
    model = Sequential()
    model.add(Conv1D(filter_conv,
                     kernel_size,
                     activation = "relu",
                     input_shape = (4034, 20)))
    model.add(MaxPooling1D(maxpool_size))
    model.add(Bidirectional(GRU(filter_gru,
                                dropout = dropout_gru,
                                recurrent_dropout = dropout_reccurent)))
    model.add(Dense(1, activation = 'sigmoid'))

    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers,
                  metrics = ['accuracy'])

    print(model.summary())

    return model

# Pass the model design to KerasClassifier() wrapper

model = KerasClassifier(build_fn = simple_CNN_GRU, verbose = 1)

# Define the parameters that will be tuned randomly
keras_param_options = {'filter_conv' : [4],
                       'kernel_size' : [1],
                       'maxpool_size' : [2],
                       'filter_gru' : [8],
                       'dropout_gru' : [0.25],
                       'dropout_reccurent' : [0.5],
                       'optimizers' : ['sgd', 'Adam', 'Adadelta'],
                       'batch_size' : [8, 16, 32],
                       'epochs' : [30]}

random_search = RandomizedSearchCV(model,
                                   param_distributions = keras_param_options,
                                   n_iter = 1,
                                   cv = 5,
                                   verbose = 10)
```

##### Results

``` bash
[{'optimizers': 'Adadelta', 'maxpool_size': 2, 'kernel_size': 1, 'filter_gru': 8, 'filter_conv': 4, 'epochs': 30, 'dropout_reccurent': 0.5, 'dropout_gru': 0.25, 'batch_size': 8}]

  8/150 [>.............................] - ETA: 42s
 16/150 [==>...........................] - ETA: 25s
 24/150 [===>..........................] - ETA: 19s
 32/150 [=====>........................] - ETA: 16s
 40/150 [=======>......................] - ETA: 14s
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
150/150 [==============================] - 15s 98ms/step
acc y_pred: 0.6266666666666667
```

#### LSTM – CNN model

``` python
def conv1d_bn(x,
              filters,
              kernel_size,
              padding = 'same',
              strides = 1,
              activation_convolution = 'relu'):
    """Utility function to apply conv + BN.
    # Arguments
        x: input tensor.
        filters: filters in `Conv1D`.
        num_row: height of the convolution kernel.
        num_col: width of the convolution kernel.
        padding: padding mode in `Conv1D`.
        strides: strides in `Conv1D`.

    # Returns
        Output tensor after applying `Conv2D` and `BatchNormalization`.
    """

    x = layers.Conv1D(filters,
                      kernel_size,
                      strides = strides,
                      padding = padding,
                      use_bias = False)(x)
    x = layers.BatchNormalization()(x)
    x = layers.Activation(activation_convolution)(x)
    return x

# Define the base model

def build_model_conv1D_lstm(filters = 48,
                            filters_LSTM = 48,
                            strides = 1,
                            padding = "valid",
                            activation_convolution = 'None',
                            activation_LSTM = 'tanh',
                            optimizers = 'sgd',
                            bn = 'yes',
                            number_hidden_units = 64):

    input = Input(shape = x_train.shape[1:])

    if bn == 'yes':
        convolution_1 = conv1d_bn(x = input,
                                  filters = filters,
                                  kernel_size = 1,
                                  strides = strides,
                                  padding = padding,
                                  activation_convolution = activation_convolution)
    else:
        convolution_1 = Conv1D(filters,
                       kernel_size = 1,
                       strides = 1,
                       padding = padding,
                       activation = activation_convolution,
                       use_bias = True,
                       kernel_initializer = 'glorot_uniform',
                       bias_initializer = 'zeros')(input)
    if bn == "yes":
        convolution_2 = conv1d_bn(x = input,
                                  filters = filters,
                                  kernel_size = 3,
                                  strides = strides,
                                  padding = padding,
                                  activation_convolution = activation_convolution)
    else:
        convolution_2 = Conv1D(filters,
                       kernel_size = 3,
                       strides = strides,
                       padding = padding,
                       activation = activation_convolution,
                       use_bias = True,
                       kernel_initializer = 'glorot_uniform',
                       bias_initializer = 'zeros')(input)

    if bn == "yes":
        convolution_3 = conv1d_bn(x = input,
                                  filters = filters,
                                  kernel_size = 5,
                                  strides = strides,
                                  padding = padding,
                                  activation_convolution = activation_convolution)
    else:
        convolution_3 = Conv1D(filters,
                       kernel_size = 5,
                       strides = strides,
                       padding = padding,
                       activation = activation_convolution,
                       use_bias = True,
                       kernel_initializer = 'glorot_uniform',
                       bias_initializer = 'zeros')(input)

    model = keras.layers.concatenate([convolution_1,
                                      convolution_2,
                                      convolution_3], axis=1)

    model = Conv1D(filters * 2,
                   kernel_size = 3,
                   strides = 1,
                   padding = padding,
                   activation = activation_convolution,
                   use_bias = True,
                   kernel_initializer = 'glorot_uniform',
                   bias_initializer = 'zeros')(model)

    model_lstm_1 = LSTM(filters_LSTM,
                        activation = activation_LSTM,
                        recurrent_activation = 'hard_sigmoid',
                        use_bias = True,
                        kernel_initializer = 'glorot_uniform',
                        recurrent_initializer = 'orthogonal',
                        bias_initializer = 'zeros',
                        go_backwards = False)(model)


    model_lstm_2 = LSTM(filters_LSTM,
                        activation = activation_LSTM,
                        recurrent_activation = 'hard_sigmoid',
                        use_bias = True,
                        kernel_initializer = 'glorot_uniform',
                        recurrent_initializer = 'orthogonal',
                        bias_initializer = 'zeros',
                        go_backwards = True)(model)

    model_final = keras.layers.concatenate([model_lstm_1,
                                            model_lstm_2])

    output = Dense(number_hidden_units, activation = 'relu')(model_final)
    output = Dense(1, activation = 'sigmoid')(model_final)

    model = Model(inputs = input, outputs = output)

    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers,
                  metrics = ['accuracy'])

    print(model.summary())
    return model
```

##### Results

``` bash
Parameters:
    strides: 1
    padding: valid
    optimizers: Adam
    number_hidden_units: 32
    filters_LSTM: 16
    filters: 16
    epochs: 30
    bn: yes
    batch_size: 16
    activation_convolution: None
    activation_LSTM: tanh

 16/150 [==>...........................] - ETA: 32s
 32/150 [=====>........................] - ETA: 24s
 48/150 [========>.....................] - ETA: 19s
 64/150 [===========>..................] - ETA: 16s
 80/150 [===============>..............] - ETA: 12s
 96/150 [==================>...........] - ETA: 9s 
112/150 [=====================>........] - ETA: 6s
128/150 [========================>.....] - ETA: 3s
144/150 [===========================>..] - ETA: 1s
150/150 [==============================] - 28s 184ms/step
acc: 0.6933333333333334
```

#### Results from Random Testing

In summary, here is the table of the results

| Model    | Accuracy |
|----------|----------|
| GRU      | 0.68     |
| LSTM     | 0.693    |
| CNN-GRU  | 0.626    |
| CNN-LSTM | 0.693    |

### Current Results from hyperparamater Scan

``` r
# Load library

library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
# Functions ------------------------------------------------

read_log_into_df_with_params_list <- function(file, params_list, numeric_params) {
  # Filter needed info from raw log, store in a vector of strings
  # lines <- system(paste("grep -E 'loss:.*acc:|Epoch'", file), intern = TRUE)
  # lines <- system(paste("grep -E 'loss:.*acc:|Epoch|\\[CV\\]'", file, "| grep -v 'total'"), intern = TRUE)
  lines <- system(paste("cat ", file, " | tr -d '\\000' | grep -E 'loss:.*acc:|Epoch|\\[CV\\]' | grep -v 'total'"), intern = TRUE)

  # Calculate size of vector
  num_lines <- length(lines)
  num_headers <- grep("Epoch", lines) %>% length()
  num_clean_lines <- num_lines - num_headers

  # We initialize the output vector of strings
  clean_lines <- rep(NA, num_clean_lines)
  count <- 1
  run_count <- 0

  # Loop for processing the lines
  params_str <- NULL
  for (i in 1:num_lines) {
    line <- lines[[i]]
    if (stringr::str_detect(line, "\\[CV\\]")) {
      params_str <- line
    } else if (stringr::str_detect(line, "Epoch")) {
      # Store epoch "header"
      epoch_str <- line
      run_count <- run_count + 1
    } else {
      # Store data/log line
      raw_str <- line

      # Paste and save processed lines
      clean_lines[[count]] <- paste(run_count, "-", params_str, "-", epoch_str, "-", raw_str)
      count <- count + 1
    }
  }

  # Transform vector of strings into data frame
  df <- data.frame(as.list(clean_lines)) %>%
    t() %>%
    as_tibble() %>%
    tibble::remove_rownames()

  # Separate single column into desired columns
  df <- df %>%
    tidyr::separate(V1, c("run", "params", "epoch", "step", "eta", "loss", "accuracy"), sep = "-") %>%
    tidyr::separate(params, params_list, sep = ", ") %>%
    dplyr::mutate_at(vars(params_list), function(x) stringr::str_split(x, "=", simplify = TRUE)[, 2]) %>%
    dplyr::mutate_at(
      numeric_params,
      as.numeric
    )

  return(df)
}
```

``` r
clean_log_df_with_params <- function(data) {

  # Use regex for getting the relevant content of each raw column
  data <- data %>%
    dplyr::mutate(
      epoch = stringr::str_extract(epoch, "[0-9]*/[0-9]*"),
      step = stringr::str_extract(step, "[0-9]*/[0-9]*"),
      loss = stringr::str_extract(loss, "[0-9]*\\.[0-9]*"),
      accuracy = stringr::str_extract(accuracy, "[0-9]*\\.[0-9]*")
    )

  # Change data types and remove useless column
  data <- data %>%
    mutate(
      run = as.numeric(run),
      loss = as.numeric(loss),
      accuracy = as.numeric(accuracy)
    ) %>%
    select(-eta)

  return(data)
}


summarise_log_data_with_params_list <- function(data, params_list) {
  data <- data %>%
    # Get last step of each single run
    group_by_at(vars(c("run", "epoch", params_list))) %>%
    slice(n()) %>%
    # Divide epoch into current and max epoch
    mutate(
      curr_epoch = stringr::str_split(epoch, "/") %>% unlist %>% .[1] %>% as.numeric(),
      max_epoch = stringr::str_split(epoch, "/") %>% unlist %>% .[2] %>% as.numeric(),
    ) %>%
    ungroup() %>%
    # Get final loss/accuracy of each epoch
    dplyr::filter(curr_epoch == max_epoch) %>%
    select(-c(step, epoch, run, curr_epoch)) %>%
    dplyr::rename(epochs = max_epoch) %>%
    mutate(epochs = as.factor(epochs)) %>%
    # Create model variable (5 runs)
    tibble::rowid_to_column(var = "run") %>%
    mutate(
      model = cut(run, breaks = seq(0,1000,5), label = 1:200)
    ) %>%
    select(-run) %>%
    # Summarise results
    group_by_at(vars(c("model", "epochs", params_list))) %>%
    summarise(
      loss_mean = mean(loss),
      loss_sd = sd(loss),
      acc_mean = mean(accuracy),
      acc_sd = sd(accuracy)
    )

  return(data)
}
```

#### CNN\_LSTM\_getting\_started

``` r
src_file <- "../../../results/logs/0004-getting-started-cnn-lstm.log"

src_params_list <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw", "bn",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)
```

``` r
data <- src_file %>%
  read_log_into_df_with_params_list(src_params_list, src_numeric_params) %>%
  clean_log_df_with_params()
```

    ## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
    ## This warning is displayed once per session.

    ## Warning: Expected 7 pieces. Additional pieces discarded in 137 rows [714,
    ## 715, 807, 808, 809, 810, 840, 841, 842, 843, 844, 845, 869, 870, 871, 872,
    ## 1520, 1521, 1522, 1523, ...].

``` r
new_data <- data %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list)
```

``` r
new_data %>% 
  arrange(desc(acc_mean)) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw| bn  |  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|   loss\_sd|  acc\_mean|    acc\_sd|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|:----|------------:|:------------------------|:-----------------|-----------:|----------:|----------:|----------:|
| 2     | 30     |        1| valid   | Adam       |                     48|             48|       48|           30| yes |           48| None                    | tanh             |      0.0084|  0.0022417|    0.99676|  0.0018393|
| 1     | 30     |        1| valid   | Adam       |                     48|             48|       48|           30| yes |           16| None                    | tanh             |      0.0091|  0.0029462|    0.99596|  0.0028815|

``` bash
Parameters:
    strides: 1
    padding: valid
    optimizers: Adam
    number_hidden_units: 64
    filters_LSTM: 48
    filters: 48
    epochs: 30
    bn: yes
    batch_size: 16
    activation_convolution: None
    activation_LSTM: tanh

 16/150 [==>...........................] - ETA: 51s
 32/150 [=====>........................] - ETA: 40s
 48/150 [========>.....................] - ETA: 34s
 64/150 [===========>..................] - ETA: 28s
 80/150 [===============>..............] - ETA: 22s
 96/150 [==================>...........] - ETA: 17s
112/150 [=====================>........] - ETA: 12s
128/150 [========================>.....] - ETA: 7s 
144/150 [===========================>..] - ETA: 1s
150/150 [==============================] - 49s 326ms/step
acc y_pred: 0.7
```

#### GRU

``` r
src_file_gru <- "../../../results/logs/0007-gru-embedding_scan.log"
src_params_list_gru <- c(
  "outputdim", "optimizers", "opt_recurrent_regs", "opt_kernel_regs",
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
  knitr::kable()
```

| model | epochs |  outputdim| optimizers | opt\_recurrent\_regs                                     | opt\_kernel\_regs                                        | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  gru\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  acc\_mean|
|:------|:-------|----------:|:-----------|:---------------------------------------------------------|:---------------------------------------------------------|:-------------------|------------------------:|-------------:|-------------------:|------------:|------------:|-----------:|----------:|
| 3     | 30     |         48| Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7fab7379ffd0&gt; | TRUE               |                     0.50|          0.50|                  64|           30|           16|     0.54952|    0.74266|
| 5     | 30     |         32| Adam       | &lt;keras.regularizers.L1L2 object at 0x7fab7379ffd0&gt; | &lt;keras.regularizers.L1L2 object at 0x7fab851888d0&gt; | FALSE              |                     0.00|          0.00|                   8|           30|           16|     0.55930|    0.73334|
| 6     | 30     |         32| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7fab7379ffd0&gt; | None                                                     | FALSE              |                     0.25|          0.25|                  32|           30|           16|     0.53730|    0.72760|
| 2     | 30     |         64| sgd        | &lt;keras.regularizers.L1L2 object at 0x7fab851888d0&gt; | None                                                     | TRUE               |                     0.50|          0.00|                   8|           30|           16|     0.79744|    0.58498|
| 1     | 30     |         48| sgd        | &lt;keras.regularizers.L1L2 object at 0x7fab851888d0&gt; | None                                                     | FALSE              |                     0.00|          0.50|                   8|           30|           16|     0.79842|    0.57032|
| 4     | 30     |         64| sgd        | &lt;keras.regularizers.L1L2 object at 0x7fab851888d0&gt; | None                                                     | FALSE              |                     0.50|          0.50|                  32|           30|           32|     1.22024|    0.54514|

#### LSTM

``` r
src_file_lstm <- "../../../results/logs/0008-lstm-embedding-scan.log"
src_params_list_lstm <- c(
  "outputdim", "optimizers", "opt_recurrent_regs", "opt_kernel_regs",
  "opt_go_backwards", "opt_dropout_recurrent", "opt_dropout",
  "lstm_hidden_units", "epochs_raw", "batch_size"
)
src_numeric_params_lstm <- c(
  "outputdim", "opt_dropout_recurrent", "opt_dropout",
  "lstm_hidden_units", "epochs_raw", "batch_size"
)
```

``` r
data_lstm <- src_file_lstm %>%
  read_log_into_df_with_params_list(src_params_list_lstm, src_numeric_params_lstm) %>%
  clean_log_df_with_params()

new_data_lstm <- data_lstm %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_lstm)
```

``` r
new_data_lstm %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  knitr::kable()
```

| model | epochs |  outputdim| optimizers | opt\_recurrent\_regs                                     | opt\_kernel\_regs                                        | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  lstm\_hidden\_units|  epochs\_raw|  batch\_size|  loss\_mean|  acc\_mean|
|:------|:-------|----------:|:-----------|:---------------------------------------------------------|:---------------------------------------------------------|:-------------------|------------------------:|-------------:|--------------------:|------------:|------------:|-----------:|----------:|
| 13    | 30     |         64| Adam       | None                                                     | None                                                     | TRUE               |                     0.00|          0.25|                   32|           30|            8|     0.33964|    0.85368|
| 9     | 30     |         64| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | None                                                     | FALSE              |                     0.00|          0.25|                   16|           30|            8|     0.39136|    0.84228|
| 18    | 30     |         64| Adam       | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.00|          0.00|                   32|           30|           16|     0.40770|    0.83130|
| 27    | 30     |         48| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | None                                                     | TRUE               |                     0.00|          0.25|                   16|           30|            8|     0.43632|    0.82518|
| 42    | 30     |         64| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | None                                                     | TRUE               |                     0.25|          0.25|                   48|           30|            8|     0.47336|    0.79878|
| 32    | 30     |         64| Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.00|          0.50|                   16|           30|           32|     0.50388|    0.79308|
| 2     | 30     |         32| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.00|          0.25|                   16|           30|           32|     0.50410|    0.79026|
| 31    | 30     |         32| Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.00|          0.50|                   16|           30|           16|     0.47670|    0.79026|
| 7     | 30     |         48| Adadelta   | None                                                     | None                                                     | TRUE               |                     0.25|          0.25|                   32|           30|           16|     0.46884|    0.78660|
| 5     | 30     |         64| Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.00|          0.50|                   32|           30|           32|     0.51540|    0.78416|
| 41    | 30     |         32| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.00|          0.00|                   16|           30|           16|     0.49694|    0.78292|
| 40    | 30     |         64| Adam       | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.25|          0.25|                   16|           30|           16|     0.50940|    0.77724|
| 8     | 30     |         64| Adam       | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.25|          0.00|                   16|           30|           16|     0.50564|    0.77560|
| 6     | 30     |         48| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.50|          0.00|                   16|           30|            8|     0.52500|    0.76586|
| 23    | 30     |         64| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | None                                                     | TRUE               |                     0.50|          0.00|                   32|           30|           32|     0.50572|    0.76382|
| 38    | 30     |         64| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.50|          0.00|                   32|           30|           16|     0.53538|    0.76340|
| 28    | 30     |         48| Adam       | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.50|          0.25|                   16|           30|            8|     0.52782|    0.76302|
| 12    | 30     |         64| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | None                                                     | TRUE               |                     0.50|          0.50|                   32|           30|           16|     0.52420|    0.76220|
| 21    | 30     |         64| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.00|          0.00|                   16|           30|           32|     0.53910|    0.75976|
| 25    | 30     |         32| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.25|          0.50|                   16|           30|            8|     0.53250|    0.75732|
| 37    | 30     |         48| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.00|          0.00|                   16|           30|           32|     0.53900|    0.75690|
| 1     | 30     |         32| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | None                                                     | TRUE               |                     0.50|          0.25|                   32|           30|           16|     0.51382|    0.75650|
| 10    | 30     |         48| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.25|          0.50|                   16|           30|           32|     0.56308|    0.75610|
| 17    | 30     |         32| Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.25|          0.50|                   32|           30|           16|     0.53450|    0.75448|
| 14    | 30     |         64| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.25|          0.50|                   48|           30|           16|     0.55132|    0.75366|
| 39    | 30     |         48| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.00|          0.50|                   32|           30|            8|     0.52840|    0.75286|
| 30    | 30     |         32| Adadelta   | None                                                     | None                                                     | TRUE               |                     0.50|          0.25|                   32|           30|           32|     0.53156|    0.74960|
| 15    | 30     |         32| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.25|          0.50|                   32|           30|            8|     0.54960|    0.74674|
| 20    | 30     |         48| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.50|          0.00|                   48|           30|           16|     0.55758|    0.74512|
| 33    | 30     |         64| Adam       | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.50|          0.25|                   48|           30|           32|     0.55600|    0.74470|
| 29    | 30     |         64| Adam       | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.50|          0.00|                   32|           30|           16|     0.57122|    0.73822|
| 16    | 30     |         32| Adam       | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.50|          0.50|                   48|           30|           16|     0.56358|    0.73700|
| 46    | 30     |         32| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.50|          0.25|                   48|           30|            8|     0.56045|    0.73680|
| 19    | 30     |         64| Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | FALSE              |                     0.25|          0.50|                   16|           30|           32|     0.59176|    0.71096|
| 26    | 30     |         32| sgd        | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.50|          0.00|                   16|           30|            8|     0.76946|    0.54108|
| 11    | 30     |         48| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | None                                                     | TRUE               |                     0.25|          0.50|                   32|           30|            8|     0.99888|    0.53294|
| 44    | 30     |         64| sgd        | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.50|          0.50|                   32|           30|           16|     1.87962|    0.52926|
| 3     | 30     |         64| sgd        | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.50|          0.00|                   32|           30|           16|     0.85548|    0.52278|
| 35    | 30     |         64| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | None                                                     | FALSE              |                     0.00|          0.00|                   16|           30|           16|     0.72200|    0.51992|
| 45    | 30     |         32| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | None                                                     | TRUE               |                     0.25|          0.00|                   48|           30|            8|     0.78138|    0.51788|
| 22    | 30     |         64| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.25|          0.25|                   32|           30|           32|     1.38884|    0.51544|
| 4     | 30     |         48| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.25|          0.25|                   32|           30|           32|     1.91122|    0.51464|
| 43    | 30     |         64| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | TRUE               |                     0.25|          0.50|                   48|           30|           32|     3.07908|    0.51424|
| 36    | 30     |         48| sgd        | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | FALSE              |                     0.50|          0.00|                   48|           30|           32|     0.84310|    0.51344|
| 24    | 30     |         32| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | FALSE              |                     0.00|          0.00|                   16|           30|           16|     1.31862|    0.51180|
| 34    | 30     |         48| sgd        | &lt;keras.regularizers.L1L2 object at 0x7f5877eddbe0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f58926fc550&gt; | TRUE               |                     0.00|          0.50|                   32|           30|           32|     1.36012|    0.51018|

#### CNN and GRU

``` r
src_file_cnn_gru <- "../../../results/logs/0009-cnn-gru-scan.log"
src_params_list_cnn_gru <- c(
  "optimizers", "opt_recurrent_regs", "opt_kernel_regs", "opt_go_backwards",
  "opt_dropout_recurrent", "opt_dropout", "maxpool_size", "kernel_size",
  "gru_hidden_units", "filter_conv", "epochs_raw", "batch_size", "activation_conv"
)
src_numeric_params_cnn_gru <- c(
  "opt_dropout_recurrent", "opt_dropout", "maxpool_size", "kernel_size",
  "gru_hidden_units", "filter_conv", "epochs_raw", "batch_size"
)
```

``` r
data_cnn_gru <- src_file_cnn_gru %>%
  read_log_into_df_with_params_list(src_params_list_cnn_gru, src_numeric_params_cnn_gru) %>%
  clean_log_df_with_params()

new_data_cnn_gru <- data_cnn_gru %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_gru)
```

``` r
new_data_cnn_gru %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  knitr::kable()
```

| model | epochs | optimizers | opt\_recurrent\_regs                                     | opt\_kernel\_regs                                        | opt\_go\_backwards |  opt\_dropout\_recurrent|  opt\_dropout|  maxpool\_size|  kernel\_size|  gru\_hidden\_units|  filter\_conv|  epochs\_raw|  batch\_size| activation\_conv |  loss\_mean|  acc\_mean|
|:------|:-------|:-----------|:---------------------------------------------------------|:---------------------------------------------------------|:-------------------|------------------------:|-------------:|--------------:|-------------:|-------------------:|-------------:|------------:|------------:|:-----------------|-----------:|----------:|
| 2     | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | TRUE               |                     0.00|          0.00|              3|             2|                  64|            48|           30|            8| relu             |     0.26848|    0.93946|
| 17    | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | None                                                     | TRUE               |                     0.00|          0.00|              3|             2|                  16|            16|           30|            8| None             |     0.19052|    0.93778|
| 14    | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | TRUE               |                     0.00|          0.00|              3|             2|                   8|            32|           30|           16| relu             |     0.28720|    0.92236|
| 18    | 30     | Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | FALSE              |                     0.00|          0.00|              3|             2|                  32|            48|           30|           32| None             |     0.37016|    0.87398|
| 15    | 30     | Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | FALSE              |                     0.50|          0.50|              3|             3|                  32|            48|           30|            8| None             |     0.36650|    0.85406|
| 21    | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | None                                                     | FALSE              |                     0.00|          0.00|              2|             3|                   8|             8|           30|           16| None             |     0.39555|    0.83640|
| 8     | 30     | Adam       | None                                                     | None                                                     | TRUE               |                     0.50|          0.00|              2|             2|                  32|            32|           30|           32| None             |     0.39266|    0.83092|
| 13    | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | TRUE               |                     0.00|          0.50|              3|             3|                  32|             8|           30|            8| relu             |     0.49548|    0.77478|
| 11    | 30     | Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | TRUE               |                     0.25|          0.50|              3|             3|                   8|            16|           30|            8| None             |     0.51360|    0.77318|
| 19    | 30     | Adadelta   | None                                                     | None                                                     | FALSE              |                     0.25|          0.25|              3|             2|                   8|            32|           30|           32| relu             |     0.48564|    0.77318|
| 3     | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | TRUE               |                     0.25|          0.25|              3|             3|                  64|            16|           30|           32| relu             |     0.53742|    0.75608|
| 12    | 30     | Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | FALSE              |                     0.25|          0.25|              2|             2|                  64|             8|           30|            8| relu             |     0.52394|    0.75326|
| 5     | 30     | Adam       | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | FALSE              |                     0.50|          0.25|              2|             3|                   8|            16|           30|           32| relu             |     0.55996|    0.73782|
| 4     | 30     | Adadelta   | None                                                     | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | FALSE              |                     0.25|          0.50|              2|             2|                  16|             8|           30|           32| None             |     0.58524|    0.70896|
| 6     | 30     | Adadelta   | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | TRUE               |                     0.00|          0.25|              2|             2|                  64|            16|           30|           16| relu             |     0.59216|    0.70568|
| 20    | 30     | sgd        | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | FALSE              |                     0.00|          0.00|              2|             3|                   8|            48|           30|            8| relu             |     1.05122|    0.68456|
| 1     | 30     | Adam       | None                                                     | None                                                     | FALSE              |                     0.50|          0.50|              3|             2|                  16|             8|           30|           16| relu             |     0.59692|    0.67966|
| 10    | 30     | sgd        | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | FALSE              |                     0.25|          0.25|              2|             3|                  16|            16|           30|            8| relu             |     0.75486|    0.57438|
| 16    | 30     | sgd        | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | None                                                     | TRUE               |                     0.00|          0.25|              3|             3|                   8|            32|           30|           32| None             |     0.70458|    0.55204|
| 7     | 30     | sgd        | None                                                     | None                                                     | FALSE              |                     0.00|          0.50|              3|             2|                  64|            16|           30|            8| None             |     0.68772|    0.54312|
| 9     | 30     | sgd        | &lt;keras.regularizers.L1L2 object at 0x7f2b15b4c8d0&gt; | &lt;keras.regularizers.L1L2 object at 0x7f2b04183080&gt; | FALSE              |                     0.00|          0.50|              2|             2|                   8|             8|           30|            8| None             |     0.78916|    0.53088|

#### CNN LSTM

``` r
# File 5
src_file_cnn_lstm <- "../../../results/logs/0010-getting-started-cnn-lstm-scan.log"
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
data_cnn_lstm <- src_file_cnn_lstm %>%
  read_log_into_df_with_params_list(src_params_list_cnn_lstm, src_numeric_params_cnn_lstm) %>%
  clean_log_df_with_params()
```

    ## Warning: Expected 7 pieces. Additional pieces discarded in 386 rows [869,
    ## 870, 871, 872, 1737, 1738, 1739, 1768, 1769, 1799, 1800, 1803, 2667, 2668,
    ## 2729, 2730, 2731, 2732, 2760, 3442, ...].

``` r
new_data_cnn_lstm <- data_cnn_lstm %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list_cnn_lstm)
```

``` r
new_data_cnn_lstm %>% 
  select(-c(loss_sd, acc_sd)) %>% 
  arrange(desc(acc_mean)) %>% 
  knitr::kable()
```

| model | epochs |  strides| padding | optimizers |  number\_hidden\_units|  filters\_LSTM|  filters|  epochs\_raw| bn  |  batch\_size| activation\_convolution | activation\_LSTM |  loss\_mean|  acc\_mean|
|:------|:-------|--------:|:--------|:-----------|----------------------:|--------------:|--------:|------------:|:----|------------:|:------------------------|:-----------------|-----------:|----------:|
| 7     | 30     |        1| valid   | Adam       |                     32|             16|       16|           30| yes |           32| None                    | tanh             |     0.02108|  0.9967400|
| 1     | 30     |        1| valid   | Adadelta   |                     32|             32|       32|           30| yes |           16| None                    | tanh             |     0.00976|  0.9963200|
| 17    | 30     |        1| valid   | Adadelta   |                     16|              8|       32|           30| yes |           16| None                    | tanh             |     0.01708|  0.9959400|
| 3     | 30     |        1| valid   | Adadelta   |                     16|             48|       48|           30| yes |           16| None                    | tanh             |     0.00884|  0.9959200|
| 15    | 30     |        1| valid   | Adam       |                     16|             16|       32|           30| yes |           16| None                    | tanh             |     0.01144|  0.9959200|
| 12    | 30     |        1| valid   | Adam       |                     32|             16|       48|           30| yes |           32| None                    | tanh             |     0.01190|  0.9955400|
| 6     | 30     |        1| valid   | Adam       |                     16|              8|       32|           30| yes |           32| None                    | tanh             |     0.02482|  0.9955200|
| 5     | 30     |        1| valid   | Adam       |                     32|             16|       32|           30| yes |           16| None                    | tanh             |     0.01002|  0.9951200|
| 10    | 30     |        1| valid   | Adam       |                      8|             32|       16|           30| yes |           16| None                    | tanh             |     0.01170|  0.9951200|
| 11    | 30     |        1| valid   | Adam       |                     16|             16|       32|           30| yes |           32| None                    | tanh             |     0.01576|  0.9951000|
| 21    | 30     |        1| valid   | Adadelta   |                      8|             32|       16|           30| yes |           32| None                    | tanh             |     0.01790|  0.9945667|
| 9     | 30     |        1| valid   | Adam       |                      8|             32|       32|           30| yes |           16| None                    | tanh             |     0.01060|  0.9943200|
| 16    | 30     |        1| valid   | Adadelta   |                     32|              8|       48|           30| yes |           32| None                    | tanh             |     0.01840|  0.9943000|
| 19    | 30     |        1| valid   | Adam       |                      8|             16|       16|           30| yes |           16| None                    | tanh             |     0.01808|  0.9943000|
| 18    | 30     |        1| valid   | Adam       |                     32|              8|       48|           30| yes |           16| None                    | tanh             |     0.03930|  0.9898400|
| 4     | 30     |        1| valid   | sgd        |                     16|             16|       48|           30| yes |           16| None                    | tanh             |     0.21496|  0.9528600|
| 8     | 30     |        1| valid   | sgd        |                     32|             16|       32|           30| yes |           16| None                    | tanh             |     0.34990|  0.8849600|
| 14    | 30     |        1| valid   | sgd        |                     32|             48|       48|           30| yes |           32| None                    | tanh             |     0.41180|  0.8699000|
| 2     | 30     |        1| valid   | sgd        |                     16|              8|       48|           30| yes |           32| None                    | tanh             |     0.47388|  0.8292600|
| 20    | 30     |        1| valid   | sgd        |                     16|             16|       16|           30| yes |           16| None                    | tanh             |     0.47260|  0.7963200|
| 13    | 30     |        1| valid   | sgd        |                      8|              8|       32|           30| yes |           32| None                    | tanh             |     0.55120|  0.7719400|
