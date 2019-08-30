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

#### Results

In summary, here is the table of the results

| Model    | Accuracy |
|----------|----------|
| GRU      | 0.68     |
| LSTM     | 0.693    |
| CNN-GRU  | 0.626    |
| CNN-LSTM | 0.693    |
