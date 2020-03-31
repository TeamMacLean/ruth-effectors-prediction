Pretrained Models (After Resularizers) that are used for the Ensemble
=====================================================================

Load library
------------

``` r
get_plot <- function(train_data, val_data) {
  
  # Get data join
  join_data <- train_data %>% 
    left_join(., val_data, by = "V1") %>% 
    `colnames<-`(c("epochs", "acc_train", "loss_train", "acc_val", "loss_val"))
  
  # Get Acc and Val data individually
  
  # Get acc data
  acc_data <- join_data %>% 
    select(c("epochs", "acc_train", "acc_val")) 

  acc_data <- melt(acc_data, id=c("epochs")) %>% 
    `colnames<-`(c("epochs", "acc", "value")) %>% 
    mutate(epochs = epochs + 1)

  # Get loss data
  loss_data <- join_data %>% 
    select(c("epochs", "loss_train", "loss_val")) 

  loss_data <- melt(loss_data, id=c("epochs")) %>% 
    `colnames<-`(c("epochs", "loss", "value")) %>% 
    mutate(epochs = epochs + 1)
  
  
  # Plot the data
  gg_acc <- ggplot(acc_data) +
    aes(x = epochs, y = value, group = acc, color = as.factor(acc)) +
    geom_line() +
    geom_point(size = 1.5) +
    labs(x = "Epochs", y = "Accuracy") +
    scale_color_manual(
    values = c("acc_train" = "salmon", "acc_val" = "seagreen3"),
    labels = c("acc_train" = "train", "acc_val" = "val")
    )

  gg_loss <- ggplot(loss_data) +
    aes(x = epochs, y = value, group = loss, color = as.factor(loss)) +
    geom_line() +
    geom_point(size = 1.5) +
    labs(x = "Epochs", y = "Loss", color = "Data") +
    scale_color_manual(
    values = c("loss_train" = "salmon", "loss_val" = "seagreen3"),
    labels = c("loss_train" = "train", "loss_val" = "val")
    )

  plot_together <- GGally::ggmatrix(
    plots = list(gg_loss, gg_acc),
    nrow = 2,
    ncol = 1,
    xlab = "Epochs",
    yAxisLabels = c("Loss", "Accuracy"),
    legend = 1
    ) +
    theme(strip.placement = "outside")
  
  return(plot_together)
}
```

CNN-LSTM Models
---------------

``` r
train_cnn_lstm <- data.table::fread("../../../../results/saved_models_for_ensemble/cnn_lstm/history_results/df_results_train_saved_model.csv")
val_cnn_lstm <- data.table::fread("../../../../results/saved_models_for_ensemble/cnn_lstm/history_results/df_results_val_cnn_lstm_saved_model.csv")
```

``` r
get_plot(train_cnn_lstm, val_cnn_lstm)
```

![](/Users/kristian/Documents/Workspace/ruth-effectors-prediction/reports/model-scripts/binary-class/0014_pretrained_models_files/figure-markdown_github/unnamed-chunk-4-1.png)

CNN-GRU
-------

``` r
train_cnn_gru <- data.table::fread("../../../../results/saved_models_for_ensemble/cnn-gru/history_results/df_results_train_cnn_gru_saved_model1.csv")
val_cnn_gru <- data.table::fread("../../../../results/saved_models_for_ensemble/cnn-gru/history_results/df_results_val_cnn_gru_saved_model1.csv")
```

``` r
get_plot(train_cnn_gru, val_cnn_gru)
```

![](/Users/kristian/Documents/Workspace/ruth-effectors-prediction/reports/model-scripts/binary-class/0014_pretrained_models_files/figure-markdown_github/unnamed-chunk-6-1.png)

LSTM - Embedding
----------------

``` r
train_lstm_emb <- data.table::fread("../../../../results/saved_models_for_ensemble/lstm_emb/history_results/df_results_train_lstm_emb_saved_model1.csv")
val_lstm_emb <- data.table::fread("../../../../results/saved_models_for_ensemble/lstm_emb/history_results/df_results_val_lstm_emb_saved_model1.csv")
```

``` r
get_plot(train_lstm_emb, val_lstm_emb)
```

![](/Users/kristian/Documents/Workspace/ruth-effectors-prediction/reports/model-scripts/binary-class/0014_pretrained_models_files/figure-markdown_github/unnamed-chunk-8-1.png)

GRU - Embedding
---------------

``` r
train_gru_emb <- data.table::fread("../../../../results/saved_models_for_ensemble/gru_emb/history_results/df_results_train_gru_emb_saved_model1.csv")
val_gru_emb <- data.table::fread("../../../../results/saved_models_for_ensemble/gru_emb/history_results/df_results_val_gru_emb_saved_model1.csv")
```

``` r
get_plot(train_gru_emb, val_gru_emb)
```

![](/Users/kristian/Documents/Workspace/ruth-effectors-prediction/reports/model-scripts/binary-class/0014_pretrained_models_files/figure-markdown_github/unnamed-chunk-10-1.png)
