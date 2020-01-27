library(tidyverse)


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


# Select data files ----

# File 1
src_file <- "data/22sept_logs/0004-getting-started-cnn-lstm.log"
src_params_list <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw", "bn",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)

# File 2
src_file <- "data/22sept_logs/0007-gru-embedding_scan.log"
src_params_list <- c(
  "outputdim", "optimizers", "opt_recurrent_regs", "opt_kernel_regs",
  "opt_go_backwards", "opt_dropout_recurrent", "opt_dropout", "gru_hidden_units",
  "epochs_raw", "batch_size"
)
src_numeric_params <- c(
  "outputdim", "opt_dropout_recurrent", "opt_dropout",
  "gru_hidden_units", "epochs_raw", "batch_size"
)

# File 3
src_file <- "data/22sept_logs/0008-lstm-embedding-scan.log"
src_params_list <- c(
  "outputdim", "optimizers", "opt_recurrent_regs", "opt_kernel_regs",
  "opt_go_backwards", "opt_dropout_recurrent", "opt_dropout",
  "lstm_hidden_units", "epochs_raw", "batch_size"
)
src_numeric_params <- c(
  "outputdim", "opt_dropout_recurrent", "opt_dropout",
  "lstm_hidden_units", "epochs_raw", "batch_size"
)

# File 4
src_file <- "data/22sept_logs/0009-cnn-gru-scan.log"
src_params_list <- c(
  "optimizers", "opt_recurrent_regs", "opt_kernel_regs", "opt_go_backwards",
  "opt_dropout_recurrent", "opt_dropout", "maxpool_size", "kernel_size",
  "gru_hidden_units", "filter_conv", "epochs_raw", "batch_size", "activation_conv"
)
src_numeric_params <- c(
  "opt_dropout_recurrent", "opt_dropout", "maxpool_size", "kernel_size",
  "gru_hidden_units", "filter_conv", "epochs_raw", "batch_size"
)

# File 5
src_file <- "data/22sept_logs/0010-getting-started-cnn-lstm-scan.log"
src_params_list <- c(
  "strides", "padding", "optimizers", "number_hidden_units",
  "filters_LSTM", "filters", "epochs_raw", "bn",
  "batch_size", "activation_convolution", "activation_LSTM"
)
src_numeric_params <- c(
  "strides", "number_hidden_units", "filters_LSTM",
  "filters", "epochs_raw", "batch_size"
)


# Calling the functions ----

data <- src_file %>%
  read_log_into_df_with_params_list(src_params_list, src_numeric_params) %>%
  clean_log_df_with_params()

new_data <- data %>%
  tidyr::drop_na() %>%
  summarise_log_data_with_params_list(src_params_list)


# Plots ----

ggplot(new_data) +
  geom_histogram(aes(acc_mean, fill = epochs), bins = 15) +
  facet_wrap(~epochs) +
  labs(x = "Mean accuracy", y = "Count", fill = "Epochs")

ggplot(new_data) +
  geom_density(aes(acc_mean, fill = epochs), alpha = 0.5) +
  labs(x = "Mean accuracy", y = "Density", fill = "Epochs")

ggplot(new_data) +
  geom_violin(aes(x = epochs, y = acc_mean, fill = epochs)) +
  geom_boxplot(aes(x = epochs, y = acc_mean), width = 0.05) +
  labs(x = "Epochs", y = "Mean accuracy", fill = "Epochs")

ggplot(new_data) +
  geom_tile(aes(x = blocks, y = epochs, fill = acc_mean)) +
  facet_wrap(~ filters) +
  viridis::scale_fill_viridis() +
  coord_equal()

block_labs <- unique(new_data$blocks) %>%
  paste("Blocks") %>%
  `names<-`(unique(new_data$blocks))

filter_labs <- unique(new_data$filters) %>%
  paste("Filters") %>%
  `names<-`(unique(new_data$filters))

new_data %>%
  dplyr::filter(acc_mean > 0.59) %>%
  ggplot() +
  geom_violin(aes(x = epochs, y = acc_mean, fill = epochs)) +
  geom_boxplot(aes(x = epochs, y = acc_mean), width = 0.05) +
  facet_grid(
    blocks ~ filters,
    labeller = labeller(
      filters = filter_labs,
      blocks = block_labs
    )
  ) +
  labs(x = "Epochs", y = "Mean accuracy", fill = "Epochs")
