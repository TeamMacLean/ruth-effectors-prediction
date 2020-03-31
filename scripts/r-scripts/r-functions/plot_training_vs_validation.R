# A function to get a acc plot from different

# Load libraries
library(tidyverse)
library(ggplot2)
library(reshape)
library(patchwork)
library(ggrepel)

# Function
get_gg_acc <- function(acc_train_data, acc_val_data, plot_title = "Accuracy for Training and Validation") {

  # Melt data
  acc_train_data_melt <- acc_train_data %>%
    dplyr::rename(epochs = V1) %>%
    tidyr::pivot_longer(
      cols = -c(epochs),
      names_to = "variable",
      values_to = "train"
    )

  acc_val_data_melt <- acc_val_data %>%
    dplyr::rename(epochs = V1) %>%
    tidyr::pivot_longer(
      cols = -c(epochs),
      names_to = "CV",
      values_to = "val"
    )

  # Rename the data column
  acc_train_data_melt <- acc_train_data_melt %>%
    `colnames<-`(c("epochs", "CV", "train"))

  acc_val_data_melt <- acc_val_data_melt %>%
    `colnames<-`(c("epochs", "CV", "val"))

  # Left join into 1 dataframe
  train_val_acc_data <- acc_train_data_melt %>%
    left_join(
      acc_val_data_melt,
      by = c("epochs", "CV")
    ) %>%
    mutate(epochs = epochs + 1) %>%
    tidyr::pivot_longer(
      cols = -c(epochs, CV),
      names_to = "dataset",
      values_to = "accuracy"
    ) %>%
    mutate(
      dataset = case_when(
        dataset == "val" ~ "validation",
        dataset == "train" ~ "training"
      )
    )

  # Plot the accuracy data
  gg_accuracy <- ggplot(train_val_acc_data) +
    aes(
      x = epochs,
      y = accuracy,
      color = dataset,
      group = dataset,
      linetype = dataset
    ) +
    geom_point(size = 1) +
    geom_line() +
    # scale_color_manual(values = c("validation" = "blue", "training" = "red")) +
    # viridis::scale_color_viridis(discrete = TRUE) +
    scale_color_viridis_d(begin = 0.65, end = 0.1) +
    scale_linetype_manual(values = c("validation" = "dashed", "training" = "solid")) +
    labs(title = plot_title, x = "Epochs", y = "Accuracy")

  return(gg_accuracy)
}




# get_gg_acc <- function(acc_train_data, acc_val_data, plot_title = "Accuracy for Training and Validation") {
#
#   # Melt data
#   acc_train_data_melt <- acc_train_data %>%
#     dplyr::rename(epochs = V1) %>%
#     tidyr::pivot_longer(
#       cols = -c(epochs),
#       names_to = "variable",
#       values_to = "train"
#     )
#
#   acc_val_data_melt <- acc_val_data %>%
#     dplyr::rename(epochs = V1) %>%
#     tidyr::pivot_longer(
#       cols = -c(epochs),
#       names_to = "CV",
#       values_to = "val"
#     )
#
#   # Rename the data column
#   acc_train_data_melt <- acc_train_data_melt %>%
#     `colnames<-`(c("epochs", "CV", "train"))
#
#   acc_val_data_melt <- acc_val_data_melt %>%
#     `colnames<-`(c("epochs", "CV", "val"))
#
#   # Left join into 1 dataframe
#   train_val_acc_data <- acc_train_data_melt %>%
#     left_join(
#       acc_val_data_melt,
#       by = c("epochs", "CV")
#     ) %>%
#     mutate(epochs = epochs + 1) %>%
#     tidyr::pivot_longer(
#       cols = -c(epochs, CV),
#       names_to = "dataset",
#       values_to = "accuracy"
#     )
#
#   # Plot the accuracy data
#   gg_accuracy <- ggplot(train_val_acc_data) +
#     aes(
#       x = epochs,
#       y = accuracy,
#       color = dataset,
#       group = dataset,
#       linetype = dataset
#     ) +
#     geom_point(size = 1) +
#     geom_line() +
#     scale_color_manual(
#       values = c("val" = "#481567FF", "train" = "#3CBB75FF"),
#       labels = c("val" = "validation", "train" = "training")
#     ) +
#     scale_linetype_manual(
#       values = c("val" = "dashed", "train" = "solid"),
#       labels = c("val" = "validation", "train" = "training")
#     ) +
#     labs(title = plot_title, x = "Epochs", y = "Accuracy")
#
#   return(gg_accuracy)
# }


plot_comparison <- function(data, x_var, y_var, group_var, show_label = TRUE, label_digits = 3) {
  plot_data <- data %>%
    dplyr::mutate_at(
      .vars = dplyr::vars({{ group_var }}),
      .funs = as.factor
    ) %>%
    dplyr::group_by({{ group_var }}) %>%
    dplyr::mutate(
      mean_x = mean({{ x_var }}, na.rm = TRUE),
      mean_y = mean({{ y_var }}, na.rm = TRUE)
    )

  gg <- plot_data %>%
    ggplot() +
    aes(x = {{ x_var }}, y = {{ y_var }}, group = {{ group_var }}, color = {{ group_var }}) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_hline(
      aes(yintercept = mean_y, col = {{ group_var }}),
      linetype = "dashed"
    ) +
    scale_color_viridis_d(begin = 0.75, end = 0)

  if (show_label) {
    gg <- gg +
      ggrepel::geom_label_repel(
        data = plot_data %>% group_by({{ group_var }}) %>% slice(1),
        aes(x = mean_x, y = mean_y, label = round(mean_y, label_digits)),
        seed = 1,
        show.legend = FALSE
      ) #+
    # facet_grid(dplyr::vars({{ group_var }}), scales = "free")
    # facet_wrap(dplyr::vars({{ group_var }}), scales = "free")
  }

  return(gg)
}
