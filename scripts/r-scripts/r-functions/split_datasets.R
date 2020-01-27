# Load library

library(docstring)
library(caret)

# Documentation of docstring https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html

#  Function to get the data


get_data_labeled_binary <- function(dataset_pos,
                                    dataset_neg) {
  #' Label the dataset for the binary classification
  #'
  #' @description this function can be used to add binary label (which either 1 or 0)
  #'
  #' @param dataset_pos dataframe. Set of dataframe that will be labeled 1
  #' @param dataset_neg dataframe. Set of dataframe that will be labeled 0

  # Labelling the data, eff = 1 and non_eff = 0
  dataset_pos_with_label <- dataset_pos %>%
    dplyr::mutate(label = 1)

  dataset_neg_with_label <- dataset_neg %>%
    dplyr::mutate(label = 0)

  full_dataset <- rbind(dataset_pos_with_label, dataset_neg_with_label)

  return(full_dataset)
}

get_data_splitted <- function(data_to_split, p1, p2, test_dataset = FALSE) {
  #' Split data into different datasets
  #'
  #' @description this function can be used to split the data into different datasets with the specific proportion
  #'
  #' @param data_to_split dataframe. Data that will be splitted, it is dataframe
  #' @param p1 int. Percentage of data goes to training
  #' @param p2 int. Percentage of data goes to validation
  #' @param test_dataset boolean. Percentage of data goes to test

  # Splitting the data
  # Create index for testing and training data
  training_id <- createDataPartition(
    y = data_to_split$label,
    p = p1,
    list = FALSE
  )

  # subset power_plant data to training
  training <- data_to_split[training_id, ]

  # subset the rest to test
  rest <- data_to_split[-training_id, ]

  # Splitting the rest into two different class of datasets
  val_id <- createDataPartition(
    y = rest$label,
    p = p2 / (1 - p1), list = FALSE
  )

  validation <- rest[val_id, ]

  dataset_list <- list(
    "training" = training,
    "validation" = validation
  )

  if (test_dataset) {
    testing <- rest[-val_id, ]
    dataset_list[["testing"]] <- testing
  }

  return(dataset_list)
}
