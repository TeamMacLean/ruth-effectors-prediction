# Define a function to get all accuracy of the models
get_all_acc <- function(data, true_label){

  # Change all of the label into factor
  data <- data %>%
    select(-c(sequence)) %>%
    mutate_each(list(as.factor))

  # Get the number of column
  num_col <- ncol(data)

  # Initialize the list of accuracy
  list_acc <- numeric(length = num_col)

  # For loop in getting acc for each models
  for (i in 1:num_col){
    pred_each_model <- data %>%
      pull(colnames(data)[i])

    tab <- table(true_label %>%
                   pull(),
                 pred_each_model)

    acc_each_model <- confusionMatrix(tab)$overall["Accuracy"]

    list_acc[i] <- acc_each_model
  }

  # Turn the list into dataframe

  df_acc  <- data.frame(matrix(unlist(list_acc), ncol=length(list_acc), byrow = F)) %>%
    `colnames<-`(c(colnames(data)))

  return(df_acc)
}
