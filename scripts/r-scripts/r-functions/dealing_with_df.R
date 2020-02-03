# Collection of a function to dp data wraggling with library(dplyr)



# Function to get the number of column in R

get_nrow <- function(df, column_to_filter, filter_name){

  nrow_organism <- df %>%
    dplyr::filter({{ column_to_filter }} == filter_name) %>%
    nrow()

  return(nrow_organism)
}
