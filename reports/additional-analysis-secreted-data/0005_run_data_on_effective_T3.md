Evaluate our data on the existed tool EffectiveT3
=================================================

Introduction
------------

Functions
---------

``` r
get_data_ready <- function(data){
  
  data <- data %>% 
  dplyr::select(-c("description", "score")) %>% 
  `colnames<-`(c("identifier", "prediction")) %>% 
  dplyr::mutate(prediction = case_when(
    prediction == "TRUE" ~ 1, 
    TRUE ~ 0
  )) %>% 
  dplyr::mutate(label = stringr::str_remove_all(identifier, ".*_"), 
                identifier =  stringr::str_extract(identifier, "[^_]*"))
  
  return(data)
}
```

``` r
calculate_accuracy <- function(data) {
  
  tab <- table(data %>% 
                 dplyr::select(prediction) %>% 
                 pull(), 
               data %>% 
                 dplyr::select(label) %>% 
                 pull())

  # Calculate acc
  acc <- caret::confusionMatrix(tab)$overall["Accuracy"]
  
  return(acc)
}
```

### Get prediction

``` r
bacteria_training_preds <- get_pred_effective_t3("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/data/secreted_data/ready_to_process/fasta_files/bacteria_training.fasta")
bacteria_validation_preds <- get_pred_effective_t3("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/data/secreted_data/ready_to_process/fasta_files/bacteria_validation.fasta")
bacteria_testing_preds <- get_pred_effective_t3("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/data/secreted_data/ready_to_process/fasta_files/bacteria_testing.fasta")
oomycete_testing_preds <- get_pred_effective_t3("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/data/secreted_data/ready_to_process/fasta_files/oomycete_testing.fasta")
fungi_testing_preds <- get_pred_effective_t3("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/data/secreted_data/ready_to_process/fasta_files/fungi_testing.fasta")
```

### Bacteria prediction results

``` r
data.frame(data = c("bacteria_train", 
                    "bacteria_val", 
                    "bacteria_test"), 
           acc_deepT3 = c(calculate_accuracy(get_data_ready(bacteria_training_preds)),
                          calculate_accuracy(get_data_ready(bacteria_validation_preds)),
                          calculate_accuracy(get_data_ready(bacteria_testing_preds)))) %>% 
  knitr::kable()
```

| data            |  acc\_deepT3|
|:----------------|------------:|
| bacteria\_train |    0.8333333|
| bacteria\_val   |    0.7763158|
| bacteria\_test  |    0.8289474|

### Fungi and oomycete results

``` r
data.frame(data = c("fungi_test", 
                    "oomycete_test"), 
           acc_deepT3 = c(calculate_accuracy(get_data_ready(fungi_testing_preds)),
                          calculate_accuracy(get_data_ready(oomycete_testing_preds)))) %>% 
  knitr::kable()
```

| data           |  acc\_deepT3|
|:---------------|------------:|
| fungi\_test    |    0.5526316|
| oomycete\_test |    0.6176471|
