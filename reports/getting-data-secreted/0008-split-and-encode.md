Data processing of new dataset (secreted protein as negative datasets)
======================================================================

Background
----------

In this report, the process in getting all of new sequences splitted
into several datsets and encoded will be shown. The mian function below
uses several different functions from `split_datasets.R` and
`blast_data.R`.

Execution
---------

### Load libraries

``` r
library(docstring)
library(tidyverse)
# library(taxize)
library(caret)
# reticulate::use_condaenv(condaenv = "tensorflow2", conda = "/anaconda3/bin/conda")
```

``` r
# Get the source of the data
source(here::here("scripts/r-scripts/r-functions", "split_datasets.R"))
source(here::here("scripts/r-scripts/r-functions", "blast_data.R"))
```

Functions
---------

``` r
# Created a function that can autmomatically get the data and blast them until there is no identical datafor each datasets
split_data_without_identical <- function(dataset, p1, p2, test_data = TRUE, while_limit = 10, path_fasta = "data/secreted_data/split-blast/fasta_files", path_blast = "data/secreted_data/split-blast/blast_files", init_seed = 11) {
  # Set seed for reproducibility
  set.seed(init_seed)
  
  # Parse dataset name
  dataset_name <- deparse(substitute(dataset)) %>%
    stringr::str_split("_") %>%
    .[[1]] %>%
    .[1]

  #cat(dataset_name)
  #cat("\n")

  # Initialze variables for while loop
  n_identical_val <- 1
  n_identical_test <- 1
  while_step <- 1

  while (!(n_identical_val == 0 & n_identical_test == 0) & while_step < while_limit) {
    cat("Step", n_identical_val, n_identical_test, while_step, "\n")

    # Get data splitted
    list_datasets <- get_data_splitted(dataset, p1, p2, test_data = TRUE)

    # Change dataframe to fasta
    for (split in c("training", "validation", "testing")) {
      get_fasta_from_df(list_datasets[[split]], fasta_name = paste0(dataset_name, "_", split), dir_path = path_fasta)
      assign(
        x = paste0("fasta_path_", split),
        value = here::here(path_fasta, paste0(dataset_name, "_", split, ".fasta"))
      )
    }

    # Blast training vs validation datasets
    get_blast_data(database_fasta_path = fasta_path_training, query_fasta_path = fasta_path_validation, dir_path = path_blast)

    # Blast training vs testing datasets
    get_blast_data(database_fasta_path = fasta_path_training, query_fasta_path = fasta_path_testing, dir_path = path_blast)

    blast_train_vs_val <- blast_results(
      result_path = here::here(dir_path = path_blast, paste0(dataset_name, "_training", "_vs_", dataset_name, "_validation.tsv")),
      percent_threshold = 90
      )[["df_identical_protein"]]
    
    blast_train_vs_test <- blast_results(
      result_path = here::here(dir_path = path_blast, paste0(dataset_name, "_training", "_vs_", dataset_name, "_testing.tsv")),
      percent_threshold = 90
      )[["df_identical_protein"]]

    # Check wether they have identical data or not
    n_identical_val <- nrow(blast_train_vs_val)
    n_identical_test <- nrow(blast_train_vs_test)

    while_step <- while_step + 1
  }

  if (n_identical_val == 0 & n_identical_test == 0) {
    print("Good dataset split has been found without identical protein!")
    return(list_datasets)
  }
}
```

Load the data
-------------

### Get effector data for each class

``` r
get_seq_each_class <- function(df_effector, class_var) {
  df_seq <- df_effector %>%
    dplyr::filter(class == class_var) %>%
    dplyr::select(c(ProteinID, Sequence))

  return(df_seq)
}
```

``` r
effector_final_after_blast <- readRDS("../../../data/secreted_data/data_processed_after_signalp/effector_final_after_blast.RDS")

effector_seq_fungi <- get_seq_each_class(effector_final_after_blast, class_var = "fungi")
effector_seq_bacteria <- get_seq_each_class(effector_final_after_blast, class_var = "bacteria")
effector_seq_oomycete <- get_seq_each_class(effector_final_after_blast, class_var = "oomycete")
```

### Get the non-effector from randomly sampling data

``` r
non_effector_seq_fungi <- readRDS("../../../data/secreted_data/updated_signalp_results/fungi_non_eff_secreted_data.RDS") %>% 
  dplyr::select(c(ID, sequence)) %>% 
  `colnames<-`(c("ProteinID", "Sequence"))
```

``` r
non_effector_seq_bacteria <- readRDS("../../../data/secreted_data/updated_signalp_results/bacteria_non_eff_secreted_data.RDS")  %>% 
  dplyr::select(c(ID, sequence)) %>% 
  `colnames<-`(c("ProteinID", "Sequence"))
```

``` r
non_effector_seq_oomycete <- readRDS("../../../data/secreted_data/updated_signalp_results/oomycete_non_eff_secreted_data.RDS")  %>% 
  dplyr::select(c(ID, sequence)) %>% 
  `colnames<-`(c("ProteinID", "Sequence"))
```

Processing the data
-------------------

``` r
effector_seq_fungi <- readRDS("effector_seq_fungi.RDS")
non_effector_seq_fungi <- readRDS("noneffector_seq_fungi.RDS")

fungi_full_datasets <- get_data_labeled_binary(effector_seq_fungi, non_effector_seq_fungi)

fungi_splitted <- get_data_splitted(fungi_full_datasets, p1 = 0.6, p2 = 0.2, test_dataset = TRUE)
```

``` r
fungi_final_split <- split_data_without_identical(
  dataset = fungi_full_datasets, 
  p1 = 0.6, 
  p2 = 0.2,
  test_data = TRUE,
  while_limit = 25, 
  init_seed = 2906
)
```

    ## Step 1 1 1 
    ## Step 2 4 2 
    ## Step 4 0 3 
    ## Step 0 2 4 
    ## Step 1 3 5 
    ## Step 2 4 6 
    ## Step 1 2 7 
    ## [1] "Good dataset split has been found without identical protein!"

``` r
bacteria_full_datasets <- get_data_labeled_binary(effector_seq_bacteria, non_effector_seq_bacteria)

bacteria_final_split <- split_data_without_identical(
  dataset = bacteria_full_datasets, 
  p1 = 0.6, 
  p2 = 0.2,
  test_data = TRUE,
  while_limit = 25, 
  init_seed = 11
)
```

    ## Step 1 1 1 
    ## Step 5 7 2 
    ## Step 7 10 3 
    ## Step 6 13 4 
    ## Step 9 5 5 
    ## Step 3 11 6 
    ## Step 5 8 7 
    ## Step 6 6 8 
    ## Step 5 8 9 
    ## Step 5 5 10 
    ## Step 8 2 11 
    ## Step 3 10 12 
    ## Step 6 8 13 
    ## Step 5 9 14 
    ## Step 4 7 15 
    ## Step 10 3 16 
    ## Step 9 5 17 
    ## Step 7 5 18 
    ## Step 8 7 19 
    ## Step 3 3 20 
    ## Step 3 4 21 
    ## Step 1 10 22 
    ## Step 6 6 23 
    ## Step 11 5 24

``` r
oomycete_full_datasets <- get_data_labeled_binary(effector_seq_oomycete, non_effector_seq_oomycete)

oomycete_final_split <- split_data_without_identical(
  dataset = oomycete_full_datasets, 
  p1 = 0.6, 
  p2 = 0.2,
  test_data = TRUE,
  while_limit = 25, 
  init_seed = 2906
)
```

    ## Step 1 1 1 
    ## Step 1 0 2 
    ## Step 1 2 3 
    ## Step 1 1 4 
    ## Step 1 0 5 
    ## Step 1 0 6 
    ## [1] "Good dataset split has been found without identical protein!"
