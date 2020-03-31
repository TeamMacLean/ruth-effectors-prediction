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
    # cat("Step", n_identical_val, n_identical_test, while_step, "\n")

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
      percent_threshold = 95
      )[["df_identical_protein"]]
    
    blast_train_vs_test <- blast_results(
      result_path = here::here(dir_path = path_blast, paste0(dataset_name, "_training", "_vs_", dataset_name, "_testing.tsv")),
      percent_threshold = 95
      )[["df_identical_protein"]]

    # Check wether they have identical data or not
    n_identical_val <- nrow(blast_train_vs_val)
    n_identical_test <- nrow(blast_train_vs_test)
    
    cat("Step", while_step, ":", "n_identical_train_vs_val", n_identical_val, ",", "n_identical_train_vs_test", n_identical_test, "\n")
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

``` r
effector_seq_bacteria %>% saveRDS("effector_seq_bacteria.RDS")

non_effector_seq_bacteria %>%  saveRDS("non_effector_seq_bacteria.RDS")
```

Processing the data
-------------------

### Fungi

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

    ## Step 1 : n_identical_train_vs_val 2 , n_identical_train_vs_test 4 
    ## Step 2 : n_identical_train_vs_val 4 , n_identical_train_vs_test 0 
    ## Step 3 : n_identical_train_vs_val 0 , n_identical_train_vs_test 2 
    ## Step 4 : n_identical_train_vs_val 1 , n_identical_train_vs_test 3 
    ## Step 5 : n_identical_train_vs_val 2 , n_identical_train_vs_test 4 
    ## Step 6 : n_identical_train_vs_val 1 , n_identical_train_vs_test 2 
    ## Step 7 : n_identical_train_vs_val 0 , n_identical_train_vs_test 0 
    ## [1] "Good dataset split has been found without identical protein!"

### Bacteria

``` r
bacteria_full_datasets <- get_data_labeled_binary(effector_seq_bacteria, non_effector_seq_bacteria)

bacteria_final_split <- split_data_without_identical(
  dataset = bacteria_full_datasets, 
  p1 = 0.6, 
  p2 = 0.2,
  test_data = TRUE,
  while_limit = 250, 
  init_seed = 30
)
```

    ## Step 1 : n_identical_train_vs_val 1 , n_identical_train_vs_test 5 
    ## Step 2 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 3 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 4 : n_identical_train_vs_val 9 , n_identical_train_vs_test 2 
    ## Step 5 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 6 : n_identical_train_vs_val 3 , n_identical_train_vs_test 8 
    ## Step 7 : n_identical_train_vs_val 1 , n_identical_train_vs_test 6 
    ## Step 8 : n_identical_train_vs_val 3 , n_identical_train_vs_test 6 
    ## Step 9 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 10 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 11 : n_identical_train_vs_val 3 , n_identical_train_vs_test 2 
    ## Step 12 : n_identical_train_vs_val 2 , n_identical_train_vs_test 9 
    ## Step 13 : n_identical_train_vs_val 7 , n_identical_train_vs_test 2 
    ## Step 14 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 15 : n_identical_train_vs_val 10 , n_identical_train_vs_test 2 
    ## Step 16 : n_identical_train_vs_val 6 , n_identical_train_vs_test 1 
    ## Step 17 : n_identical_train_vs_val 8 , n_identical_train_vs_test 0 
    ## Step 18 : n_identical_train_vs_val 5 , n_identical_train_vs_test 0 
    ## Step 19 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 20 : n_identical_train_vs_val 4 , n_identical_train_vs_test 7 
    ## Step 21 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 22 : n_identical_train_vs_val 5 , n_identical_train_vs_test 2 
    ## Step 23 : n_identical_train_vs_val 2 , n_identical_train_vs_test 11 
    ## Step 24 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 25 : n_identical_train_vs_val 3 , n_identical_train_vs_test 4 
    ## Step 26 : n_identical_train_vs_val 2 , n_identical_train_vs_test 5 
    ## Step 27 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 28 : n_identical_train_vs_val 9 , n_identical_train_vs_test 4 
    ## Step 29 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 30 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 31 : n_identical_train_vs_val 7 , n_identical_train_vs_test 2 
    ## Step 32 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 33 : n_identical_train_vs_val 1 , n_identical_train_vs_test 5 
    ## Step 34 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 35 : n_identical_train_vs_val 0 , n_identical_train_vs_test 9 
    ## Step 36 : n_identical_train_vs_val 1 , n_identical_train_vs_test 6 
    ## Step 37 : n_identical_train_vs_val 6 , n_identical_train_vs_test 6 
    ## Step 38 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 39 : n_identical_train_vs_val 2 , n_identical_train_vs_test 5 
    ## Step 40 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 41 : n_identical_train_vs_val 4 , n_identical_train_vs_test 7 
    ## Step 42 : n_identical_train_vs_val 0 , n_identical_train_vs_test 6 
    ## Step 43 : n_identical_train_vs_val 6 , n_identical_train_vs_test 5 
    ## Step 44 : n_identical_train_vs_val 7 , n_identical_train_vs_test 6 
    ## Step 45 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 46 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 47 : n_identical_train_vs_val 8 , n_identical_train_vs_test 3 
    ## Step 48 : n_identical_train_vs_val 3 , n_identical_train_vs_test 4 
    ## Step 49 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 50 : n_identical_train_vs_val 3 , n_identical_train_vs_test 3 
    ## Step 51 : n_identical_train_vs_val 2 , n_identical_train_vs_test 7 
    ## Step 52 : n_identical_train_vs_val 8 , n_identical_train_vs_test 4 
    ## Step 53 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 54 : n_identical_train_vs_val 9 , n_identical_train_vs_test 1 
    ## Step 55 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 56 : n_identical_train_vs_val 3 , n_identical_train_vs_test 1 
    ## Step 57 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 58 : n_identical_train_vs_val 2 , n_identical_train_vs_test 1 
    ## Step 59 : n_identical_train_vs_val 3 , n_identical_train_vs_test 8 
    ## Step 60 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 61 : n_identical_train_vs_val 3 , n_identical_train_vs_test 3 
    ## Step 62 : n_identical_train_vs_val 5 , n_identical_train_vs_test 3 
    ## Step 63 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 64 : n_identical_train_vs_val 7 , n_identical_train_vs_test 3 
    ## Step 65 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 66 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 67 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 68 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 69 : n_identical_train_vs_val 8 , n_identical_train_vs_test 3 
    ## Step 70 : n_identical_train_vs_val 6 , n_identical_train_vs_test 8 
    ## Step 71 : n_identical_train_vs_val 2 , n_identical_train_vs_test 4 
    ## Step 72 : n_identical_train_vs_val 2 , n_identical_train_vs_test 3 
    ## Step 73 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 74 : n_identical_train_vs_val 7 , n_identical_train_vs_test 1 
    ## Step 75 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 76 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 77 : n_identical_train_vs_val 1 , n_identical_train_vs_test 6 
    ## Step 78 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 79 : n_identical_train_vs_val 3 , n_identical_train_vs_test 8 
    ## Step 80 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 81 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 82 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 83 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 84 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 85 : n_identical_train_vs_val 3 , n_identical_train_vs_test 2 
    ## Step 86 : n_identical_train_vs_val 6 , n_identical_train_vs_test 5 
    ## Step 87 : n_identical_train_vs_val 3 , n_identical_train_vs_test 8 
    ## Step 88 : n_identical_train_vs_val 2 , n_identical_train_vs_test 8 
    ## Step 89 : n_identical_train_vs_val 3 , n_identical_train_vs_test 6 
    ## Step 90 : n_identical_train_vs_val 6 , n_identical_train_vs_test 4 
    ## Step 91 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 92 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 93 : n_identical_train_vs_val 9 , n_identical_train_vs_test 5 
    ## Step 94 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 95 : n_identical_train_vs_val 7 , n_identical_train_vs_test 5 
    ## Step 96 : n_identical_train_vs_val 2 , n_identical_train_vs_test 6 
    ## Step 97 : n_identical_train_vs_val 8 , n_identical_train_vs_test 2 
    ## Step 98 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 99 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 100 : n_identical_train_vs_val 2 , n_identical_train_vs_test 4 
    ## Step 101 : n_identical_train_vs_val 3 , n_identical_train_vs_test 6 
    ## Step 102 : n_identical_train_vs_val 8 , n_identical_train_vs_test 5 
    ## Step 103 : n_identical_train_vs_val 3 , n_identical_train_vs_test 10 
    ## Step 104 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 105 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 106 : n_identical_train_vs_val 7 , n_identical_train_vs_test 6 
    ## Step 107 : n_identical_train_vs_val 3 , n_identical_train_vs_test 2 
    ## Step 108 : n_identical_train_vs_val 8 , n_identical_train_vs_test 2 
    ## Step 109 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 110 : n_identical_train_vs_val 3 , n_identical_train_vs_test 4 
    ## Step 111 : n_identical_train_vs_val 5 , n_identical_train_vs_test 2 
    ## Step 112 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 113 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 114 : n_identical_train_vs_val 4 , n_identical_train_vs_test 7 
    ## Step 115 : n_identical_train_vs_val 4 , n_identical_train_vs_test 7 
    ## Step 116 : n_identical_train_vs_val 4 , n_identical_train_vs_test 2 
    ## Step 117 : n_identical_train_vs_val 10 , n_identical_train_vs_test 4 
    ## Step 118 : n_identical_train_vs_val 1 , n_identical_train_vs_test 4 
    ## Step 119 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 120 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 121 : n_identical_train_vs_val 1 , n_identical_train_vs_test 4 
    ## Step 122 : n_identical_train_vs_val 10 , n_identical_train_vs_test 2 
    ## Step 123 : n_identical_train_vs_val 2 , n_identical_train_vs_test 3 
    ## Step 124 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 125 : n_identical_train_vs_val 3 , n_identical_train_vs_test 8 
    ## Step 126 : n_identical_train_vs_val 11 , n_identical_train_vs_test 2 
    ## Step 127 : n_identical_train_vs_val 7 , n_identical_train_vs_test 3 
    ## Step 128 : n_identical_train_vs_val 3 , n_identical_train_vs_test 6 
    ## Step 129 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 130 : n_identical_train_vs_val 2 , n_identical_train_vs_test 5 
    ## Step 131 : n_identical_train_vs_val 6 , n_identical_train_vs_test 4 
    ## Step 132 : n_identical_train_vs_val 4 , n_identical_train_vs_test 8 
    ## Step 133 : n_identical_train_vs_val 2 , n_identical_train_vs_test 1 
    ## Step 134 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 135 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 136 : n_identical_train_vs_val 1 , n_identical_train_vs_test 7 
    ## Step 137 : n_identical_train_vs_val 4 , n_identical_train_vs_test 8 
    ## Step 138 : n_identical_train_vs_val 4 , n_identical_train_vs_test 1 
    ## Step 139 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 140 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 141 : n_identical_train_vs_val 2 , n_identical_train_vs_test 5 
    ## Step 142 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 143 : n_identical_train_vs_val 6 , n_identical_train_vs_test 7 
    ## Step 144 : n_identical_train_vs_val 3 , n_identical_train_vs_test 3 
    ## Step 145 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 146 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 147 : n_identical_train_vs_val 5 , n_identical_train_vs_test 2 
    ## Step 148 : n_identical_train_vs_val 7 , n_identical_train_vs_test 6 
    ## Step 149 : n_identical_train_vs_val 5 , n_identical_train_vs_test 2 
    ## Step 150 : n_identical_train_vs_val 2 , n_identical_train_vs_test 3 
    ## Step 151 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 152 : n_identical_train_vs_val 1 , n_identical_train_vs_test 9 
    ## Step 153 : n_identical_train_vs_val 6 , n_identical_train_vs_test 5 
    ## Step 154 : n_identical_train_vs_val 6 , n_identical_train_vs_test 4 
    ## Step 155 : n_identical_train_vs_val 3 , n_identical_train_vs_test 4 
    ## Step 156 : n_identical_train_vs_val 4 , n_identical_train_vs_test 2 
    ## Step 157 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 158 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 159 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 160 : n_identical_train_vs_val 5 , n_identical_train_vs_test 3 
    ## Step 161 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 162 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 163 : n_identical_train_vs_val 7 , n_identical_train_vs_test 8 
    ## Step 164 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 165 : n_identical_train_vs_val 4 , n_identical_train_vs_test 2 
    ## Step 166 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 167 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 168 : n_identical_train_vs_val 3 , n_identical_train_vs_test 6 
    ## Step 169 : n_identical_train_vs_val 6 , n_identical_train_vs_test 1 
    ## Step 170 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 171 : n_identical_train_vs_val 3 , n_identical_train_vs_test 4 
    ## Step 172 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 173 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 174 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 175 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 176 : n_identical_train_vs_val 7 , n_identical_train_vs_test 3 
    ## Step 177 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 178 : n_identical_train_vs_val 8 , n_identical_train_vs_test 3 
    ## Step 179 : n_identical_train_vs_val 9 , n_identical_train_vs_test 3 
    ## Step 180 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 181 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 182 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 183 : n_identical_train_vs_val 7 , n_identical_train_vs_test 7 
    ## Step 184 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 185 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 186 : n_identical_train_vs_val 1 , n_identical_train_vs_test 7 
    ## Step 187 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 188 : n_identical_train_vs_val 7 , n_identical_train_vs_test 5 
    ## Step 189 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 190 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 191 : n_identical_train_vs_val 8 , n_identical_train_vs_test 3 
    ## Step 192 : n_identical_train_vs_val 6 , n_identical_train_vs_test 4 
    ## Step 193 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 194 : n_identical_train_vs_val 3 , n_identical_train_vs_test 6 
    ## Step 195 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 196 : n_identical_train_vs_val 7 , n_identical_train_vs_test 4 
    ## Step 197 : n_identical_train_vs_val 7 , n_identical_train_vs_test 3 
    ## Step 198 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 199 : n_identical_train_vs_val 8 , n_identical_train_vs_test 4 
    ## Step 200 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 201 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 202 : n_identical_train_vs_val 3 , n_identical_train_vs_test 3 
    ## Step 203 : n_identical_train_vs_val 8 , n_identical_train_vs_test 0 
    ## Step 204 : n_identical_train_vs_val 3 , n_identical_train_vs_test 4 
    ## Step 205 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 206 : n_identical_train_vs_val 6 , n_identical_train_vs_test 3 
    ## Step 207 : n_identical_train_vs_val 5 , n_identical_train_vs_test 6 
    ## Step 208 : n_identical_train_vs_val 6 , n_identical_train_vs_test 4 
    ## Step 209 : n_identical_train_vs_val 6 , n_identical_train_vs_test 4 
    ## Step 210 : n_identical_train_vs_val 7 , n_identical_train_vs_test 6 
    ## Step 211 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4 
    ## Step 212 : n_identical_train_vs_val 7 , n_identical_train_vs_test 6 
    ## Step 213 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 214 : n_identical_train_vs_val 5 , n_identical_train_vs_test 7 
    ## Step 215 : n_identical_train_vs_val 4 , n_identical_train_vs_test 1 
    ## Step 216 : n_identical_train_vs_val 4 , n_identical_train_vs_test 8 
    ## Step 217 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 218 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 219 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 220 : n_identical_train_vs_val 4 , n_identical_train_vs_test 4 
    ## Step 221 : n_identical_train_vs_val 2 , n_identical_train_vs_test 9 
    ## Step 222 : n_identical_train_vs_val 4 , n_identical_train_vs_test 3 
    ## Step 223 : n_identical_train_vs_val 7 , n_identical_train_vs_test 7 
    ## Step 224 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 225 : n_identical_train_vs_val 7 , n_identical_train_vs_test 5 
    ## Step 226 : n_identical_train_vs_val 7 , n_identical_train_vs_test 5 
    ## Step 227 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 228 : n_identical_train_vs_val 9 , n_identical_train_vs_test 7 
    ## Step 229 : n_identical_train_vs_val 5 , n_identical_train_vs_test 3 
    ## Step 230 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 231 : n_identical_train_vs_val 3 , n_identical_train_vs_test 5 
    ## Step 232 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 233 : n_identical_train_vs_val 8 , n_identical_train_vs_test 3 
    ## Step 234 : n_identical_train_vs_val 4 , n_identical_train_vs_test 6 
    ## Step 235 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 236 : n_identical_train_vs_val 3 , n_identical_train_vs_test 0 
    ## Step 237 : n_identical_train_vs_val 6 , n_identical_train_vs_test 5 
    ## Step 238 : n_identical_train_vs_val 2 , n_identical_train_vs_test 5 
    ## Step 239 : n_identical_train_vs_val 3 , n_identical_train_vs_test 7 
    ## Step 240 : n_identical_train_vs_val 2 , n_identical_train_vs_test 6 
    ## Step 241 : n_identical_train_vs_val 5 , n_identical_train_vs_test 5 
    ## Step 242 : n_identical_train_vs_val 4 , n_identical_train_vs_test 5 
    ## Step 243 : n_identical_train_vs_val 7 , n_identical_train_vs_test 3 
    ## Step 244 : n_identical_train_vs_val 1 , n_identical_train_vs_test 7 
    ## Step 245 : n_identical_train_vs_val 6 , n_identical_train_vs_test 2 
    ## Step 246 : n_identical_train_vs_val 2 , n_identical_train_vs_test 8 
    ## Step 247 : n_identical_train_vs_val 5 , n_identical_train_vs_test 3 
    ## Step 248 : n_identical_train_vs_val 2 , n_identical_train_vs_test 6 
    ## Step 249 : n_identical_train_vs_val 5 , n_identical_train_vs_test 4

### Oomycete

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

    ## Step 1 : n_identical_train_vs_val 0 , n_identical_train_vs_test 0 
    ## [1] "Good dataset split has been found without identical protein!"

``` r
fungi_final_split[[1]] %>% 
  data.table::fwrite("secreted_fungi_training.csv")

fungi_final_split[[2]] %>% 
  data.table::fwrite("secreted_fungi_validation.csv")

fungi_final_split[[3]] %>% 
  data.table::fwrite("secreted_fungi_testing.csv")
```

``` r
# Save the results

oomycete_final_split[[1]] %>% 
  data.table::fwrite("secreted_oomycete_training.csv")

oomycete_final_split[[2]] %>% 
  data.table::fwrite("secreted_oomycete_validation.csv")

oomycete_final_split[[3]] %>% 
  data.table::fwrite("secreted_oomycete_testing.csv")

```
