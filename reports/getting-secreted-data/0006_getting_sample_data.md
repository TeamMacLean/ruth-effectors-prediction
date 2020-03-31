Selecting the secreted data randomly from the SignalP-Prediction results
========================================================================

In this report, a process to get the secreted data randomly will be
shown. First,

Function
--------

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
sample_from_lookup_table <- function(data, data_var, lookup, lookup_var) {
  # Create lookup pattern for base_name
  lookup_pattern <- lookup %>%
    # Arrange by level of "specificity"
    select(patterns = {{ lookup_var }}) %>%
    mutate(
      specificity_level = stringr::str_count(patterns, "_")
    ) %>%
    arrange(desc(specificity_level)) %>%
    # Get vector
    pull(patterns) %>%
    as.character() %>%
    # Create regex pattern
    paste0("^", ., collapse = "|")

  merged_table <- dplyr::left_join(
    # Data with new base_names column
    data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        base_names =
          stringr::str_match(
            string = {{ data_var }},
            pattern = lookup_pattern
          ) %>%
            .[1] %>%
            as.character()
      ),
    # Lookup table with base_names
    lookup %>%
      dplyr::mutate_at(
        .vars = dplyr::vars({{ lookup_var }}),
        .funs = as.character
      ) %>%
      dplyr::select(base_names = {{ lookup_var }}, total_count),
    by = "base_names"
  ) %>%
    # Take samples
    dplyr::filter(!is.na(base_names)) %>%
    dplyr::group_by(base_names) %>%
    dplyr::mutate(
      available_count = n(),
      enough_data = available_count >= total_count
    ) %>%
    dplyr::select(-available_count)

  # Take samples without replacement for available data
  data_with_enough <- merged_table %>%
    dplyr::filter(enough_data) %>%
    dplyr::group_by(base_names) %>%
    dplyr::sample_n(
      size = total_count,
      replace = FALSE
    )

  # Take samples with replacement when missing data
  data_without_enough <- merged_table %>%
    dplyr::filter(!enough_data) %>%
    dplyr::group_by(base_names) %>%
    dplyr::sample_n(
      size = total_count,
      replace = TRUE
    )

  # Merge again
  merged_table <- rbind(
    data_with_enough,
    data_without_enough
  ) %>%
    # Remove intermediary base_names column
    dplyr::ungroup() %>%
    dplyr::mutate(oversampled = !enough_data) %>%
    dplyr::select(
      # -base_names,
      -enough_data
    )

  return(merged_table)
}
```

Bacteria
--------

``` r
# bacteria_lookup_table <- data.table::fread("../../../data/secreted_data/dataset-download-scripts/pathogen_species_final.csv") %>%
#   `colnames<-`(c("Pathogenspecies", "Count")) %>%
#   group_by(Pathogenspecies) %>%
#   summarise(total_count = sum(Count))
#
# # Read the data file results
# bacteria_full_table <-
#   data.table::fread("../../../data/secreted_data/signalp-pipeline/bacteria_full_table.csv") %>%
#   # Filter only data with signal peptide and with prediction > 0.9
#   dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
# Read the full table pf the prediction result of the SignalP

bacteria_full_table <- data.table::fread("../../../data/secreted_data/updated_signalp_results/bacteria_full_table.csv") %>%
  # Make the name of organisms consistent with the lookup table
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("^_", "")  %>%
      stringr::str_to_lower()
  ) %>% 
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)


bacteria_lookup_table <- data.table::fread("../../../data/secreted_data/dataset-download-scripts/pathogen_species_final.csv") %>%
  `colnames<-`(c("Pathogenspecies", "Count")) %>%
  group_by(Pathogenspecies) %>%
  summarise(total_count = sum(Count))

bacteria_lookup_table <- bacteria_lookup_table %>%
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_to_lower() %>%
      stringr::str_remove_all("\\.") %>%
      stringr::str_replace_all(" ", "_")
  )

bacteria_lookup_table <- bacteria_lookup_table %>%
  ungroup() %>% 
  # Manual name fixes
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_replace_all("\\(", "") %>%
      stringr::str_replace_all("\\)", "") %>%
      stringr::str_replace_all("-", "_"),
    Pathogenspecies = ifelse(
      Pathogenspecies == "xanthomonas_citri_subsp_malvacearum",
      "xanthomonas_citri_pv_malvacearum",
      Pathogenspecies
    ),
    Pathogenspecies = ifelse(
      Pathogenspecies == "pantoea_stewartii_subsp_stewartii",
      "pantoea_stewartii",
      Pathogenspecies
    )
  ) %>%
  # Fix repeated Pathogenspecies
  group_by(Pathogenspecies) %>%
  summarise(total_count = sum(total_count))
```

``` r
bacteria_non_eff_secreted_data  <- sample_from_lookup_table(
  data = bacteria_full_table,
  data_var = organism_name,
  lookup = bacteria_lookup_table,
  lookup_var = Pathogenspecies
)
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

``` r
bacteria_non_eff_secreted_data  %>% 
  saveRDS("../../../data/secreted_data/updated_signalp_results/bacteria_non_eff_secreted_data.RDS")
```

``` r
bacteria_lookup_table %>%
  knitr::kable()
```

| Pathogenspecies                                              |  total\_count|
|:-------------------------------------------------------------|-------------:|
| clavibacter\_michiganensis\_subsp\_michiganensis\_ncppb\_382 |             2|
| cystobacter\_fuscus\_dsm\_2262                               |             1|
| erwinia\_amylovora\_cfbp1430                                 |            10|
| pantoea\_stewartii                                           |             1|
| pantoea\_stewartii\_subsp\_stewartii\_dc283                  |             1|
| pseudomonas\_amygdali\_pv\_morsprunorum                      |             1|
| pseudomonas\_amygdali\_pv\_tabaci                            |             1|
| pseudomonas\_cichorii\_jbc1                                  |             1|
| pseudomonas\_savastanoi\_pv\_phaseolicola                    |             5|
| pseudomonas\_savastanoi\_pv\_phaseolicola\_1448a             |             6|
| pseudomonas\_savastanoi\_pv\_savastanoi\_ncppb\_3335         |             1|
| pseudomonas\_syringae                                        |             5|
| pseudomonas\_syringae\_pv\_maculicola                        |             2|
| pseudomonas\_syringae\_pv\_maculicola\_str\_es4326           |             2|
| pseudomonas\_syringae\_pv\_maculicola\_str\_m6               |             1|
| pseudomonas\_syringae\_pv\_pisi\_str\_1704b                  |             2|
| pseudomonas\_syringae\_pv\_spinaceae                         |             1|
| pseudomonas\_syringae\_pv\_syringae                          |             2|
| pseudomonas\_syringae\_pv\_syringae\_b728a                   |            10|
| pseudomonas\_syringae\_pv\_syringae\_gca\_001401075          |             1|
| pseudomonas\_syringae\_pv\_tomato                            |             4|
| pseudomonas\_syringae\_pv\_tomato\_str\_dc3000               |            29|
| ralstonia\_solanacearum                                      |             5|
| ralstonia\_solanacearum\_gmi1000                             |            40|
| ralstonia\_solanacearum\_psi07                               |             1|
| xanthomonas\_axonopodis\_pv\_citri\_str\_306                 |             7|
| xanthomonas\_axonopodis\_pv\_manihotis                       |             1|
| xanthomonas\_campestris\_pv\_campestris\_str\_8004           |             3|
| xanthomonas\_campestris\_pv\_campestris\_str\_atcc\_33913    |             3|
| xanthomonas\_campestris\_pv\_vesicatoria\_str\_85\_10        |            13|
| xanthomonas\_citri\_pv\_malvacearum                          |             1|
| xanthomonas\_oryzae\_pv\_oryzae                              |             7|
| xanthomonas\_oryzae\_pv\_oryzae\_pxo99                       |             1|
| xanthomonas\_oryzae\_pv\_oryzicola                           |             1|
| xanthomonas\_oryzae\_pv\_oryzicola\_bls256                   |            15|
| xylella\_fastidiosa\_temecula1                               |             3|

Getting info from the sample data we got

``` r
bacteria_sample_data_count <- bacteria_non_eff_secreted_data  %>%
  group_by(organism_name) %>%
  summarise(count = n())

bacteria_sample_data_count
```

    ## # A tibble: 36 x 2
    ##    organism_name                                                     count
    ##    <chr>                                                             <int>
    ##  1 clavibacter_michiganensis_subsp_michiganensis_ncppb_382.asm6348v1     2
    ##  2 cystobacter_fuscus_dsm_2262.asm33547v2                                1
    ##  3 erwinia_amylovora_cfbp1430.asm9156v1                                 10
    ##  4 pantoea_stewartii_subsp_stewartii_dc283.asm24839v2                    1
    ##  5 pantoea_stewartii.asm78625v1                                          1
    ##  6 pseudomonas_amygdali_pv_morsprunorum.pmpftrs_u7805                    1
    ##  7 pseudomonas_amygdali_pv_tabaci.asm93464v1                             1
    ##  8 pseudomonas_cichorii_jbc1.asm51730v1                                  1
    ##  9 pseudomonas_savastanoi_pv_phaseolicola_1448a.asm1220v1                6
    ## 10 pseudomonas_savastanoi_pv_phaseolicola.pphy5_2                        5
    ## # … with 26 more rows

### Check the debug on the sampling data

We can check here if there is missing data when we sampling the data
randomly. This can due to:

-   There is no enough predicted sequence with signalp
-   Or this can due to the

``` r
bacteria_matched_organisms_name <- left_join(
  bacteria_lookup_table %>%
    dplyr::mutate_at(
      .vars = dplyr::vars(Pathogenspecies),
      .funs = as.character
    ) %>%
    dplyr::select(base_names = Pathogenspecies, total_count_src = total_count),
  bacteria_non_eff_secreted_data,
  by = "base_names"
)

# Check if there is NA pathogen organims that can not be mapped
bacteria_matched_organisms_name  %>%
  group_by(base_names) %>%
  slice(1) %>%
  filter(is.na(organism_name)) %>%
  select(base_names, total_count_src, organism_name) %>% 
  knitr::kable()
```

| base\_names |  total\_count\_src| organism\_name |
|:------------|------------------:|:---------------|

From the results above, as we can see, there is no pathogen organism
with missing data. It means that organisms names are mapped.

Fungi
-----

``` r
fungi_lookup_table <- data.table::fread("../../../data/secreted_data/signalp-pipeline/fungi_lookup_table.csv", drop = "V1") %>%
  `colnames<-`(c("Pathogenspecies", "total_count"))

# Lookup table
fungi_lookup_table <- fungi_lookup_table %>%
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_to_lower() %>%
      # stringr::str_remove_all("\\.") %>%
      stringr::str_replace_all(" ", "_"), 
    Pathogenspecies = ifelse(
      Pathogenspecies == "pyrenophora_tritici-repentis",
      "pyrenophora_triticirepentis",
      Pathogenspecies
      )
  )

# Filter for the 
# not_available_fungi <- c("cercospora_apii", "cercospora_beticola", "pyrenophora_tritici-repentis", "melampsora_lini", "parastagonospora_nodorum", "passalora_fulva", "pseudocercospora_fuligena")

# fungi_lookup_table <- fungi_lookup_table %>% 
#   dplyr::filter(!(Pathogenspecies %in% not_available_fungi))

# Read the data file of signalP prediction results
fungi_full_table <-
  data.table::fread("../../../data/secreted_data/updated_signalp_results/fungi_full_table.csv")  %>%
  # Filter only data with signal peptide and with prediction > 0.9
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
# Getting all of the organism name consistent

# Results table
fungi_full_table <- fungi_full_table %>%
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("\\.", "_") %>%
      stringr::str_to_lower()
  )
```

``` r
fungi_non_eff_secreted_data <- sample_from_lookup_table(
  data = fungi_full_table,
  data_var = organism_name,
  lookup = fungi_lookup_table,
  lookup_var = Pathogenspecies
)
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

``` r
fungi_non_eff_secreted_data %>% 
  saveRDS("../../../data/secreted_data/updated_signalp_results/fungi_non_eff_secreted_data.RDS")
```

### Checking the debug

``` r
fungi_non_eff_secreted_count <- fungi_non_eff_secreted_data %>% 
  group_by(organism_name) %>% 
  summarise(count = n())

fungi_non_eff_secreted_count %>% 
  dplyr::select(count) %>% 
  sum()
```

    ## [1] 97

``` r
# Compare the lookup table and the result of sampling
data.frame(table = c("result", "lookup"), 
           count_seq = c(fungi_non_eff_secreted_count %>% select(count) %>% 
                           sum(), 
                         fungi_lookup_table %>% select(total_count) %>% 
                           sum()))
```

    ##    table count_seq
    ## 1 result        97
    ## 2 lookup       113

``` r
fungi_debug <- left_join(
  fungi_lookup_table %>%
    dplyr::mutate_at(
      .vars = dplyr::vars(Pathogenspecies),
      .funs = as.character
    ) %>%
    dplyr::select(base_names = Pathogenspecies, total_count_src = total_count),
  fungi_non_eff_secreted_data,
  by = "base_names"
)

fungi_debug %>% 
  group_by(base_names) %>% 
  slice(1) %>% 
  filter(is.na(organism_name)) %>%
  select(base_names, total_count_src, organism_name)
```

    ## # A tibble: 5 x 3
    ## # Groups:   base_names [5]
    ##   base_names                total_count_src organism_name
    ##   <chr>                               <int> <chr>        
    ## 1 cercospora_apii                         1 <NA>         
    ## 2 melampsora_lini                         5 <NA>         
    ## 3 parastagonospora_nodorum                2 <NA>         
    ## 4 passalora_fulva                         7 <NA>         
    ## 5 pseudocercospora_fuligena               1 <NA>

``` r
not_available_fungi <- fungi_debug %>% 
  group_by(organism_name, base_names) %>% 
  summarise(count_organism = n()) %>% 
  filter(is.na(organism_name)) 

not_available_fungi_list <- not_available_fungi %>% 
  pull(base_names)

not_available_fungi %>% 
  knitr::kable()
```

| organism\_name | base\_names                |  count\_organism|
|:---------------|:---------------------------|----------------:|
| NA             | cercospora\_apii           |                1|
| NA             | melampsora\_lini           |                1|
| NA             | parastagonospora\_nodorum  |                1|
| NA             | passalora\_fulva           |                1|
| NA             | pseudocercospora\_fuligena |                1|

Oomycete
--------

``` r
oomycete_lookup_table <- data.table::fread("../../../data/secreted_data/signalp-pipeline/oomycete_lookup_table.csv") %>% `colnames<-`(c("Pathogenspecies", "total_count"))


# Read the data file results
oomycete_full_table <-
  data.table::fread("../../../data/secreted_data/updated_signalp_results/protist_full_table.csv") %>%
  # Filter only data with signal peptide and with prediction > 0.9
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
# Getting all of the organism name consistent

# Lookup table
oomycete_lookup_table <- oomycete_lookup_table %>%
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_to_lower() %>%
      # stringr::str_remove_all("\\.") %>%
      stringr::str_replace_all(" ", "_")
  )

# Results table
oomycete_full_table <- oomycete_full_table %>%
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("\\.", "_") %>%
      stringr::str_to_lower()
  )
```

``` r
oomycete_non_eff_secreted_data <- sample_from_lookup_table(
  data = oomycete_full_table,
  data_var = organism_name,
  lookup = oomycete_lookup_table,
  lookup_var = Pathogenspecies
)
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

``` r
oomycete_non_eff_secreted_data %>% 
  saveRDS("../../../data/secreted_data/updated_signalp_results/oomycete_non_eff_secreted_data.RDS")
```

``` r
oomycete_debug <- left_join(
  oomycete_lookup_table %>%
    dplyr::mutate_at(
      .vars = dplyr::vars(Pathogenspecies),
      .funs = as.character
    ) %>%
    dplyr::select(base_names = Pathogenspecies, total_count_src = total_count),
  oomycete_non_eff_secreted_data,
  by = "base_names"
)

not_available_oomycete_list <- 
  oomycete_debug %>% 
  group_by(base_names) %>% 
  slice(1) %>% 
  filter(is.na(organism_name)) %>%
  select(base_names)
```

Summary
-------

Checking how many sequence data each pathogen has
-------------------------------------------------

``` r
effector_data_info <- readRDS("../../../data/secreted_data/data_processed_after_signalp/effector_data_info.RDS")

effector_data_info %>%
  group_by(class) %>%
  summarise(count_class = n())
```

    ## # A tibble: 5 x 2
    ##   class    count_class
    ##   <fct>          <int>
    ## 1 bacteria         190
    ## 2 fungi            113
    ## 3 insecta            2
    ## 4 nematoda           2
    ## 5 oomycete          95

After getting the info above, now needs to check the sum of total
sequence in lookup table, to make sure that the total is consistent.

``` r
data.frame(
  pathogen = c("bacteria", "fungi", "oomycete"),
  num_seq = c(
    bacteria_lookup_table %>% dplyr::select(total_count) %>% sum(),
    fungi_lookup_table %>% dplyr::select(total_count) %>% sum(),
    oomycete_lookup_table %>% dplyr::select(total_count) %>% sum()
  )
)
```

    ##   pathogen num_seq
    ## 1 bacteria     190
    ## 2    fungi     113
    ## 3 oomycete      95

#### Bacteria

``` r
num_seq_retrieved <- bacteria_non_eff_secreted_data %>% 
  nrow()

num_seq_needed <- bacteria_lookup_table %>% 
  select(total_count) %>% 
  sum()

# Summary of all the total count of the sequences needed, not available, and successfully needed
data.frame(desc = c("retrieved_seq", "lookup_table_seq"), 
              count = c(num_seq_retrieved, num_seq_needed)) %>% 
  knitr::kable()
```

| desc               |  count|
|:-------------------|------:|
| retrieved\_seq     |    190|
| lookup\_table\_seq |    190|

#### Fungi

List of fungi that does not exist

``` r
num_seq_retrieved <- fungi_non_eff_secreted_count %>% 
  select(count) %>% 
  sum()

num_seq_needed <- fungi_lookup_table %>% 
  select(total_count) %>% 
  sum()

num_seq_na <- fungi_lookup_table %>%
  dplyr::filter(Pathogenspecies %in% not_available_fungi_list) %>% 
  select(total_count) %>% 
  sum()

# Summary of all the total count of the sequences needed, not available, and successfully needed
data.frame(desc = c("retrieved_seq", "lookup_table_seq", "na_seq"), 
              count = c(num_seq_retrieved, num_seq_needed, num_seq_na)) %>% 
  knitr::kable()
```

| desc               |  count|
|:-------------------|------:|
| retrieved\_seq     |     97|
| lookup\_table\_seq |    113|
| na\_seq            |     16|

#### Oomycete

``` r
num_seq_retrieved <- oomycete_non_eff_secreted_data %>% 
  nrow()

num_seq_needed <- oomycete_lookup_table %>% 
  select(total_count) %>% 
  sum()

num_seq_na <- oomycete_lookup_table %>%
  dplyr::filter(Pathogenspecies %in% not_available_oomycete_list) %>% 
  select(total_count) %>% 
  sum()

data.frame(desc = c("retrieved_seq", "lookup_table_seq", "na_seq"), 
           count = c(num_seq_retrieved, num_seq_needed, num_seq_na)) %>% 
  knitr::kable()
```

| desc               |  count|
|:-------------------|------:|
| retrieved\_seq     |     85|
| lookup\_table\_seq |     95|
| na\_seq            |     10|
