Predicting the presence of signal peptide in fatsa data using SignalP
=====================================================================

Aim:
----

In order to have information whether a protein is secreted or not, we
can check the presence of signal peptide in the N-terminal sequence.
Using SignalP, we can do so.

Execution
---------

``` r
get_signalp_pred <- function(fasta_filename, verbose = FALSE) {
  signalp_path <- "/Users/kristian/Documents/Workspace/Software/signalp/bin"
  verbose_string <- tolower(deparse(substitute(verbose)))
  fasta_path <- here::here("scripts/r-scripts/getting-data-new/multi-class", fasta_filename)

  command_string <- paste0(
    "export PATH=", signalp_path, ":$PATH ;",
    " signalp -fasta ", fasta_path, " -org euk -format short ",
    "-verbose=", verbose_string, " -stdout"
  )

  data_lines <- command_string %>%
    system(intern = TRUE)

  data <- data_lines %>%
    .[-c(1, 2)] %>%
    as.list() %>%
    data.frame() %>%
    t() %>%
    as_tibble() %>%
    tidyr::separate(col = V1, into = c("ID", "Prediction", "SP(Sec/SPI)", "OTHER", "CS Position"), sep = "\t")

  return(data)
}
```

``` r
bacteria_test_pred <- get_signalp_pred("bacteria_test_non_eff.fasta")
```

    ## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
    ## This warning is displayed once per session.

``` r
bacteria_train_pred <- get_signalp_pred("bacteria_train_non_eff.fasta")
bacteria_val_pred <- get_signalp_pred("bacteria_val_non_eff.fasta")

fungi_test_pred <- get_signalp_pred("fungi_test_non_eff.fasta")
fungi_train_pred <- get_signalp_pred("fungi_train_non_eff.fasta")
fungi_val_pred <- get_signalp_pred("fungi_val_non_eff.fasta")

oomycete_test_pred <- get_signalp_pred("oomycete_test_non_eff.fasta")
oomycete_train_pred <- get_signalp_pred("oomycete_train_non_eff.fasta")
oomycete_val_pred <- get_signalp_pred("oomycete_val_non_eff.fasta")
```

``` r
# Bacteria 

summarise_count_sp <- function(data){
  sum_data <- data %>% 
  group_by(Prediction) %>% 
  summarise(count = n())
  # %>% 
  # tibble::rownames_to_column() %>%
  # .[,-1] %>% 
  # pivot_wider(names_from = "Prediction", values_from = "count")
 
  
  return(sum_data)
}
```

``` r
# Get the bacteria summary
bacteria_test_summary <- summarise_count_sp(bacteria_test_pred)
bacteria_train_summary <- summarise_count_sp(bacteria_train_pred)
bacteria_val_summary <- summarise_count_sp(bacteria_val_pred)

# Get the fungi summary
fungi_test_summary <- summarise_count_sp(fungi_test_pred)
fungi_val_summary <- summarise_count_sp(fungi_val_pred)
fungi_train_summary <- summarise_count_sp(fungi_train_pred)

# Get the oomycete summary
oomycete_test_summary <- summarise_count_sp(oomycete_test_pred)
oomycete_val_summary <- summarise_count_sp(oomycete_val_pred)
oomycete_train_summary <- summarise_count_sp(oomycete_train_pred)
```

``` r
# Combine all of the data based on the pathogen 
bacteria_all <- rbind(
  cbind(bacteria_train_summary, dataset = "train", pathogen = "bacteria"), 
  cbind(bacteria_val_summary, dataset = "val", pathogen = "bacteria"),  
  cbind(bacteria_test_summary, dataset = "test", pathogen = "bacteria")
)

fungi_all <- rbind(
  cbind(fungi_train_summary, dataset = "train", pathogen = "fungi"), 
  cbind(fungi_val_summary, dataset = "val", pathogen = "fungi"),  
  cbind(fungi_test_summary, dataset = "test", pathogen = "fungi")
)

oomycete_all <- rbind(
  cbind(oomycete_train_summary, dataset = "train", pathogen = "oomycete"), 
  cbind(oomycete_val_summary, dataset = "val", pathogen = "oomycete"),  
  cbind(oomycete_test_summary, dataset = "test", pathogen = "oomycete")
)

non_eff_all <- bacteria_all %>% rbind(., fungi_all, oomycete_all)

non_eff_all %>% 
  pivot_wider(id_cols = c("pathogen", "dataset"), names_from = "Prediction", values_from = "count")
```

    ## # A tibble: 9 x 4
    ##   pathogen dataset OTHER `SP(Sec/SPI)`
    ##   <fct>    <fct>   <int>         <int>
    ## 1 bacteria train     105             9
    ## 2 bacteria val        38            NA
    ## 3 bacteria test       31             7
    ## 4 fungi    train      55            17
    ## 5 fungi    val        13             4
    ## 6 fungi    test       17             5
    ## 7 oomycete train      36            21
    ## 8 oomycete val        11             8
    ## 9 oomycete test       12             7

``` r
# All bacteria protein that have Signal peptide
rbind(bacteria_test_pred, bacteria_train_pred, bacteria_val_pred) %>% 
  dplyr::filter(Prediction == "SP(Sec/SPI)") %>% 
  knitr::kable()
```

| ID  | Prediction  | SP(Sec/SPI) | OTHER    | CS Position                       |
|:----|:------------|:------------|:---------|:----------------------------------|
| 2   | SP(Sec/SPI) | 0.923110    | 0.076890 | CS pos: 20-21. GFA-AT. Pr: 0.5289 |
| 7   | SP(Sec/SPI) | 0.972276    | 0.027724 | CS pos: 27-28. VHA-AV. Pr: 0.8946 |
| 12  | SP(Sec/SPI) | 0.526236    | 0.473764 | CS pos: 19-20. ADC-SQ. Pr: 0.4209 |
| 13  | SP(Sec/SPI) | 0.963586    | 0.036414 | CS pos: 29-30. VLA-NS. Pr: 0.7967 |
| 24  | SP(Sec/SPI) | 0.691623    | 0.308377 | CS pos: 24-25. SFS-IR. Pr: 0.2003 |
| 28  | SP(Sec/SPI) | 0.979053    | 0.020947 | CS pos: 24-25. VHA-AE. Pr: 0.8717 |
| 30  | SP(Sec/SPI) | 0.986228    | 0.013772 | CS pos: 18-19. SGC-SL. Pr: 0.3877 |
| 18  | SP(Sec/SPI) | 0.755491    | 0.244509 | CS pos: 30-31. AAA-NA. Pr: 0.3632 |
| 35  | SP(Sec/SPI) | 0.511312    | 0.488688 | CS pos: 44-45. ANA-AY. Pr: 0.3834 |
| 37  | SP(Sec/SPI) | 0.929118    | 0.070882 | CS pos: 30-31. CSS-HA. Pr: 0.6143 |
| 42  | SP(Sec/SPI) | 0.905464    | 0.094536 | CS pos: 23-24. AVA-EA. Pr: 0.7269 |
| 49  | SP(Sec/SPI) | 0.846630    | 0.153370 | CS pos: 20-21. ASA-CS. Pr: 0.2930 |
| 62  | SP(Sec/SPI) | 0.922936    | 0.077064 | CS pos: 15-16. ARP-DL. Pr: 0.5439 |
| 67  | SP(Sec/SPI) | 0.996715    | 0.003285 | CS pos: 22-23. GHA-SQ. Pr: 0.9490 |
| 87  | SP(Sec/SPI) | 0.749137    | 0.250863 | CS pos: 30-31. CAA-TG. Pr: 0.3115 |
| 109 | SP(Sec/SPI) | 0.616773    | 0.383227 | CS pos: 20-21. AAA-GG. Pr: 0.2716 |

``` r
# All protein that have Signal peptide
rbind(fungi_test_pred, fungi_train_pred, fungi_val_pred) %>% 
  dplyr::filter(Prediction == "SP(Sec/SPI)") %>% 
  knitr::kable()
```

| ID  | Prediction  | SP(Sec/SPI) | OTHER    | CS Position                       |
|:----|:------------|:------------|:---------|:----------------------------------|
| 3   | SP(Sec/SPI) | 0.995877    | 0.004123 | CS pos: 18-19. SFA-DL. Pr: 0.9681 |
| 4   | SP(Sec/SPI) | 0.998760    | 0.001240 | CS pos: 21-22. ATA-AK. Pr: 0.7254 |
| 18  | SP(Sec/SPI) | 0.932384    | 0.067616 | CS pos: 21-22. AAA-FG. Pr: 0.5978 |
| 19  | SP(Sec/SPI) | 0.990190    | 0.009810 | CS pos: 17-18. AFA-AV. Pr: 0.3604 |
| 20  | SP(Sec/SPI) | 0.995759    | 0.004241 | CS pos: 20-21. VLC-AG. Pr: 0.6922 |
| 3   | SP(Sec/SPI) | 0.995734    | 0.004266 | CS pos: 14-15. ALC-AP. Pr: 0.7187 |
| 6   | SP(Sec/SPI) | 0.972945    | 0.027055 | CS pos: 16-17. AVA-VP. Pr: 0.6771 |
| 7   | SP(Sec/SPI) | 0.698142    | 0.301858 | CS pos: 34-35. TLA-AD. Pr: 0.4687 |
| 11  | SP(Sec/SPI) | 0.979235    | 0.020765 | CS pos: 17-18. STA-TP. Pr: 0.7041 |
| 15  | SP(Sec/SPI) | 0.974366    | 0.025634 | CS pos: 21-22. VNA-ET. Pr: 0.9137 |
| 16  | SP(Sec/SPI) | 0.972880    | 0.027120 | CS pos: 21-22. VNA-ET. Pr: 0.9127 |
| 36  | SP(Sec/SPI) | 0.985289    | 0.014711 | CS pos: 19-20. VLG-RP. Pr: 0.8138 |
| 43  | SP(Sec/SPI) | 0.975586    | 0.024414 | CS pos: 21-22. VNA-ET. Pr: 0.9145 |
| 44  | SP(Sec/SPI) | 0.973824    | 0.026176 | CS pos: 21-22. VNA-ET. Pr: 0.9129 |
| 45  | SP(Sec/SPI) | 0.975072    | 0.024928 | CS pos: 21-22. VNA-ET. Pr: 0.9138 |
| 48  | SP(Sec/SPI) | 0.978301    | 0.021699 | CS pos: 18-19. AVA-VS. Pr: 0.6589 |
| 57  | SP(Sec/SPI) | 0.987746    | 0.012254 | CS pos: 16-17. ALA-VP. Pr: 0.8529 |
| 58  | SP(Sec/SPI) | 0.852741    | 0.147259 | CS pos: 20-21. VQA-MP. Pr: 0.8054 |
| 59  | SP(Sec/SPI) | 0.993669    | 0.006331 | CS pos: 20-21. TLA-VP. Pr: 0.8207 |
| 65  | SP(Sec/SPI) | 0.806833    | 0.193167 | CS pos: 36-37. IVS-LP. Pr: 0.4104 |
| 66  | SP(Sec/SPI) | 0.758296    | 0.241704 | CS pos: 36-37. IVS-LP. Pr: 0.3044 |
| 70  | SP(Sec/SPI) | 0.999686    | 0.000314 | CS pos: 16-17. ALA-SP. Pr: 0.9418 |
| 9   | SP(Sec/SPI) | 0.982845    | 0.017155 | CS pos: 18-19. ISA-QY. Pr: 0.9179 |
| 11  | SP(Sec/SPI) | 0.995013    | 0.004987 | CS pos: 24-25. AQG-EA. Pr: 0.8119 |
| 12  | SP(Sec/SPI) | 0.998621    | 0.001379 | CS pos: 24-25. AYG-EE. Pr: 0.9682 |
| 14  | SP(Sec/SPI) | 0.560519    | 0.439481 | CS pos: 36-37. IVS-LP. Pr: 0.1991 |

``` r
# All protein that have Signal peptide
rbind(oomycete_test_pred, oomycete_train_pred, oomycete_val_pred) %>% 
  dplyr::filter(Prediction == "SP(Sec/SPI)") %>% 
  knitr::kable()
```

| ID  | Prediction  | SP(Sec/SPI) | OTHER    | CS Position                       |
|:----|:------------|:------------|:---------|:----------------------------------|
| 4   | SP(Sec/SPI) | 0.995169    | 0.004831 | CS pos: 22-23. ADA-QY. Pr: 0.9601 |
| 7   | SP(Sec/SPI) | 0.999070    | 0.000930 | CS pos: 19-20. VRA-QD. Pr: 0.9727 |
| 8   | SP(Sec/SPI) | 0.992980    | 0.007020 | CS pos: 19-20. VDA-SL. Pr: 0.8396 |
| 9   | SP(Sec/SPI) | 0.995370    | 0.004630 | CS pos: 20-21. TDA-TQ. Pr: 0.8813 |
| 10  | SP(Sec/SPI) | 0.995885    | 0.004115 | CS pos: 21-22. VTS-QP. Pr: 0.8284 |
| 11  | SP(Sec/SPI) | 0.994922    | 0.005078 | CS pos: 21-22. VTS-QP. Pr: 0.8324 |
| 13  | SP(Sec/SPI) | 0.974768    | 0.025232 | CS pos: 21-22. VAG-AI. Pr: 0.5748 |
| 10  | SP(Sec/SPI) | 0.992178    | 0.007822 | CS pos: 21-22. CLA-AG. Pr: 0.3358 |
| 24  | SP(Sec/SPI) | 0.961206    | 0.038794 | CS pos: 22-23. VVS-ES. Pr: 0.3810 |
| 28  | SP(Sec/SPI) | 0.995271    | 0.004729 | CS pos: 19-20. THA-QE. Pr: 0.9508 |
| 29  | SP(Sec/SPI) | 0.913547    | 0.086453 | CS pos: 21-22. SSP-QE. Pr: 0.4142 |
| 30  | SP(Sec/SPI) | 0.996660    | 0.003340 | CS pos: 20-21. VKS-DQ. Pr: 0.9505 |
| 31  | SP(Sec/SPI) | 0.961910    | 0.038090 | CS pos: 20-21. VSS-PD. Pr: 0.3148 |
| 32  | SP(Sec/SPI) | 0.998757    | 0.001243 | CS pos: 16-17. ANA-ND. Pr: 0.9410 |
| 33  | SP(Sec/SPI) | 0.999520    | 0.000480 | CS pos: 21-22. VSA-QR. Pr: 0.8358 |
| 34  | SP(Sec/SPI) | 0.981863    | 0.018137 | CS pos: 20-21. ANA-EQ. Pr: 0.8289 |
| 35  | SP(Sec/SPI) | 0.998009    | 0.001991 | CS pos: 24-25. AHA-RS. Pr: 0.8121 |
| 36  | SP(Sec/SPI) | 0.988736    | 0.011264 | CS pos: 19-20. TTA-DT. Pr: 0.9186 |
| 37  | SP(Sec/SPI) | 0.999380    | 0.000620 | CS pos: 25-26. AEG-QS. Pr: 0.4484 |
| 38  | SP(Sec/SPI) | 0.991375    | 0.008625 | CS pos: 18-19. AAA-MS. Pr: 0.5337 |
| 39  | SP(Sec/SPI) | 0.995828    | 0.004172 | CS pos: 18-19. ANG-GT. Pr: 0.9346 |
| 40  | SP(Sec/SPI) | 0.987656    | 0.012344 | CS pos: 22-23. VVG-RL. Pr: 0.4803 |
| 41  | SP(Sec/SPI) | 0.996922    | 0.003078 | CS pos: 21-22. VTS-QP. Pr: 0.8477 |
| 42  | SP(Sec/SPI) | 0.996315    | 0.003685 | CS pos: 21-22. AMG-DG. Pr: 0.6901 |
| 43  | SP(Sec/SPI) | 0.997944    | 0.002056 | CS pos: 24-25. SSA-WD. Pr: 0.4712 |
| 44  | SP(Sec/SPI) | 0.996044    | 0.003956 | CS pos: 19-20. GSA-ME. Pr: 0.6258 |
| 45  | SP(Sec/SPI) | 0.995919    | 0.004081 | CS pos: 19-20. GSA-TE. Pr: 0.5516 |
| 47  | SP(Sec/SPI) | 0.996018    | 0.003982 | CS pos: 17-18. ASA-AP. Pr: 0.8388 |
| 6   | SP(Sec/SPI) | 0.984681    | 0.015319 | CS pos: 19-20. VQA-QE. Pr: 0.9236 |
| 7   | SP(Sec/SPI) | 0.897876    | 0.102124 | CS pos: 19-20. VRA-QE. Pr: 0.7945 |
| 8   | SP(Sec/SPI) | 0.996659    | 0.003341 | CS pos: 20-21. VKS-DQ. Pr: 0.9518 |
| 10  | SP(Sec/SPI) | 0.984648    | 0.015352 | CS pos: 20-21. AKS-FQ. Pr: 0.8413 |
| 11  | SP(Sec/SPI) | 0.992903    | 0.007097 | CS pos: 20-21. ARG-AN. Pr: 0.4894 |
| 12  | SP(Sec/SPI) | 0.994210    | 0.005790 | CS pos: 20-21. VNG-QG. Pr: 0.9058 |
| 14  | SP(Sec/SPI) | 0.973570    | 0.026430 | CS pos: 21-22. CIA-QS. Pr: 0.3089 |
| 18  | SP(Sec/SPI) | 0.970780    | 0.029220 | CS pos: 19-20. AIA-DP. Pr: 0.6493 |
