Analysis of the Effector and Noneffector Sequence Data
======================================================

Introduction
------------

Before getting to know further whether the protein are identical or not
(we can determine by using BLAST), first of all, we can do simple
stastical analysis on the sequences themselves.

Distribution of the length of sequence
--------------------------------------

``` r
# Read CSV of the effector data
effector_data <- data.table::fread("effector_data.csv")

# Read CSV of the noneffector data 
non_effector_data <- data.table::fread("non_effector_data.csv")
```

``` r
effector_seq <- effector_data %>% 
  dplyr::select(Sequence) %>% 
  rename(., sequence = Sequence) %>% 
  mutate(data = 1)

non_effector_seq <- non_effector_data %>% 
  dplyr::select(sequence) %>% 
  mutate(data = 0)

seq_all <- effector_seq %>% 
  rbind(non_effector_seq)
```

``` r
# Count the length of the string

library(stringr)
seq_info <- seq_all %>% 
  rowwise() %>% 
  dplyr::mutate(len = stringr::str_length(sequence))
```

``` r
seq_info %>% 
  summary(len)
```

    ##    sequence              data             len        
    ##  Length:800         Min.   :0.0000   Min.   :  32.0  
    ##  Class :character   1st Qu.:0.0000   1st Qu.: 174.5  
    ##  Mode  :character   Median :1.0000   Median : 311.0  
    ##                     Mean   :0.5025   Mean   : 402.3  
    ##                     3rd Qu.:1.0000   3rd Qu.: 483.5  
    ##                     Max.   :1.0000   Max.   :4034.0

Plot using ggplot()

``` r
ggplot(seq_info) +
  geom_histogram(aes(len, fill = as.factor(data)), bins = 50) +
  facet_wrap(~data) +
  labs(x = "Length of Sequence", y = "Count", fill = "Data")
```

![](/Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/reports/0002-simple-analysis-prot-sequence_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
# Show how many data have length for more than 2000 chars

seq_info %>% 
  dplyr::filter(len > 2000) %>% 
  # Show only the beginning of sequence
  mutate(sequence = substr(sequence, 1, 30)) %>% 
  knitr::kable()
```

| sequence                       |  data|   len|
|:-------------------------------|-----:|-----:|
| MRDEMWNTATEPIAIIGSGCKFPGGSTTPS |     1|  4034|
| MKFNRTHRNLPPAHVGHTTRGHAPAPAGPQ |     1|  2338|
| MPSRMGYSRISSGLNASRGASPAPQPDTPP |     1|  2574|
| MTLFNGSNGANGTSSGHGAHPSANGFHNAA |     0|  2410|
| MFLKQQQQQPVEADFVDCQSTPRAGGPVAG |     0|  2758|

Analyzing Amino Acids in Protein Sequence data
----------------------------------------------

``` r
seq_with_all_aa_counts <- seq_info %>% 
  rowwise() %>% 
  mutate(
  G_count = stringr::str_count(sequence, "G"),
  A_count = stringr::str_count(sequence, "A"),
  L_count = stringr::str_count(sequence, "L"),
  M_count = stringr::str_count(sequence, "M"),
  F_count = stringr::str_count(sequence, "F"),
  W_count = stringr::str_count(sequence, "W"),
  K_count = stringr::str_count(sequence, "K"),
  Q_count = stringr::str_count(sequence, "Q"),
  E_count = stringr::str_count(sequence, "E"),
  S_count = stringr::str_count(sequence, "S"),
  P_count = stringr::str_count(sequence, "P"),
  V_count = stringr::str_count(sequence, "V"),
  I_count = stringr::str_count(sequence, "I"),
  C_count = stringr::str_count(sequence, "C"),
  Y_count = stringr::str_count(sequence, "Y"),
  H_count = stringr::str_count(sequence, "H"),
  R_count = stringr::str_count(sequence, "R"),
  N_count = stringr::str_count(sequence, "N"),
  D_count = stringr::str_count(sequence, "D"),
  T_count = stringr::str_count(sequence, "T")
)
```

``` r
test <- seq_with_all_aa_counts %>% 
  select(ends_with("_count")) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(
    sum_each = rowSums(.),
    mean = rowMeans(.)
  ) %>% 
  `rownames<-`(seq_with_all_aa_counts %>% names() %>% grep("_count", ., value = TRUE)) %>% 
  select(sum_each, mean)

total_amino <- test %>% 
   tibble::rownames_to_column() %>% 
   summarise(sum(sum_each)) %>% 
   unlist()

test_summary <- test %>% 
  tibble::rownames_to_column(var = "aminoacid") %>% 
  mutate(percent = sum_each / total_amino * 100) %>% 
  arrange(desc(percent)) %>% 
  mutate(aminoacid = stringr::str_remove_all(aminoacid, "_count"))
# saveRDS(test, "test.RDS")

test_summary %>% 
  knitr::kable()
```

| aminoacid |  sum\_each|      mean|    percent|
|:----------|----------:|---------:|----------:|
| A         |      35667|  44.58375|  11.083661|
| L         |      31156|  38.94500|   9.681850|
| G         |      23807|  29.75875|   7.398119|
| S         |      22804|  28.50500|   7.086433|
| V         |      21216|  26.52000|   6.592956|
| R         |      20368|  25.46000|   6.329436|
| D         |      18086|  22.60750|   5.620296|
| T         |      18032|  22.54000|   5.603515|
| E         |      17903|  22.37875|   5.563428|
| P         |      17662|  22.07750|   5.488536|
| Q         |      14565|  18.20625|   4.526131|
| I         |      14125|  17.65625|   4.389400|
| K         |      13985|  17.48125|   4.345894|
| N         |      11492|  14.36500|   3.571184|
| F         |      10576|  13.22000|   3.286534|
| H         |       8123|  10.15375|   2.524254|
| Y         |       7541|   9.42625|   2.343396|
| M         |       6930|   8.66250|   2.153525|
| W         |       3906|   4.88250|   1.213805|
| C         |       3854|   4.81750|   1.197646|

``` r
library(wesanderson)
ggplot(test_summary) +
  geom_col(aes(x = reorder(aminoacid, desc(percent)), y = percent), fill = "darksalmon") +
  scale_y_continuous(breaks = seq(0, 12, by = 1)) +
  labs(x = "Amino Acids", y = "Percentage")
```

![](/Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/reports/0002-simple-analysis-prot-sequence_files/figure-markdown_github/unnamed-chunk-10-1.png)
