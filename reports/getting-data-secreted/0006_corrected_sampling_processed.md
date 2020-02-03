Re-do the resampling data with BLAST
====================================

Background
----------

Based on the innefective sampling methods I have done before, I was not
effective, since I did not do BLAST while sampling the data, this
results on the additional effort to do BLAST while splitting the data
which there was no point since if we have identical protein sequences on
our dataset, no matter how hard we blast to avoiud identical data
between datasets, it will still be identical protein sequences inside
(this is useless).

Approach
--------

A new way to tackle this issue (not to have identical protein data
inside the datasets) are doing BLAST while doing the sampling process.

Execution
---------

### Load libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(dplyr)

# Source functions
source(here::here("scripts/r-scripts/r-functions", "blast_data.R"))
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colMeans,
    ##     colnames, colSums, dirname, do.call, duplicated, eval, evalq,
    ##     Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax, pmax.int,
    ##     pmin, pmin.int, Position, rank, rbind, Reduce, rowMeans, rownames,
    ##     rowSums, sapply, setdiff, sort, table, tapply, union, unique,
    ##     unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'XVector'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
# source(here::here("scripts/r-scripts/r-functions", "dealing_with_df.R"))
source(here::here("scripts/r-scripts/r-functions", "getting_sample_data.R"))
```

### Run with bacteria

#### Load sequence and lookup data

``` r
# Read the full table of the prediction result of the SignalP
bacteria_full_table <- data.table::fread(here::here("data/secreted_data/updated_signalp_results", "bacteria_full_table.csv")) %>%
  # Make the name of organisms consistent with the lookup table
  dplyr::mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("^_", "") %>%
      stringr::str_to_lower()
  ) %>%
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
bacteria_lookup_table <- data.table::fread(here::here("data/secreted_data/signalp-pipeline", "bacteria_lookup_table.csv")) %>%
  `colnames<-`(c("Pathogenspecies", "Count")) %>%
  group_by(Pathogenspecies) %>%
  summarise(samples_needed = sum(Count)) %>%
  ungroup()

bacteria_lookup_table <- bacteria_lookup_table %>%
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_to_lower() %>%
      stringr::str_remove_all("\\.") %>%
      stringr::str_replace_all(" ", "_")
  )

bacteria_lookup_table <- bacteria_lookup_table %>%
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
  summarise(samples_needed = sum(samples_needed))
```

#### Merge sequence data with lookup

``` r
bacteria_merged_table <- merge_data_with_lookup_table(
  bag_data = bacteria_full_table,
  bag_var = organism_name,
  lookup_data = bacteria_lookup_table,
  lookup_var = Pathogenspecies
)

# Check if the merged table doesn't have any NAs
if (
  length(
    bacteria_merged_table %>%
      filter(samples_needed %>% is.na()) %>%
      pull(organism_name) %>% unique()
  ) == 0
) {
  rm(bacteria_full_table, bacteria_lookup_table)
}
```

#### Sample sequences

``` r
bacteria_sampled_table <- bacteria_merged_table %>%
  map_sampleing_function(
    col_organism = organism_name, 
    col_id = ID, 
    col_seq = sequence
  )
```

#### Check quality of samples

We check first if enough samples were obtained compares to the needed
samples from the lookup table.

``` r
# Read the results 

bacteria_sampled_table <- readRDS("../../../data/secreted_data/ready_to_process/sampled_data_without_identical/bacteria_sampled_table_good.RDS")
```

``` r
bacteria_sampled_table %>% 
  dplyr::group_by(organism_name) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(
    needed_samples = sum(num_samples),
    actual_samples = sum(sampled_number)
  )
```

    ## # A tibble: 1 x 2
    ##   needed_samples actual_samples
    ##            <int>          <int>
    ## 1            190            190

Then we check if there are any duplicate sequences within the samples.

``` r
bacteria_sampled_table %>% blast_with_ifself(
  col_id = ID,
  col_seq = sequence,
  percent_threshold = 90
) %>% 
  knitr::kable()
```

    ## The data frame has been saved in /var/folders/c2/s7xgl0x53ds16553lpzktdp5hb4zny/T//RtmpM9ugOc/fasta_self.fasta

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 19%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">qseqid</th>
<th style="text-align: right;">qlen</th>
<th style="text-align: left;">sseqid</th>
<th style="text-align: right;">slen</th>
<th style="text-align: right;">length</th>
<th style="text-align: right;">nident</th>
<th style="text-align: right;">mismatch</th>
<th style="text-align: right;">positive</th>
<th style="text-align: right;">percent_identical</th>
<th style="text-align: right;">percent_positive</th>
</tr>
</thead>
<tbody>
</tbody>
</table>

``` r
bacteria_sampled_table %>% blast_with_ifself(
  col_id = ID,
  col_seq = sequence,
  percent_threshold = 10
) %>% 
  knitr::kable()
```

    ## The data frame has been saved in /var/folders/c2/s7xgl0x53ds16553lpzktdp5hb4zny/T//RtmpM9ugOc/fasta_self.fasta

| qseqid   |  qlen| sseqid   |  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|:---------|-----:|:---------|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
| CBA20126 |   430| AAO57431 |   433|     435|     189|       238|       265|            43.64896|           61.20092|
| AAO57431 |   433| CBA20126 |   430|     435|     189|       238|       265|            43.64896|           61.20092|
| CBA20245 |   402| AAZ36322 |   395|     397|     166|       223|       248|            41.29353|           61.69154|
| AAZ36322 |   395| CBA20245 |   402|     370|     163|       199|       237|            40.54726|           58.95522|
| AAZ34209 |   321| KPY81876 |   324|     292|     129|       153|       183|            39.81481|           56.48148|
| KPY81876 |   324| AAZ34209 |   321|     292|     129|       153|       183|            39.81481|           56.48148|
| CAD17488 |   351| CAD16800 |   376|     387|     118|       193|       177|            31.38298|           47.07447|
| AAO56997 |   642| ALS93752 |   657|     633|     199|       414|       288|            30.28919|           43.83562|
| CAD16800 |   376| CAD17488 |   351|     374|     113|       186|       171|            30.05319|           45.47872|
| ALS93752 |   657| AAO56997 |   642|     620|     197|       402|       284|            29.98478|           43.22679|
| CAJ23883 |   252| AAO55380 |   270|     224|      80|       112|       111|            29.62963|           41.11111|
| CBA20264 |   243| KFE54762 |   250|     244|      74|       157|       120|            29.60000|           48.00000|
| KFE54762 |   250| CBA20264 |   243|     244|      74|       157|       120|            29.60000|           48.00000|
| AAO55380 |   270| CAJ23883 |   252|     196|      74|        92|        98|            27.40741|           36.29630|
| CBA19959 |   397| KPB79686 |   388|     342|     107|       210|       171|            26.95214|           43.07305|
| KPB79686 |   388| CBA19959 |   397|     340|     105|       214|       169|            26.44836|           42.56927|
| AAO56109 |   685| ALS94050 |   625|     692|     174|       401|       270|            25.40146|           39.41606|
| ALS94050 |   625| AAO56109 |   685|     692|     174|       401|       270|            25.40146|           39.41606|
| KPY70815 |   123| AEQ97890 |   138|     128|      34|        79|        62|            24.63768|           44.92754|
| AEQ97890 |   138| KPY70815 |   123|     128|      34|        79|        62|            24.63768|           44.92754|
| AAO58528 |   204| CAJ25915 |   215|     185|      48|       125|        86|            22.32558|           40.00000|
| CAD13699 |   246| CAJ23883 |   252|     193|      55|       124|        77|            21.82540|           30.55556|
| CAJ23883 |   252| CAD13699 |   246|     193|      54|       125|        77|            21.42857|           30.55556|
| CAJ25915 |   215| AAO58528 |   204|     181|      45|       129|        84|            20.93023|           39.06977|
| CAJ22591 |   469| ALS95173 |   441|     328|      98|       168|       144|            20.89552|           30.70362|
| AAO56109 |   685| CAD14502 |   705|     638|     147|       391|       243|            20.85106|           34.46809|
| CAD14502 |   705| AAO56109 |   685|     638|     147|       391|       243|            20.85106|           34.46809|
| ALS95173 |   441| CAJ22591 |   469|     328|      95|       171|       143|            20.25586|           30.49041|
| AAO57168 |   543| CAD17454 |   515|     359|     107|       222|       182|            19.70534|           33.51750|
| CAD17454 |   515| AAO57168 |   543|     359|     107|       222|       182|            19.70534|           33.51750|
| AAZ34192 |   177| KFE54371 |   141|     148|      34|        81|        57|            19.20904|           32.20339|
| KFE54371 |   141| AAZ34192 |   177|     148|      34|        81|        57|            19.20904|           32.20339|
| KPY81876 |   324| KFE54762 |   250|     225|      58|       131|        97|            17.90123|           29.93827|
| CBA19959 |   397| AAO54363 |   371|     238|      71|       138|       106|            17.88413|           26.70025|
| AAO54363 |   371| CBA19959 |   397|     238|      71|       138|       106|            17.88413|           26.70025|
| AAO55380 |   270| CAD13699 |   246|     187|      46|       126|        75|            17.03704|           27.77778|
| CAD13699 |   246| AAO55380 |   270|     187|      46|       126|        75|            17.03704|           27.77778|
| AAO58528 |   204| AAO56826 |   207|     154|      35|       104|        54|            16.90821|           26.08696|
| CAD16424 |   198| AAZ34192 |   177|     135|      32|        82|        52|            16.16162|           26.26263|
| CAD18216 |   187| AAO54866 |   190|      87|      26|        52|        34|            13.68421|           17.89474|
| CAN01605 |   258| ACD58249 |   180|     151|      35|        93|        60|            13.56589|           23.25581|
| ACD58249 |   180| CAN01605 |   258|     151|      34|        94|        58|            13.17829|           22.48062|
| AAY35834 |   314| AEQ94279 |   268|     164|      40|        96|        63|            12.73885|           20.06369|
| KPB79686 |   388| AAM39135 |   314|     135|      48|        75|        72|            12.37113|           18.55670|
| AAM39135 |   314| KPB79686 |   388|     135|      48|        75|        72|            12.37113|           18.55670|
| EGH45942 |   184| CBA20614 |   333|     112|      41|        66|        68|            12.31231|           20.42042|
| AAO54363 |   371| AAM39135 |   314|     138|      45|        89|        68|            12.12938|           18.32884|
| AAM39135 |   314| AAO54363 |   371|     138|      45|        89|        68|            12.12938|           18.32884|
| AAO54866 |   190| CAD18216 |   187|      75|      23|        45|        30|            12.10526|           15.78947|
| EPX56555 |   176| AAO56826 |   207|      87|      25|        45|        39|            12.07729|           18.84058|
| AAO56826 |   207| EPX56555 |   176|      87|      25|        45|        39|            12.07729|           18.84058|
| KPC04104 |    90| KGK93493 |   207|     102|      24|        50|        41|            11.59420|           19.80676|
| AEQ96119 |   111| AAY35709 |   188|      78|      21|        43|        33|            11.17021|           17.55319|
| EHT99880 |   272| AAM40940 |   332|     147|      37|        69|        52|            11.14458|           15.66265|
| AAM40940 |   332| EHT99880 |   272|     147|      37|        69|        52|            11.14458|           15.66265|
| CAD17319 |   200| ALS95397 |   198|      78|      22|        38|        29|            11.00000|           14.50000|
| CAJ22785 |   108| AAO57894 |   155|      54|      17|        27|        27|            10.96774|           17.41935|
| KPB37766 |   322| KFE54762 |   250|     129|      35|        77|        55|            10.86957|           17.08075|
| KPB79686 |   388| AAO54363 |   371|     163|      42|       105|        67|            10.82474|           17.26804|
| AAO54363 |   371| KPB79686 |   388|     163|      42|       105|        67|            10.82474|           17.26804|
| AEQ94816 |   125| AAO58251 |   342|     123|      36|        61|        46|            10.52632|           13.45029|
| KFE56648 |   301| AAY36061 |   229|      74|      31|        43|        52|            10.29900|           17.27575|
| AAY36061 |   229| KFE56648 |   301|      74|      31|        43|        52|            10.29900|           17.27575|
| AAO58528 |   204| AAO54667 |   193|      57|      21|        25|        28|            10.29412|           13.72549|
| EGH58888 |   102| CAJ22785 |   108|      40|      11|        29|        20|            10.18519|           18.51852|
| CAJ22785 |   108| EGH58888 |   102|      40|      11|        29|        20|            10.18519|           18.51852|

#### Save results

``` r
# bacteria_sampled_table %>%
#   saveRDS(here::here("scripts/r-scripts/getting-secreted-data", "bacteria_sampled_table_good02.RDS"))
```

### Run with fungi

#### Load sequence and lookup data

``` r
# Read the full table of the prediction result of the SignalP
fungi_full_table <- data.table::fread(here::here("data/secreted_data/updated_signalp_results", "fungi_full_table.csv")) %>%
  # Make the name of organisms consistent with the lookup table
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("^_", "") %>%
      stringr::str_to_lower()
  ) %>%
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
fungi_lookup_table <- data.table::fread(here::here("data/secreted_data/signalp-pipeline", "fungi_lookup_table.csv")) %>%
  # `colnames<-`(c("Pathogenspecies", "Count")) %>%
  group_by(Pathogenspecies) %>%
  summarise(samples_needed = sum(Count)) %>%
  ungroup()

fungi_lookup_table <- fungi_lookup_table %>%
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_to_lower() %>%
      stringr::str_remove_all("\\.") %>%
      stringr::str_replace_all(" ", "_")
  )

fungi_lookup_table <- fungi_lookup_table %>%
  # Manual name fixes
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_replace_all("\\(", "") %>%
      stringr::str_replace_all("\\)", "") %>%
      stringr::str_replace_all("-", "_"),
    Pathogenspecies = ifelse(
      Pathogenspecies == "pyrenophora_tritici_repentis",
      "pyrenophora_triticirepentis",
      Pathogenspecies
    ),
  ) %>%
  # Fix repeated Pathogenspecies
  group_by(Pathogenspecies) %>%
  summarise(samples_needed = sum(samples_needed))
```

#### Merge sequence data with lookup

``` r
fungi_merged_table <- merge_data_with_lookup_table(
  bag_data = fungi_full_table,
  bag_var = organism_name,
  lookup_data = fungi_lookup_table,
  lookup_var = Pathogenspecies
)

# Check if the merged table doesn't have any NAs
if (
  length(
    fungi_merged_table %>%
      filter(samples_needed %>% is.na()) %>%
      pull(organism_name) %>% unique()
  ) == 0
) {
  rm(fungi_full_table, fungi_lookup_table)
}
```

#### Sample sequences

``` r
fungi_sampled_table <- fungi_merged_table %>%
  map_sampleing_function(
    col_organism = organism_name, 
    col_id = ID, 
    col_seq = sequence
  )
```

#### Check quality of samples

We check first if enough samples were obtained compares to the needed
samples from the lookup table.

``` r
fungi_sampled_table <- readRDS("../../../data/secreted_data/ready_to_process/sampled_data_without_identical/fungi_sampled_table_good.RDS")
```

``` r
fungi_sampled_table %>% 
  dplyr::group_by(organism_name) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(
    needed_samples = sum(samples_needed),
    actual_samples = sum(samples_obtained)
  )
```

    ## # A tibble: 1 x 2
    ##   needed_samples actual_samples
    ##            <int>          <int>
    ## 1             97             97

Then we check if there are any duplicate sequences within the samples.

``` r
fungi_sampled_table %>% blast_with_ifself(
  col_id = ID,
  col_seq = sequence,
  percent_threshold = 90
) %>% 
  knitr::kable()
```

    ## The data frame has been saved in /var/folders/c2/s7xgl0x53ds16553lpzktdp5hb4zny/T//RtmpM9ugOc/fasta_self.fasta

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 19%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">qseqid</th>
<th style="text-align: right;">qlen</th>
<th style="text-align: left;">sseqid</th>
<th style="text-align: right;">slen</th>
<th style="text-align: right;">length</th>
<th style="text-align: right;">nident</th>
<th style="text-align: right;">mismatch</th>
<th style="text-align: right;">positive</th>
<th style="text-align: right;">percent_identical</th>
<th style="text-align: right;">percent_positive</th>
</tr>
</thead>
<tbody>
</tbody>
</table>

``` r
fungi_sampled_table %>% blast_with_ifself(
  col_id = ID,
  col_seq = sequence,
  percent_threshold = 10
) %>% 
  knitr::kable()
```

    ## The data frame has been saved in /var/folders/c2/s7xgl0x53ds16553lpzktdp5hb4zny/T//RtmpM9ugOc/fasta_self.fasta

| qseqid             |  qlen| sseqid             |  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|:-------------------|-----:|:-------------------|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
| BLGH\_00605-mRNA-1 |   443| KGO36348           |   265|     247|     115|       118|       159|            25.95937|           35.89165|
| KGO36348           |   265| BLGH\_00605-mRNA-1 |   443|     247|     115|       118|       159|            25.95937|           35.89165|
| BLGH\_06958-mRNA-1 |   118| BLGH\_06852-mRNA-1 |   116|     113|      30|        83|        48|            25.42373|           40.67797|
| BLGH\_06852-mRNA-1 |   116| BLGH\_06958-mRNA-1 |   118|     113|      30|        83|        48|            25.42373|           40.67797|
| MGG\_10699T0       |   326| KGO36348           |   265|     289|      79|       164|       121|            24.23313|           37.11656|
| KGO36348           |   265| MGG\_10699T0       |   326|     289|      79|       164|       121|            24.23313|           37.11656|
| BLGH\_06852-mRNA-1 |   116| KGO42382           |    95|      77|      22|        36|        30|            18.96552|           25.86207|
| KGO42382           |    95| BLGH\_06852-mRNA-1 |   116|      77|      22|        36|        30|            18.96552|           25.86207|
| MGG\_10699T0       |   326| BLGH\_00605-mRNA-1 |   443|     265|      77|       138|       119|            17.38149|           26.86230|
| BLGH\_00605-mRNA-1 |   443| MGG\_10699T0       |   326|     242|      71|       126|       110|            16.02709|           24.83070|
| BLGH\_06958-mRNA-1 |   118| BLGH\_05020-mRNA-1 |   106|      42|      15|        26|        21|            12.71186|           17.79661|
| BLGH\_05020-mRNA-1 |   106| BLGH\_06958-mRNA-1 |   118|      42|      15|        26|        21|            12.71186|           17.79661|
| BLGH\_05020-mRNA-1 |   106| CBX91396           |   131|      71|      16|        52|        25|            12.21374|           19.08397|
| EME46666           |    92| MGG\_10280T0       |   110|      39|      13|        25|        19|            11.81818|           17.27273|
| BLGH\_05020-mRNA-1 |   106| BLGH\_06852-mRNA-1 |   116|      41|      13|        25|        18|            11.20690|           15.51724|
| BLGH\_06852-mRNA-1 |   116| BLGH\_05020-mRNA-1 |   106|      41|      13|        25|        18|            11.20690|           15.51724|
| MGG\_05344T0       |   137| CBX92994           |    98|      42|      15|        24|        17|            10.94891|           12.40876|
| BLGH\_06875-mRNA-1 |   100| MGG\_08401T0       |   229|      92|      23|        35|        34|            10.04367|           14.84716|
| MGG\_08401T0       |   229| BLGH\_06875-mRNA-1 |   100|      92|      23|        35|        34|            10.04367|           14.84716|

#### Save results

``` r
# fungi_sampled_table %>%
  # saveRDS(here::here("scripts/r-scripts/getting-secreted-data", "fungi_sampled_table_good01.RDS"))
```

### Run with oomycete

#### Load sequence and lookup data

``` r
# Read the full table of the prediction result of the SignalP
oomycete_full_table <- data.table::fread(here::here("data/secreted_data/updated_signalp_results", "protists_full_table.csv")) %>%
  # Make the name of organisms consistent with the lookup table
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("^_", "") %>%
      stringr::str_to_lower()
  ) %>%
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
oomycete_lookup_table <- data.table::fread(here::here("data/secreted_data/signalp-pipeline", "oomycete_lookup_table.csv")) %>%
  # `colnames<-`(c("Pathogenspecies", "Count")) %>%
  group_by(Pathogenspecies) %>%
  summarise(samples_needed = sum(Count)) %>%
  ungroup()

oomycete_lookup_table <- oomycete_lookup_table %>%
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_to_lower() %>%
      stringr::str_remove_all("\\.") %>%
      stringr::str_replace_all(" ", "_")
  )

oomycete_lookup_table <- oomycete_lookup_table %>%
  # Manual name fixes
  mutate(
    Pathogenspecies = Pathogenspecies %>%
      stringr::str_replace_all("\\(", "") %>%
      stringr::str_replace_all("\\)", "") %>%
      stringr::str_replace_all("-", "_")
  ) %>%
  # Fix repeated Pathogenspecies
  group_by(Pathogenspecies) %>%
  summarise(samples_needed = sum(samples_needed))
```

#### Merge sequence data with lookup

``` r
oomycete_merged_table <- merge_data_with_lookup_table(
  bag_data = oomycete_full_table,
  bag_var = organism_name,
  lookup_data = oomycete_lookup_table,
  lookup_var = Pathogenspecies
)

# Check if the merged table doesn't have any NAs
if (
  length(
    oomycete_merged_table %>%
      filter(samples_needed %>% is.na()) %>%
      pull(organism_name) %>% unique()
  ) == 0
) {
  rm(oomycete_full_table, oomycete_lookup_table)
}
```

#### Sample sequences

``` r
oomycete_sampled_table <- oomycete_merged_table %>%
  map_sampleing_function(
    col_organism = organism_name, 
    col_id = ID, 
    col_seq = sequence
  )
```

#### Check quality of samples

We check first if enough samples were obtained compares to the needed
samples from the lookup table.

``` r
oomycete_sampled_table <- readRDS("../../../scripts/r-scripts/getting-secreted-data/oomycete_sampled_table_good.RDS")
```

``` r
oomycete_sampled_table %>% 
  dplyr::group_by(organism_name) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(
    needed_samples = sum(samples_needed),
    actual_samples = sum(samples_obtained)
  )
```

    ## # A tibble: 1 x 2
    ##   needed_samples actual_samples
    ##            <int>          <int>
    ## 1             85             85

Then we check if there are any duplicate sequences within the samples.

``` r
oomycete_sampled_table %>% blast_with_ifself(
  col_id = ID,
  col_seq = sequence,
  percent_threshold = 90
) %>% 
  knitr::kable()
```

    ## The data frame has been saved in /var/folders/c2/s7xgl0x53ds16553lpzktdp5hb4zny/T//RtmpM9ugOc/fasta_self.fasta

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 9%" />
<col style="width: 19%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">qseqid</th>
<th style="text-align: right;">qlen</th>
<th style="text-align: left;">sseqid</th>
<th style="text-align: right;">slen</th>
<th style="text-align: right;">length</th>
<th style="text-align: right;">nident</th>
<th style="text-align: right;">mismatch</th>
<th style="text-align: right;">positive</th>
<th style="text-align: right;">percent_identical</th>
<th style="text-align: right;">percent_positive</th>
</tr>
</thead>
<tbody>
</tbody>
</table>

``` r
oomycete_sampled_table %>% blast_with_ifself(
  col_id = ID,
  col_seq = sequence,
  percent_threshold = 10
) %>% 
  knitr::kable()
```

    ## The data frame has been saved in /var/folders/c2/s7xgl0x53ds16553lpzktdp5hb4zny/T//RtmpM9ugOc/fasta_self.fasta

| qseqid        |  qlen| sseqid        |  slen|  length|  nident|  mismatch|  positive|  percent\_identical|  percent\_positive|
|:--------------|-----:|:--------------|-----:|-------:|-------:|---------:|---------:|-------------------:|------------------:|
| HpaP805901    |   429| PITG\_12145T0 |   455|     438|     177|       192|       241|            38.90110|           52.96703|
| PITG\_12145T0 |   455| HpaP805901    |   429|     435|     176|       190|       240|            38.68132|           52.74725|
| HpaP800407    |   389| PITG\_14711T0 |   296|     275|     138|       131|       189|            35.47558|           48.58612|
| PITG\_14711T0 |   296| HpaP800407    |   389|     275|     138|       131|       189|            35.47558|           48.58612|
| HpaP811009    |   344| HpaP807735    |   329|     341|     117|       193|       172|            34.01163|           50.00000|
| HpaP807735    |   329| HpaP811009    |   344|     341|     117|       193|       172|            34.01163|           50.00000|
| HpaP807013    |   467| PITG\_01985T0 |   513|     507|     140|       278|       218|            27.29045|           42.49513|
| HpaP807906    |   147| HpaP808319    |   124|      85|      39|        46|        48|            26.53061|           32.65306|
| HpaP808319    |   124| HpaP807906    |   147|      85|      39|        46|        48|            26.53061|           32.65306|
| PITG\_01985T0 |   513| HpaP807013    |   467|     461|     129|       235|       198|            25.14620|           38.59649|
| RAW23090      |   157| PITG\_09213T0 |   151|     156|      39|        78|        61|            24.84076|           38.85350|
| EGZ29639      |   126| HpaP803993    |   129|     133|      31|        84|        52|            24.03101|           40.31008|
| EGZ29639      |   126| EGZ18917      |   127|      99|      29|        51|        43|            22.83465|           33.85827|
| HpaP807906    |   147| EGZ29639      |   126|     122|      31|        80|        52|            21.08844|           35.37415|
| EGZ29639      |   126| HpaP807906    |   147|     122|      31|        80|        52|            21.08844|           35.37415|
| PITG\_09213T0 |   151| RAW23090      |   157|     126|      32|        67|        51|            20.38217|           32.48408|
| PITG\_12779T0 |   302| PITG\_12751T0 |   180|     121|      61|        48|        68|            20.19868|           22.51656|
| EGZ29639      |   126| PITG\_11350T0 |   121|     115|      25|        75|        50|            19.84127|           39.68254|
| HpaP808850    |   190| EGZ16303      |   236|      73|      32|        41|        41|            13.55932|           17.37288|
| EGZ16303      |   236| HpaP808850    |   190|      73|      32|        41|        41|            13.55932|           17.37288|
| PITG\_07488T0 |   104| HpaP807413    |    83|      49|      13|        36|        24|            12.50000|           23.07692|
| RAW23090      |   157| HpaP810226    |   171|      69|      21|        29|        29|            12.28070|           16.95906|
| HpaP807413    |    83| PITG\_07488T0 |   104|      42|      12|        30|        21|            11.53846|           20.19231|
| HpaP808850    |   190| HpaP800963    |   200|      90|      23|        44|        35|            11.50000|           17.50000|
| HpaP800963    |   200| HpaP808850    |   190|      90|      23|        44|        35|            11.50000|           17.50000|
| RAW23090      |   157| HpaP807435    |    98|      49|      16|        21|        26|            10.19108|           16.56051|

#### Save results

``` r
# oomycete_sampled_table %>%
#   saveRDS(here::here("scripts/r-scripts/getting-secreted-data", "oomycete_sampled_table_good01.RDS"))
```
