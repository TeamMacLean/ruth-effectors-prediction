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

bacteria_full_table <- data.table::fread("../../../data/secreted_data/signalp-pipeline/bacteria_full_table.csv") %>%
  # Make the name of organisms consistent with the lookup table
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("^_", "")  %>%
      stringr::str_to_lower()
  )


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
  saveRDS("../../../data/secreted_data/data_processed_after_signalp/bacteria_non_eff_secreted_data.RDS")
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
```

### Check the debug on the sampling data

``` r
bacteria_missing_organisms <- left_join(
  bacteria_lookup_table %>%
    dplyr::mutate_at(
      .vars = dplyr::vars(Pathogenspecies),
      .funs = as.character
    ) %>%
    dplyr::select(base_names = Pathogenspecies, total_count_src = total_count),
  bacteria_merged_table,
  by = "base_names"
)

bacteria_missing_organisms %>%
  group_by(base_names) %>%
  slice(1) %>%
  # filter(is.na(organism_name)) %>%
  select(base_names, total_count_src, organism_name)
```

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
      stringr::str_replace_all(" ", "_")
  )

# Filter for the 
not_available_fungi <- c("cercospora_apii", "cercospora_beticola", "pyrenophora_tritici-repentis", "melampsora_lini", "parastagonospora_nodorum", "passalora_fulva", "pseudocercospora_fuligena")

fungi_lookup_table <- fungi_lookup_table %>% 
  dplyr::filter(!(Pathogenspecies %in% not_available_fungi))

# Read the data file results
fungi_full_table <-
  data.table::fread("../../../data/secreted_data/signalp-pipeline/fungi_full_table.csv")  %>%
  # Filter only data with signal peptide and with prediction > 0.9
  dplyr::filter(prediction == "Signal peptide", signalp_prob > 0.9)
```

``` r
# Check teh availability of "cercospora_beticola" signalP prediction

fungi_full_table <-
  data.table::fread("../../../data/secreted_data/signalp-pipeline/fungi_full_table.csv") 
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
  saveRDS("../../../data/secreted_data/data_processed_after_signalp/fungi_non_eff_secreted_data.RDS")
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

    ## [1] 94

``` r
# Compare the lookup table and the result of sampling
data.frame(table = c("result", "lookup"), 
           count_seq = c(fungi_non_eff_secreted_count %>% select(count) %>% 
                           sum(), 
                         fungi_lookup_table %>% select(total_count) %>% 
                           sum()))
```

    ##    table count_seq
    ## 1 result        94
    ## 2 lookup        94

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

    ## # A tibble: 0 x 3
    ## # Groups:   base_names [0]
    ## # … with 3 variables: base_names <chr>, total_count_src <int>,
    ## #   organism_name <chr>

``` r
fungi_debug %>% 
  group_by(organism_name) %>% 
  summarise(count_organism = n()) %>% 
  knitr::kable()
```

| organism\_name                                      |  count\_organism|
|:----------------------------------------------------|----------------:|
| blumeria\_graminis\_ef2                             |               14|
| botrytis\_cinerea\_asm83294v1                       |                1|
| colletotrichum\_orbiculare\_gca\_000350065\_1       |                1|
| dothistroma\_septosporum\_gca\_000340195\_1         |                1|
| fusarium\_graminearum\_rr1                          |                1|
| fusarium\_oxysporum\_fo2                            |                9|
| leptosphaeria\_maculans\_asm23037v1                 |                6|
| magnaporthe\_oryzae\_mg8                            |               19|
| penicillium\_expansum\_gca\_000769735\_asm76973v1   |                4|
| puccinia\_striiformis\_pst-130\_1\_0                |                2|
| rhynchosporium\_commune\_gca\_900074885\_version\_1 |                3|
| ustilago\_maydis\_umaydis521\_2\_0                  |               25|
| verticillium\_dahliae\_asm15067v2                   |                6|
| zymoseptoria\_tritici\_mg2                          |                2|

Oomycete
--------

``` r
oomycete_lookup_table <- data.table::fread("../../../data/secreted_data/signalp-pipeline/oomycete_lookup_table.csv") %>% `colnames<-`(c("Pathogenspecies", "total_count"))


# Read the data file results
oomycete_full_table <-
  data.table::fread("../../../data/secreted_data/signalp-pipeline/protist_full_table.csv") %>%
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
oomyecte_full_table <- oomycete_full_table %>%
  mutate(
    organism_name = organism_name %>%
      stringr::str_replace_all("\\.", "_") %>%
      stringr::str_to_lower()
  )
```

``` r
oomycete_non_eff_secreted_data <- sample_from_lookup_table(
  data = oomyecte_full_table,
  data_var = organism_name,
  lookup = oomycete_lookup_table,
  lookup_var = Pathogenspecies
)
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

``` r
oomycete_non_eff_secreted_data %>% 
  saveRDS("../../../data/secreted_data/data_processed_after_signalp/oomycete_non_eff_secreted_data.RDS")
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
    ## 2    fungi      94
    ## 3 oomycete      95

``` r
all_non_eff_secreted <- rbind(
  bacteria_non_eff_secreted_data %>% dplyr::select(sequence),
  fungi_non_eff_secreted_data %>% dplyr::select(sequence),
  oomycete_non_eff_secreted_data %>% dplyr::select(sequence)
)
```

``` r
all_non_eff_secreted %>%
  head(10) %>%
  knitr::kable()
```

<table>
<colgroup>
<col style="width: 100%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">sequence</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">MSTLARTAHVTRQTSESTIDLQLDLDGTGASEISTSVPFYDHMLTAFAKHSLTDLRVTATGDTHIDVHHTVEDVGIVLGQAIREALGDKSGIARFGDALVPLDEALVQSVVDISGRPFLVHSGEPAGFEMHLIGGHFTGSMVRHVFEAITFHAGLTVHVTVLGGRDPHHIAEAEFKSFARAFRQAKELDPRVSGIPSTKGAL</td>
</tr>
<tr class="even">
<td style="text-align: left;">MDYDDNDTVGGYAIPVDPMELLQCDSCQ</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MALAKHPPEELQVLLFSHDADMPAVETFLGGPPDPALHLRLDAGKRAAHAFGVDTLPTSILVVDGRLVARFQGPREWDSRAMRRLLEKLTEEHPARDPAPVH</td>
</tr>
<tr class="even">
<td style="text-align: left;">MIPDLSQALLWLEAHPDALKGIGRGIERETLRVSPDGFLATTAHPASLGSALTHKWITTDFAEALLEFITPVDRDIDHLLAFLRDIHRHVARELGEERMWPMSMPCRVDPEQDIELAQYGSSNIGQMKTLYRQGLKNRYGALMQIISGVHYNFSLPLSFWQAWAGVKDAESGKEAISAGYLQLIRNYYRFGWVIPYLFGASPAVCSSFLQGKESTLPFERSDKGLLSLPFATSLRLSDLGYTNKSQSSLEMTFNHLPDYIAGLKAAIKTPSEEFAAMGVKNSAGDWLQLNTNVLQIENELYAPVRPKRVTRSGEAPSDALQRGGIEYIEVRSLDINPFSPTGVDADQVRFLDLFLIWCTLADAPEMNSAELACTRKNWDLVITEGRKPGLTLVMGCSETQYALVDLGKMLFNDLRRIAETLDSQGDNQLYQQVCDRLIASFDDPDLTCSARFLQMLKENGIEATGLALAGQYRAQLCAEPLEVLTGQRFSDEAQRSRFSQQETEESDTLSFEAFLQSKS</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MKTTLSQSFIINKLSIHVKPELNSSGKVVFEPNPQQKPYIVFDDHRDSPVGFGVKVSLTKKTYVIQRRVSSGTRGVKEGKKPSSVLKVKVGNVSDFPSIDQAREAARQLVQTMIVTKRNPNRIKRETEASELTISEVFAQYRQHLLGRSKPAKPNTLSVLDKAENRLKEWEALRVKDLTGNEILRKFDEIASRARTAAEQTFRWINVAVKHAIEIEAGNAQTQQRPPTLSYNPFNILKVQKKFRTRSELEDSYRAKGVRNPLSPKDTLGRFLTALHNKRSFNRLGCDYLLLRWPNSLCHYSDIIHIYVAD</td>
</tr>
<tr class="even">
<td style="text-align: left;">MTSPVRSTPAIRLAFALQMVLNAGLIVLACILIIFLGKETMHLGNVLLNTGEQTSSYLLIDGIVIYFLYFEFIALIIKYFQSGYHFPLRYFVYIGITAIIRLIIVDHKNPFDTLAYAIAILILVITLWLANTNRLKRE</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MTTVIHLLGADIPHHNQTVLGFFDRTLYHEVPTAAARQFWLVSGQPERAEQFPRLHIRVFADKAALASAVIDAGCQQRELRFFCHGQFNPRLWLALLSGKLRRNQVYWHIWGADLYEDAAGLKFRLFYLMRRLAQKRVAHVFATRGDIHRFHQRNPQIPASLLYFPTRLAQHAVTSEQPAGPLTILLGNSGDASNRHIQALQQIQRQFGSEIKVVVPLGYPLNNEAYIQRVRAVADDLFAPGTVQLLTEKLAFSDYLQLLSRCHLGYFLFQRQQGIGTLCLLMQANVPFVVSRKNPFWRDLVEQQVPVLFTSDRLDVATIGEARRQMQLLDKRHIAFFEPGYLAGWRQALRLAEGDNA</td>
</tr>
<tr class="even">
<td style="text-align: left;">MENAKQSFQDVLEFVRLFRRKNKLQREIHDNEKKIRDNQKRVLLLDNLSEYIKPGMSIEAVQEIITSMRSDYEDRVDEYIIKVADLSKERRDLSKKLKTLGEVKG</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MAKGSSDFIPLNVAVMTVSDRHSADSDSSGSYLQEALAAVGHQVLDRAIVPDNRYRIRAVVSSWIASTDVQVVLINGGTGFNDTNATPEAIGPLFDREVAGFGELFRMISWEDIATSTLQSRALAGIANGTLIFAIPGSTAACQLAWERIIQHQLDARTRPCHFVSHLKKT</td>
</tr>
<tr class="even">
<td style="text-align: left;">MLRILCLWMALISPLAAAQPPLTLAFDNAPVAQVLQALADYQQLNLVVAPGVEGHLSLRLKDVPWRQALQLVIKTGKLSMEQQGNVLMVFPASWQQEEQRKAAVQREEQEKMLPLFDRILTLVHADAASVNTSLQSERARLMSARGSITLDSRTNSLLLRDTEKALAQTEHWVRQLDVPLEQIELAAQIVTMSEESLRELGVTWGMSGEAQVADALRTSQLRVDLASGRPAGVAGLTLAKLDGRLLELELSALESEHQADIIASPRLFTSHQQTASIKQGTEIPYEVSTGNSGSTTMEFKEAMLAMEVTPLVQSNGRILLKLHISQNVPGRNMRSGENEVLTIDKQEIDTQVMLKDGQTLALGGIFQQQSATGSTRVPWLGDIPLLGSLFRHDVRQQKRRELVIFITPRLVRDD</td>
</tr>
</tbody>
</table>

List of fungi that does not exist

``` r
not_available_fungi <- c("cercospora_apii", "cercospora_beticola", "melampsora_lini", "parastagonospora_nodorum", "passalora_fulva", "pseudocercospora_fuligena")
```

``` r
fungi_lookup_table %>%
  dplyr::filter(Pathogenspecies %in% not_available_fungi) 
```

    ## [1] Pathogenspecies total_count    
    ## <0 rows> (or 0-length row.names)

``` r
not_available_oomycete <- c("phytophthora_capsici")
```

``` r
oomycete_lookup_table %>%
  dplyr::filter(Pathogenspecies %in% not_available_oomycete)
```

    ##        Pathogenspecies total_count
    ## 1 phytophthora_capsici          10
