Getting data for Control sets from NCBI
=======================================

In order to understand the data sets that we will train, we can do
several investigations including getting the control set data and do the
bootstrap resample for the data using the model we developed.

Getting the effector data information
-------------------------------------

``` r
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
phi_base <- data.table::fread("../../../data/getting-data-old/phi-base-main.csv", header = TRUE)


# filter all of the data with 'plant avirulence determinant' information
phi_plant_effector <- phi_base %>%
  dplyr::filter_all(any_vars(str_detect(., 'plant avirulence determinant'))) %>% 
  dplyr::select(`Protein ID`,
                `Pathogen species`, 
                `Pathogen strain`)

phi_plant_effector_to_print <- phi_base %>%
  dplyr::filter_all(any_vars(str_detect(., 'plant avirulence determinant')))

phi_plant_effector_to_print <- phi_plant_effector_to_print %>% 
  head(30)

write.csv(phi_plant_effector_to_print, "phi_plant_effector_to_print.csv")
```

In this `phi_plant_effector` data, we have several cases that we need to
consider:

### Many duplicate protein IDS

``` r
phi_protein_with_duplicated_pathogen_species <- phi_plant_effector %>% 
  group_by(`Protein ID`) %>% 
  arrange(`Protein ID`) %>% 
  dplyr::filter(n() > 1) 

phi_protein_with_duplicated_pathogen_species %>% 
  head(50) %>% 
  knitr::kable()
```

| Protein ID | Pathogen species       | Pathogen strain                                   |
|:-----------|:-----------------------|:--------------------------------------------------|
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A023UJQ9 | Passalora fulva        | no data found                                     |
| A0A0M5K865 | Phytophthora sojae     | P6497                                             |
| A0A0M5K865 | Phytophthora sojae     | P6497                                             |
| A0A0M5K865 | Phytophthora sojae     | P6497                                             |
| A0A0Q0BGR4 | Pseudomonas syringae   | pv. Tomato str DC3000                             |
| A0A0Q0BGR4 | Pseudomonas syringae   | pv. Tomato str DC3000                             |
| A0A0U1JRP7 | Salmonella enterica    | subsp. enterica serovar Typhimurium strain SL1344 |
| A0A0U1JRP7 | Salmonella enterica    | subsp. enterica serovar Typhimurium strain SL1344 |
| A0A2R2Z552 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z552 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z560 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z560 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z561 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z561 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z564 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z564 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z570 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z570 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z572 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z572 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z574 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z574 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z575 | Phytophthora capsici   | Pc537                                             |
| A0A2R2Z575 | Phytophthora capsici   | Pc537                                             |
| A4F4L2     | Leptosphaeria maculans | v23.1.3                                           |
| A4F4L2     | Leptosphaeria maculans | 00-100                                            |
| A4L9T6     | Phytophthora sojae     | P7076                                             |
| A4L9T6     | Phytophthora sojae     | P7076                                             |
| A4LAA7     | Phytophthora sojae     | P7076                                             |
| A4LAA7     | Phytophthora sojae     | P7076                                             |
| A5YTY8     | Phytophthora sojae     | P6497                                             |
| A5YTY8     | Phytophthora sojae     | P6497                                             |
| A6V3V7     | Pseudomonas aeruginosa | PAO1                                              |
| A6V3V7     | Pseudomonas aeruginosa | PAO1                                              |
| A7L812     | Phytophthora sojae     | P7076                                             |
| A7L812     | Phytophthora sojae     | P7076                                             |
| B2C6F3     | Phytophthora sojae     | P6497                                             |
| B2C6F3     | Phytophthora sojae     | P6497                                             |
| B2C6F5     | Phytophthora sojae     | P6954                                             |
| B2C6F5     | Phytophthora sojae     | P6954                                             |
| B2C6F5     | Phytophthora sojae     | P6954                                             |
| B2C6F5     | Phytophthora sojae     | P6954                                             |
| B2C6F5     | Phytophthora sojae     | P6954                                             |
| B2C6F5     | Phytophthora sojae     | P6954                                             |

### Pathogen species that have same protein IDs

``` r
phi_effector_same_ID_diff_pathogen <- phi_plant_effector %>%
  dplyr::filter_all(all_vars(!str_detect(., 'no data found'))) %>% 
  group_by(`Protein ID`) %>%
  dplyr::filter(n_distinct(`Pathogen species`) > 1) %>% 
  dplyr::group_by(`Protein ID`, `Pathogen species`, `Pathogen strain`)  %>%
  slice(1)

phi_effector_same_ID_diff_pathogen %>% 
  knitr::kable()
```

| Protein ID | Pathogen species     | Pathogen strain                   |
|:-----------|:---------------------|:----------------------------------|
| B3VBK9     | Fusarium oxysporum   | f. sp. Lycopersici                |
| B3VBK9     | Passalora fulva      | race 5                            |
| Q2A0P1     | Botrytis cinerea     | B05.10                            |
| Q2A0P1     | Fusarium oxysporum   | f. sp. Lycopersici, strain Fol007 |
| Q2A0P1     | Fusarium oxysporum   | f. sp. Lycopersici 4287           |
| Q2A0P1     | Fusarium oxysporum   | f. sp. lycopersici Fol007         |
| Q2A0P1     | Fusarium oxysporum   | Fo5176                            |
| Q2A0P1     | Pseudomonas syringae | pv. Tomato str DC3000             |
| Q2A0P1     | Verticillium dahliae | JR2                               |

Previous sloopy step
====================

![Flow Chart of Getting the
Data](../data/images/getting-data-flowchart.png)

This is how I get the unique protein ID (to retrieve sequence from and
take the first pathogen species and neglect others.

``` r
# since there are many duplicate pathogen species on 
phi_plant_effector_proteinID_unique <- phi_plant_effector %>%
  group_by(`Protein ID`) %>% 
  summarise(`Pathogen species` = first(`Pathogen species`)) %>% 
  dplyr::filter_all(all_vars(!str_detect(., 'no data found')))
```

However, for the control set, I tried to get all of the list of pathogen
name and strain to get the whole genome sequence IDs. Now, we need to
take the pathogen name and strain for each pathogen species and strain,
so that we can get all of the samples of the protein from NCBI.

``` r
# select only the protein ID data
phi_proteinID <- phi_plant_effector %>%
  dplyr::select(`Protein ID`) %>% 
  unique() %>% 
  dplyr::filter_all(any_vars(!str_detect(., 'no data found')))
```

``` r
# Take any pathogen species without strain
phi_pathogen_strain_no_data_found <- phi_plant_effector %>% 
  dplyr::filter(`Pathogen strain` == 'no data found') %>% 
  group_by(`Pathogen species`) %>% 
  slice(1)

# Take only unique strains
phi_pathogen_strain_unique <- phi_plant_effector%>% 
  group_by(`Pathogen strain`) %>% 
  slice(1) %>% 
  bind_rows(phi_pathogen_strain_no_data_found)
```

``` r
effector_data <- data.table::fread("../../../data/getting-data-old/effector_with_IDs_organism.csv") %>% 
  rename(`Protein ID` = protein_id, `Pathogen species` = pathogen_short)

effector_data %>% 
  mutate(sequence = substr(sequence, 1, 30)) %>% 
  head(20) %>% 
  knitr::kable()
```

| Protein ID | sequence                       | Pathogen species               | name\_src |
|:-----------|:-------------------------------|:-------------------------------|:----------|
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva                | uniprot   |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva                | uniprot   |
| Q258K5     | MVQFKTIFLSTALAALFSTGSSSPATKNNV | Leptosphaeria maculans         | uniprot   |
| Q02039     | MKFLVLPLSLAFLQIGLVFSTPDRCRYTLC | Rhynchosporium commune         | phi\_base |
| Q01145     | MKFNKTIPLYILAFFSTAVIAGGRKWTNKV | Magnaporthe oryzae             | phi\_base |
| Q01144     | MKCNNIILPFALVFFSTTVTAGGGWTNKQF | Magnaporthe oryzae             | phi\_base |
| Q00365     | MLFNAAAAAVFAPLLVMGNVLPRNAGNSPG | Passalora fulva                | uniprot   |
| Q01905     | MNFRALFAATVAALVGSTSATTCTTSQQTV | Phytophthora infestans         | uniprot   |
| Q9C478     | MLFYSLFFFHTVAISAFTNIGTFSHPVYDY | Magnaporthe oryzae             | phi\_base |
| Q6ZX14     | MRDEMWNTATEPIAIIGSGCKFPGGSTTPS | Magnaporthe oryzae             | phi\_base |
| Q5G7L3     | MRLVHAVLLPGIIVFVSNGNLLHAHALHED | Hyaloperonospora arabidopsidis | phi\_base |
| Q709E2     | MAPYSMVLLGALSILGFGAYAQEAAVREPQ | Fusarium oxysporum             | uniprot   |
| Q8NID8     | MKLFILTFIWLLTASEVIAAAKKLPGCDKD | Passalora fulva                | uniprot   |
| Q572D3     | MRLAIMLSATAVAINFATSSAIDQTKVLVY | Phytophthora infestans         | uniprot   |
| Q5VB04     | MQFPTPHLLLTTLLATITTADFSRDCPPGS | Passalora fulva                | uniprot   |
| Q8L552     | MRLSFVLSLVVAIGYVVTCNATEYSDETNI | Phytophthora sojae             | uniprot   |
| Q4VKJ6     | MRVCYFVLVPSVALAVIATESSETSGTIVH | Hyaloperonospora arabidopsidis | phi\_base |
| Q6R661     | MKINLSWKTYAVLSIVSLGGIHAMEHVPAE | Melampsora lini                | uniprot   |
| Q6R659     | MKINLSWKTYAVLSIVSLGGIHAMEHVPAE | Melampsora lini                | uniprot   |
| Q2MV46     | MMMISTKSKLRNIFLLLAIANVNFCVKGHP | Melampsora lini                | uniprot   |

Now, we can map the effector data with all of the strains via protein
IDs

``` r
effector_sequence_with_metadata <- effector_data %>% 
  left_join(., phi_plant_effector, by = 'Protein ID') %>% 
  select(-c('name_src', 'Pathogen species.x')) %>% 
  rename(`Pathogen species` = `Pathogen species.y`)

effector_sequence_with_metadata <- effector_sequence_with_metadata %>% 
  group_by(`Pathogen strain`, `Pathogen species`)

# write.csv(effector_sequence_with_metadata , "../../../data/getting-data-old/effector_sequence_with_metadata.csv", row.names = FALSE)
effector_sequence_with_metadata %>% 
  mutate(sequence = substr(sequence, 1, 30)) %>% 
  head(30) %>% 
  knitr::kable()
```

| Protein ID | sequence                       | Pathogen species | Pathogen strain                        |
|:-----------|:-------------------------------|:-----------------|:---------------------------------------|
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | Race 5                                 |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | Race 0                                 |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | Race 2                                 |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | Race 4                                 |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | Race 5                                 |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | Race 9                                 |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| P22287     | MKLSLLSVELALLIATTLPLCWAAALPVGL | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | Multiple alleles from numerous strains |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | race 4                                 |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |
| Q00363     | MHYTTLLLSTLLVGTALAQPTNPPAKTPKK | Passalora fulva  | no data found                          |

Get the unique pathogen species and pathogen strains:

``` r
effector_sequence_with_metadata_uniq <- effector_sequence_with_metadata %>% 
  dplyr::group_by(`Protein ID`, sequence, `Pathogen species`, `Pathogen strain`) %>% 
  slice(1) %>% 
  ungroup()


effector_sequence_with_metadata_uniq <- effector_sequence_with_metadata_uniq %>% 
  dplyr::select(`Pathogen species`, `Pathogen strain`) %>%
  dplyr::group_by(`Pathogen species`, `Pathogen strain`) %>% 
  slice(1)

effector_sequence_with_metadata_uniq %>% 
  knitr::kable()
```

| Pathogen species               | Pathogen strain                                     |
|:-------------------------------|:----------------------------------------------------|
| Acinetobacter baumannii        | ATCC 17978                                          |
| Aeromonas hydrophila           | SSU                                                 |
| Aeromonas salmonicida          | A449                                                |
| Beauveria bassiana             | ARSEF2860                                           |
| Blumeria graminis              | f. sp. Hordei                                       |
| Blumeria graminis              | f. sp. hordei DH14                                  |
| Botrytis cinerea               | B05.10                                              |
| Botrytis cinerea               | no data found                                       |
| Brucella abortus               | 2308                                                |
| Burkholderia glumae            | 106619                                              |
| Burkholderia pseudomallei      | 10276                                               |
| Burkholderia pseudomallei      | DD503                                               |
| Burkholderia pseudomallei      | K96243                                              |
| Campylobacter jejuni           | 11168                                               |
| Cercospora apii                | no data found                                       |
| Cercospora beticola            | no data found                                       |
| Citrobacter rodentium          | ICC169                                              |
| Clavibacter michiganensis      | Cm15-2.0 sm                                         |
| Colletotrichum orbiculare      | 104-T                                               |
| Coxiella burnetii              | Nine Mile (RSA493)                                  |
| Coxiella burnetii              | Nine Mile Phase II Clone 4 (RSA439)                 |
| Coxiella burnetii              | NMII                                                |
| Coxiella burnetii              | Tn1832                                              |
| Cystobacter fuscus             | DSMZ2262                                            |
| Dothistroma septosporum        | no data found                                       |
| Edwardsiella ictaluri          | no data found                                       |
| Erwinia amylovora              | CFBP 1430                                           |
| Erwinia amylovora              | NCPPB1665                                           |
| Erwinia amylovora              | no data found                                       |
| Escherichia coli               | EDL933                                              |
| Escherichia coli               | TW-XM                                               |
| Francisella tularensis         | SCHU S4                                             |
| Fusarium oxysporum             | f. sp. Lycopersici, strain Fol007                   |
| Fusarium oxysporum             | f. sp. Lycopersici                                  |
| Fusarium oxysporum             | f. sp. Lycopersici 007                              |
| Fusarium oxysporum             | f. sp. Lycopersici 4287                             |
| Fusarium oxysporum             | f. sp. Lycopersici Chiba-5                          |
| Fusarium oxysporum             | f. sp. lycopersici Fol007                           |
| Fusarium oxysporum             | f. sp. Lycopersici Fol007                           |
| Fusarium oxysporum             | f. sp.melonis 26406                                 |
| Fusarium oxysporum             | Fo5176                                              |
| Fusarium oxysporum             | Fol                                                 |
| Globodera rostochiensis        | Ro1-Mierenbos                                       |
| Helicobacter pylori            | 51B                                                 |
| Heterodera glycines            | no data found                                       |
| Hyaloperonospora arabidopsidis | Cala2                                               |
| Hyaloperonospora arabidopsidis | Emco5                                               |
| Hyaloperonospora arabidopsidis | Emoy2                                               |
| Hyaloperonospora arabidopsidis | Emoy2?                                              |
| Hyaloperonospora arabidopsidis | Emwa1                                               |
| Hyaloperonospora arabidopsidis | Maks9                                               |
| Hyaloperonospora arabidopsidis | no data found                                       |
| Hyaloperonospora arabidopsidis | Noks1                                               |
| Hyaloperonospora arabidopsidis | Waco9                                               |
| Legionella pneumophila         | 130b                                                |
| Legionella pneumophila         | 130b serogroup 1 (ATCC BAA-74)                      |
| Legionella pneumophila         | no data found                                       |
| Leptosphaeria maculans         | 00-100                                              |
| Leptosphaeria maculans         | IBCN 14                                             |
| Leptosphaeria maculans         | no data found                                       |
| Leptosphaeria maculans         | Nzt-4                                               |
| Leptosphaeria maculans         | v23.1.3                                             |
| Leptosphaeria maculans         | v23.1.3Xv37.1. 4                                    |
| Listeria monocytogenes         | EGD-e                                               |
| Macrosiphum euphorbiae         | no data found                                       |
| Magnaporthe oryzae             | 4091-5-8                                            |
| Magnaporthe oryzae             | 4360-R-62                                           |
| Magnaporthe oryzae             | 4392-1-6                                            |
| Magnaporthe oryzae             | 70-15                                               |
| Magnaporthe oryzae             | 84R-62B                                             |
| Magnaporthe oryzae             | 85-14B1                                             |
| Magnaporthe oryzae             | 98-06                                               |
| Magnaporthe oryzae             | CP987                                               |
| Magnaporthe oryzae             | Guy11                                               |
| Magnaporthe oryzae             | Ina168                                              |
| Magnaporthe oryzae             | INA168                                              |
| Magnaporthe oryzae             | Ina72                                               |
| Magnaporthe oryzae             | Ina86-137                                           |
| Magnaporthe oryzae             | JS153                                               |
| Magnaporthe oryzae             | Ku80                                                |
| Magnaporthe oryzae             | KV103                                               |
| Magnaporthe oryzae             | no data found                                       |
| Magnaporthe oryzae             | O-137                                               |
| Magnaporthe oryzae             | OS99-G-7a                                           |
| Magnaporthe oryzae             | PWL2                                                |
| Magnaporthe oryzae             | R88-002                                             |
| Magnaporthe oryzae             | RB22                                                |
| Magnaporthe oryzae             | RO1-1                                               |
| Magnaporthe oryzae             | sasa2                                               |
| Magnaporthe oryzae             | Sasa2                                               |
| Magnaporthe oryzae             | TH68-141                                            |
| Melampsora lini                | no data found                                       |
| Mycobacterium tuberculosis     | H37Rv                                               |
| Pantoea stewartii              | subsp. Stewartii                                    |
| Pantoea stewartii              | subsp. Stewartii DC283                              |
| Parastagonospora nodorum       | no data found                                       |
| Parastagonospora nodorum       | Sn2000                                              |
| Passalora fulva                | Multiple alleles from numerous strains              |
| Passalora fulva                | no data found                                       |
| Passalora fulva                | Race 0                                              |
| Passalora fulva                | race 2                                              |
| Passalora fulva                | Race 2                                              |
| Passalora fulva                | race 4                                              |
| Passalora fulva                | Race 4                                              |
| Passalora fulva                | race 5                                              |
| Passalora fulva                | Race 5                                              |
| Passalora fulva                | Race 9                                              |
| Penicillium expansum           | PE 100/PEX2                                         |
| Phytophthora cactorum          | P381                                                |
| Phytophthora capsici           | Pc537                                               |
| Phytophthora infestans         | 88069                                               |
| Phytophthora infestans         | 90128                                               |
| Phytophthora infestans         | HP10-45                                             |
| Phytophthora infestans         | no data found                                       |
| Phytophthora infestans         | T30-4                                               |
| Phytophthora parasitica        | 1828                                                |
| Phytophthora parasitica        | var. nicotianae race 0                              |
| Phytophthora sojae             | no data found                                       |
| Phytophthora sojae             | P6497                                               |
| Phytophthora sojae             | P6497 (race 2)                                      |
| Phytophthora sojae             | P6954                                               |
| Phytophthora sojae             | P7076                                               |
| Phytophthora sojae             | race 2 strain P6497                                 |
| Phytophthora sojae             | race 7 strain P7064                                 |
| Pseudocercospora fuligena      | no data found                                       |
| Pseudomonas aeruginosa         | PAO1                                                |
| Pseudomonas cichorii           | JBC1                                                |
| Pseudomonas savastanoi         | pv. savastanoi NCPPB 3335                           |
| Pseudomonas syringae           | no data found                                       |
| Pseudomonas syringae           | pv. glycinea race 1                                 |
| Pseudomonas syringae           | pv. glycinea race 3                                 |
| Pseudomonas syringae           | pv. glycinea race 4                                 |
| Pseudomonas syringae           | pv. glycinea race 6                                 |
| Pseudomonas syringae           | pv. maculicola                                      |
| Pseudomonas syringae           | pv. maculicola M2                                   |
| Pseudomonas syringae           | pv. maculicola str. ES4326                          |
| Pseudomonas syringae           | pv. maculicola str. M6                              |
| Pseudomonas syringae           | pv. phaseolicola NPS3121                            |
| Pseudomonas syringae           | pv. phaseolicola 1302A (race 4)                     |
| Pseudomonas syringae           | pv. phaseolicola 1448A                              |
| Pseudomonas syringae           | pv. phaseolicola 1449A race5 & 7                    |
| Pseudomonas syringae           | pv. phaseolicola 1449B (race7)                      |
| Pseudomonas syringae           | pv. phaseolicola race 3                             |
| Pseudomonas syringae           | pv. pisi 151                                        |
| Pseudomonas syringae           | pv. pisi 870A                                       |
| Pseudomonas syringae           | pv. syringae                                        |
| Pseudomonas syringae           | pv. syringae 61                                     |
| Pseudomonas syringae           | pv. syringae B728a                                  |
| Pseudomonas syringae           | pv. syringae B728A                                  |
| Pseudomonas syringae           | pv. tomato JL1065                                   |
| Pseudomonas syringae           | pv. tomato PT23                                     |
| Pseudomonas syringae           | pv. Tomato str DC3000                               |
| Puccinia striiformis           | f. sp. tritici                                      |
| Pythium aphanidermatum         | no data found                                       |
| Ralstonia solanacearum         | GMI1000                                             |
| Ralstonia solanacearum         | no data found                                       |
| Ralstonia solanacearum         | RS1000                                              |
| Ralstonia solanacearum         | SL2029                                              |
| Ralstonia solanacearum         | SL341                                               |
| Ralstonia solanacearum         | UW551                                               |
| Rhynchosporium commune         | no data found                                       |
| Rhynchosporium commune         | UK7                                                 |
| Salmonella enterica            | subsp. enterica serovar Typhi str. Ty2              |
| Salmonella enterica            | subsp. enterica serovar Typhimurium                 |
| Salmonella enterica            | subsp. enterica serovar Typhimurium str. 14028s     |
| Salmonella enterica            | subsp. enterica serovar Typhimurium str. ATCC 14028 |
| Salmonella enterica            | subsp. enterica serovar Typhimurium str. LT2        |
| Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028    |
| Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344   |
| Shigella flexneri              | M90T-Sm                                             |
| Staphylococcus aureus          | SA564                                               |
| Toxoplasma gondii              | RH                                                  |
| Ustilago maydis                | SG200                                               |
| Verticillium dahliae           | 12008                                               |
| Verticillium dahliae           | JR2                                                 |
| Verticillium dahliae           | V592                                                |
| Verticillium dahliae           | VdLs17                                              |
| Vibrio parahaemolyticus        | RIMD2210633                                         |
| Xanthomonas axonopodis         | pv. citri                                           |
| Xanthomonas axonopodis         | pv. manihotis CFBP1851                              |
| Xanthomonas campestris         | pv. Campestris B186                                 |
| Xanthomonas campestris         | pv. Campestris str. 8004                            |
| Xanthomonas campestris         | pv. vesicatoria str. 85-10                          |
| Xanthomonas citri              | subsp. citri 306                                    |
| Xanthomonas citri              | subsp. Citri 306                                    |
| Xanthomonas citri              | subsp. Malvacearum, strain H1005                    |
| Xanthomonas oryzae             | pv. oryzae BXO1                                     |
| Xanthomonas oryzae             | pv. oryzae PXO99                                    |
| Xanthomonas oryzae             | pv. oryzae PXO99A                                   |
| Xanthomonas oryzae             | pv. oryzae X11-5A                                   |
| Xanthomonas oryzae             | pv. oryzicola BLS256                                |
| Xanthomonas oryzae             | pv. Oryzicola BLS256                                |
| Xanthomonas oryzae             | pv. Oryzicola MAI10                                 |
| Xanthomonas oryzae             | pv. oryzicola RS105                                 |
| Xanthomonas oryzae             | X11-5A                                              |
| Xylella fastidiosa             | EB92-1                                              |
| Yersinia enterocolitica        | JB580v                                              |
| Yersinia pestis                | 91001                                               |
| Yersinia pestis                | KIM1001                                             |
| Yersinia pseudotuberculosis    | 32777                                               |
| Yersinia pseudotuberculosis    | YPIII(pIB102)                                       |
| Zymoseptoria tritici           | IPO323 deltaKu80                                    |
