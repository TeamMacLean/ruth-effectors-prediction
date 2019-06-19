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

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.1       ✔ purrr   0.3.2  
    ## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
phi_base <- data.table::fread("../../data/phi-base-main.csv", header = TRUE)

# filter all of the data with 'plant avirulence determinant' information
phi_plant_effector <- phi_base %>%
  dplyr::filter_all(any_vars(str_detect(., 'plant avirulence determinant'))) %>% 
  dplyr::select(`Protein ID`,
                `Pathogen species`, 
                `Pathogen strain`)
```

In this `phi_plant_effector` data, we have several cases that we need to
consider:

1.  Many duplicate protein IDS

``` r
phi_protein_with_duplicated_pathogen_species <- phi_plant_effector %>% 
  group_by(`Protein ID`) %>% 
  arrange(`Protein ID`) %>% 
  dplyr::filter(n() > 1) 

phi_protein_with_duplicated_pathogen_species %>% 
  knitr::kable()
```

| Protein ID    | Pathogen species               | Pathogen strain                                   |
|:--------------|:-------------------------------|:--------------------------------------------------|
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A023UJQ9    | Passalora fulva                | no data found                                     |
| A0A0M5K865    | Phytophthora sojae             | P6497                                             |
| A0A0M5K865    | Phytophthora sojae             | P6497                                             |
| A0A0M5K865    | Phytophthora sojae             | P6497                                             |
| A0A0Q0BGR4    | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| A0A0Q0BGR4    | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| A0A0U1JRP7    | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| A0A0U1JRP7    | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| A0A2R2Z552    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z552    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z560    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z560    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z561    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z561    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z564    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z564    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z570    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z570    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z572    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z572    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z574    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z574    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z575    | Phytophthora capsici           | Pc537                                             |
| A0A2R2Z575    | Phytophthora capsici           | Pc537                                             |
| A4F4L2        | Leptosphaeria maculans         | v23.1.3                                           |
| A4F4L2        | Leptosphaeria maculans         | 00-100                                            |
| A4L9T6        | Phytophthora sojae             | P7076                                             |
| A4L9T6        | Phytophthora sojae             | P7076                                             |
| A4LAA7        | Phytophthora sojae             | P7076                                             |
| A4LAA7        | Phytophthora sojae             | P7076                                             |
| A5YTY8        | Phytophthora sojae             | P6497                                             |
| A5YTY8        | Phytophthora sojae             | P6497                                             |
| A6V3V7        | Pseudomonas aeruginosa         | PAO1                                              |
| A6V3V7        | Pseudomonas aeruginosa         | PAO1                                              |
| A7L812        | Phytophthora sojae             | P7076                                             |
| A7L812        | Phytophthora sojae             | P7076                                             |
| B2C6F3        | Phytophthora sojae             | P6497                                             |
| B2C6F3        | Phytophthora sojae             | P6497                                             |
| B2C6F5        | Phytophthora sojae             | P6954                                             |
| B2C6F5        | Phytophthora sojae             | P6954                                             |
| B2C6F5        | Phytophthora sojae             | P6954                                             |
| B2C6F5        | Phytophthora sojae             | P6954                                             |
| B2C6F5        | Phytophthora sojae             | P6954                                             |
| B2C6F5        | Phytophthora sojae             | P6954                                             |
| B2C6F8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| B2C6F8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| B2SJG5        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| B2SJG5        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| B2SJG5        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| B2SU53        | Xanthomonas oryzae             | pv. oryzae PXO99                                  |
| B3VBK9        | Passalora fulva                | race 5                                            |
| B3VBK9        | Fusarium oxysporum             | f. sp. Lycopersici                                |
| B3VBK9        | Passalora fulva                | no data found                                     |
| B3VBK9        | Passalora fulva                | no data found                                     |
| B9A9V1        | Magnaporthe oryzae             | CP987                                             |
| B9A9V1        | Magnaporthe oryzae             | 4091-5-8                                          |
| B9WZW9        | Magnaporthe oryzae             | JS153                                             |
| B9WZW9        | Magnaporthe oryzae             | Ina168                                            |
| B9WZW9        | Magnaporthe oryzae             | Ina86-137                                         |
| B9WZW9        | Magnaporthe oryzae             | TH68-141                                          |
| B9ZUL0        | Leptosphaeria maculans         | v23.1.3                                           |
| B9ZUL0        | Leptosphaeria maculans         | v23.1.3Xv37.1. 4                                  |
| B9ZUL0        | Leptosphaeria maculans         | v23.1.3Xv37.1. 4                                  |
| B9ZUL0        | Leptosphaeria maculans         | Nzt-4                                             |
| C0LT61        | Fusarium oxysporum             | f. sp. Lycopersici 4287                           |
| C0LT61        | Fusarium oxysporum             | f. sp. lycopersici Fol007                         |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C0SPN9        | Ralstonia solanacearum         | no data found                                     |
| C1KJH9        | Magnaporthe oryzae             | O-137                                             |
| C1KJH9        | Magnaporthe oryzae             | O-137                                             |
| C4B8B7        | Magnaporthe oryzae             | Ina168                                            |
| C4B8B7        | Magnaporthe oryzae             | sasa2                                             |
| C4B8B7        | Magnaporthe oryzae             | Ina86-137                                         |
| C4B8B7        | Magnaporthe oryzae             | sasa2                                             |
| C4B8B7        | Magnaporthe oryzae             | Sasa2                                             |
| C4B8B7        | Magnaporthe oryzae             | Sasa2                                             |
| C4B8B7        | Magnaporthe oryzae             | INA168                                            |
| C4B8B7        | Magnaporthe oryzae             | INA168                                            |
| C4B8B7        | Magnaporthe oryzae             | INA168                                            |
| C4B8B7        | Magnaporthe oryzae             | INA168                                            |
| C4B8B7        | Magnaporthe oryzae             | INA168                                            |
| C4B8B7        | Magnaporthe oryzae             | INA168                                            |
| C4B8B8        | Magnaporthe oryzae             | 84R-62B                                           |
| C4B8B8        | Magnaporthe oryzae             | Sasa2                                             |
| C4B8B9        | Magnaporthe oryzae             | Ina168                                            |
| C4B8B9        | Magnaporthe oryzae             | sasa2                                             |
| C4B8B9        | Magnaporthe oryzae             | Ina72                                             |
| C4B8B9        | Magnaporthe oryzae             | sasa2                                             |
| C4B8B9        | Magnaporthe oryzae             | Ina72                                             |
| C4B8B9        | Magnaporthe oryzae             | sasa2                                             |
| C4B8B9        | Magnaporthe oryzae             | Ina72                                             |
| C4B8B9        | Magnaporthe oryzae             | 84R-62B                                           |
| C5A9K2        | Burkholderia glumae            | 106619                                            |
| C5A9K2        | Burkholderia glumae            | 106619                                            |
| C5A9K2        | Burkholderia glumae            | 106619                                            |
| C6ZEZ6        | Magnaporthe oryzae             | Guy11                                             |
| C6ZEZ6        | Magnaporthe oryzae             | RO1-1                                             |
| C6ZEZ6        | Magnaporthe oryzae             | RB22                                              |
| C9WMG8        | Fusarium oxysporum             | f. sp. Lycopersici 007                            |
| C9WMG8        | Fusarium oxysporum             | f. sp. lycopersici Fol007                         |
| D0MZL5        | Phytophthora infestans         | 88069                                             |
| D0MZL5        | Phytophthora infestans         | 88069                                             |
| D0N0I4        | Phytophthora infestans         | no data found                                     |
| D0N0I4        | Phytophthora infestans         | 88069                                             |
| D0U2D1        | Fusarium oxysporum             | f. sp. Lycopersici Fol007                         |
| D0U2D1        | Fusarium oxysporum             | f. sp. Lycopersici Fol007                         |
| D0U2D1        | Fusarium oxysporum             | f. sp. Lycopersici Fol007                         |
| D2TT37        | Citrobacter rodentium          | ICC169                                            |
| D2TT37        | Citrobacter rodentium          | ICC169                                            |
| D2Z000        | Ralstonia solanacearum         | RS1000                                            |
| D2Z000        | Ralstonia solanacearum         | RS1000                                            |
| D2Z001        | Ralstonia solanacearum         | RS1000                                            |
| D2Z001        | Ralstonia solanacearum         | RS1000                                            |
| D4I0P2        | Erwinia amylovora              | NCPPB1665                                         |
| D4I0P2        | Erwinia amylovora              | CFBP 1430                                         |
| D4I0Q8        | Erwinia amylovora              | NCPPB1665                                         |
| D4I0Q8        | Erwinia amylovora              | CFBP 1430                                         |
| D4I218        | Erwinia amylovora              | NCPPB1665                                         |
| D4I218        | Erwinia amylovora              | CFBP 1430                                         |
| D4I227        | Erwinia amylovora              | NCPPB1665                                         |
| D4I227        | Erwinia amylovora              | CFBP 1430                                         |
| D4I230        | Erwinia amylovora              | NCPPB1665                                         |
| D4I230        | Erwinia amylovora              | CFBP 1430                                         |
| D6BK26        | Ralstonia solanacearum         | SL341                                             |
| D6BK26        | Ralstonia solanacearum         | SL2029                                            |
| D6BK26        | Ralstonia solanacearum         | SL341                                             |
| D6BK26        | Ralstonia solanacearum         | SL2029                                            |
| E0W4Y2        | Phytophthora sojae             | P6497 (race 2)                                    |
| E0W4Y2        | Phytophthora sojae             | P6497 (race 2)                                    |
| E0W5A0        | Phytophthora sojae             | P6497                                             |
| E0W5A0        | Phytophthora sojae             | P6497                                             |
| E0W5A0        | Phytophthora sojae             | P6497                                             |
| E2MM96        | Pseudomonas syringae           | no data found                                     |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MM96        | Pseudomonas syringae           | pv. glycinea race 4                               |
| E2MMQ4        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| E2MMQ4        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| F8UNP0        | Phytophthora sojae             | P6497                                             |
| F8UNP0        | Phytophthora sojae             | P6497                                             |
| F8UNP0        | Phytophthora sojae             | P6497                                             |
| F8UNP0        | Phytophthora sojae             | P6497                                             |
| F8UNP0        | Phytophthora sojae             | P6497                                             |
| G0X840        | Ustilago maydis                | SG200                                             |
| G0X840        | Ustilago maydis                | SG200                                             |
| G0X840        | Ustilago maydis                | SG200                                             |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1FM79        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G1JUH0        | Globodera rostochiensis        | Ro1-Mierenbos                                     |
| G1JUH0        | Globodera rostochiensis        | Ro1-Mierenbos                                     |
| G1JUH0        | Globodera rostochiensis        | Ro1-Mierenbos                                     |
| G2WSX9        | Verticillium dahliae           | JR2                                               |
| G2WSX9        | Verticillium dahliae           | JR2                                               |
| G2WSX9        | Verticillium dahliae           | JR2                                               |
| G2WSX9        | Verticillium dahliae           | JR2                                               |
| G2X444        | Verticillium dahliae           | JR2                                               |
| G2X444        | Verticillium dahliae           | JR2                                               |
| G2X444        | Verticillium dahliae           | JR2                                               |
| G2X444        | Verticillium dahliae           | JR2                                               |
| G2X4U8        | Verticillium dahliae           | VdLs17                                            |
| G2X4U8        | Verticillium dahliae           | VdLs17                                            |
| G2X4U8        | Verticillium dahliae           | VdLs17                                            |
| G2X4U8        | Verticillium dahliae           | VdLs17                                            |
| G2XA95        | Verticillium dahliae           | JR2                                               |
| G2XA95        | Verticillium dahliae           | JR2                                               |
| G2XA95        | Verticillium dahliae           | JR2                                               |
| G2XA95        | Verticillium dahliae           | JR2                                               |
| G3C9N5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9N9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P1        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P1        | Hyaloperonospora arabidopsidis | Emoy2?                                            |
| G3C9P1        | Hyaloperonospora arabidopsidis | Emoy2?                                            |
| G3C9P1        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P1        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P1        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P9        | Hyaloperonospora arabidopsidis | Waco9                                             |
| G3C9P9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9P9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9Q5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9Q5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9Q5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9Q5        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9Q9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9Q9        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R4        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R4        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R4        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R4        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9R4        | Hyaloperonospora arabidopsidis | no data found                                     |
| G3C9R4        | Hyaloperonospora arabidopsidis | no data found                                     |
| G3C9R4        | Hyaloperonospora arabidopsidis | no data found                                     |
| G3C9R4        | Hyaloperonospora arabidopsidis | no data found                                     |
| G3C9S3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9S3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9S3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9S3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9S3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T3        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3C9T8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| G3XDC5        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| G3XDC5        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| G4MVX4        | Magnaporthe oryzae             | 70-15                                             |
| G4MVX4        | Magnaporthe oryzae             | 70-15                                             |
| G4MX34        | Magnaporthe oryzae             | 70-15                                             |
| G4MX34        | Magnaporthe oryzae             | 70-15                                             |
| G4N8Y3        | Magnaporthe oryzae             | 70-15                                             |
| G4N8Y3        | Magnaporthe oryzae             | 70-15                                             |
| G4ZRQ8        | Phytophthora sojae             | P6954                                             |
| G4ZRQ8        | Phytophthora sojae             | P6954                                             |
| G4ZRQ8        | Phytophthora sojae             | P6954                                             |
| G4ZRQ8        | Phytophthora sojae             | P6954                                             |
| G4ZRQ8        | Phytophthora sojae             | P6954                                             |
| G4ZRQ8        | Phytophthora sojae             | P6954                                             |
| G5A0G9        | Phytophthora sojae             | P6497                                             |
| G5A0G9        | Phytophthora sojae             | P6497                                             |
| G5EI17        | Magnaporthe oryzae             | Ina72                                             |
| G5EI17        | Magnaporthe oryzae             | Ina72                                             |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7T9N6        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TAC5        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TAC5        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TFJ1        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| G7TFJ1        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| G7TIV7        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TIV7        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TIV7        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TJZ8        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| G7TK09        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| H9DUR1        | Verticillium dahliae           | JR2                                               |
| H9DUR1        | Verticillium dahliae           | JR2                                               |
| H9DUR1        | Verticillium dahliae           | JR2                                               |
| H9DUR1        | Verticillium dahliae           | JR2                                               |
| H9DUR1        | Verticillium dahliae           | 12008                                             |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT7        | Rhynchosporium commune         | UK7                                               |
| I1VGT8        | Rhynchosporium commune         | UK7                                               |
| I1VGT8        | Rhynchosporium commune         | UK7                                               |
| I1VGT8        | Rhynchosporium commune         | UK7                                               |
| I1VGT8        | Rhynchosporium commune         | UK7                                               |
| K4LAY3        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| K4LAY3        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| K4LAY3        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| L7NCR0        | Phytophthora capsici           | Pc537                                             |
| L7NCR0        | Phytophthora capsici           | Pc537                                             |
| L7NCS1        | Phytophthora capsici           | Pc537                                             |
| L7NCS1        | Phytophthora capsici           | Pc537                                             |
| L7WWG5        | Phytophthora sojae             | P6497                                             |
| L7WWG5        | Phytophthora sojae             | P6497                                             |
| L7WWG5        | Phytophthora sojae             | P6497                                             |
| L7WWG5        | Phytophthora sojae             | P6497                                             |
| M1T205        | Macrosiphum euphorbiae         | no data found                                     |
| M1T205        | Macrosiphum euphorbiae         | no data found                                     |
| M1T205        | Macrosiphum euphorbiae         | no data found                                     |
| N1J7E2        | Blumeria graminis              | f. sp. Hordei                                     |
| N1J7E2        | Blumeria graminis              | f. sp. hordei DH14                                |
| N1JJH4        | Blumeria graminis              | f. sp. Hordei                                     |
| N1JJH4        | Blumeria graminis              | f. sp. hordei DH14                                |
| no data found | Myzus persicae                 | RRes (genotype O)                                 |
| no data found | Myzus persicae                 | RRes (genotype O)                                 |
| no data found | Myzus persicae                 | RRes (genotype O)                                 |
| no data found | Meloidogyne javanica           | VW4                                               |
| no data found | Myzus persicae                 | RRes (genotype O)                                 |
| no data found | Myzus persicae                 | RRes (genotype O)                                 |
| no data found | Myzus persicae                 | RRes (genotype O)                                 |
| no data found | Meloidogyne javanica           | VW4                                               |
| no data found | Phytophthora parasitica        | 310                                               |
| no data found | Magnaporthe oryzae             | 70-15                                             |
| no data found | Magnaporthe oryzae             | 70-15                                             |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Melampsora larici-populina     | 98AG31                                            |
| no data found | Phytophthora capsici           | no data found                                     |
| no data found | Phytophthora capsici           | no data found                                     |
| no data found | Phytophthora capsici           | no data found                                     |
| no data found | Phytophthora capsici           | no data found                                     |
| no data found | Phytophthora capsici           | no data found                                     |
| no data found | Phytophthora capsici           | no data found                                     |
| no data found | Bremia lactucae                | no data found                                     |
| no data found | Bremia lactucae                | no data found                                     |
| no data found | Pseudomonas syringae           | pv. lachrymans M301315                            |
| no data found | Phytophthora infestans         | 88069                                             |
| no data found | Hyaloperonospora arabidopsidis | Emoy2                                             |
| no data found | Phytophthora capsici           | LT263                                             |
| no data found | Phytophthora capsici           | LT263                                             |
| no data found | Phytophthora capsici           | LT263                                             |
| no data found | Phytophthora capsici           | LT263                                             |
| no data found | Phytophthora capsici           | LT263                                             |
| no data found | Phytophthora capsici           | LT263                                             |
| no data found | Pseudocercospora fijiensis     | CIRAD86                                           |
| no data found | Pseudocercospora fijiensis     | CIRAD86                                           |
| no data found | Pseudocercospora fijiensis     | CIRAD86                                           |
| no data found | Phytophthora cactorum          | 10300                                             |
| no data found | Phytophthora cactorum          | 10300                                             |
| no data found | Ralstonia solanacearum         | GMI1000                                           |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Ralstonia solanacearum         | no data found                                     |
| no data found | Burkholderia pseudomallei      | DD503                                             |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Fusarium oxysporum             | f. sp. cepae                                      |
| no data found | Magnaporthe oryzae             | 4091-5-8                                          |
| no data found | Magnaporthe oryzae             | 4091-5-8                                          |
| no data found | Magnaporthe oryzae             | 4091-5-8                                          |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Acinetobacter nosocomialis     | M2                                                |
| no data found | Magnaporthe oryzae             | INA72                                             |
| no data found | Magnaporthe oryzae             | INA72                                             |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Puccinia striiformis           | f. sp. tritici                                    |
| no data found | Edwardsiella ictaluri          | no data found                                     |
| no data found | Staphylococcus aureus          | Geraldine (HT20030749)                            |
| no data found | Phytophthora capsici           | LT1534                                            |
| no data found | Xanthomonas axonopodis         | pv. Punicae str. ITCCBD 0003                      |
| no data found | Porphyromonas gingivalis       | ATCC 33277                                        |
| no data found | Puccinia graminis              | f. sp. tritici                                    |
| no data found | Puccinia graminis              | f. sp. tritici                                    |
| no data found | Puccinia graminis              | f. sp. tritici                                    |
| no data found | Puccinia graminis              | f. sp. tritici                                    |
| O08242        | Pseudomonas syringae           | pv. tomato PT23                                   |
| O08242        | Pseudomonas syringae           | pv. glycinea race 1                               |
| O08243        | Pseudomonas syringae           | pv. glycinea race 3                               |
| O08243        | Pseudomonas syringae           | pv. glycinea race 6                               |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium str. LT2      |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| O30916        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| O42718        | Phytophthora infestans         | 88069                                             |
| O42718        | Phytophthora infestans         | 88069                                             |
| O42719        | Phytophthora infestans         | 88069                                             |
| O42719        | Phytophthora infestans         | 88069                                             |
| O52623        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| O52623        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| O52623        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| P0CL52        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| P11437        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P11437        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P11437        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P11437        | Pseudomonas syringae           | pv. glycinea race 6                               |
| P11437        | Pseudomonas syringae           | pv. glycinea race 6                               |
| P11437        | Pseudomonas syringae           | pv. glycinea race 6                               |
| P11437        | Pseudomonas syringae           | pv. glycinea race 6                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. glycinea race 4                               |
| P13835        | Pseudomonas syringae           | pv. syringae                                      |
| P17778        | Yersinia pestis                | KIM1001                                           |
| P17778        | Yersinia pestis                | KIM1001                                           |
| P17778        | Yersinia pestis                | 91001                                             |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | Race 5                                            |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | Race 0                                            |
| P22287        | Passalora fulva                | Race 2                                            |
| P22287        | Passalora fulva                | Race 4                                            |
| P22287        | Passalora fulva                | Race 5                                            |
| P22287        | Passalora fulva                | Race 9                                            |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| P22287        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | Multiple alleles from numerous strains            |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | race 4                                            |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00363        | Passalora fulva                | no data found                                     |
| Q00365        | Passalora fulva                | race 5                                            |
| Q00365        | Passalora fulva                | no data found                                     |
| Q00365        | Passalora fulva                | no data found                                     |
| Q01144        | Magnaporthe oryzae             | 4392-1-6                                          |
| Q01144        | Magnaporthe oryzae             | CP987                                             |
| Q01144        | Magnaporthe oryzae             | 4091-5-8                                          |
| Q01905        | Phytophthora infestans         | 88069                                             |
| Q01905        | Phytophthora infestans         | 88069                                             |
| Q01905        | Phytophthora infestans         | 88069                                             |
| Q02039        | Rhynchosporium commune         | no data found                                     |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q02039        | Rhynchosporium commune         | UK7                                               |
| Q08242        | Pseudomonas syringae           | pv. tomato JL1065                                 |
| Q08242        | Pseudomonas syringae           | no data found                                     |
| Q08242        | Pseudomonas syringae           | no data found                                     |
| Q258K5        | Leptosphaeria maculans         | v23.1.3                                           |
| Q258K5        | Leptosphaeria maculans         | no data found                                     |
| Q258K5        | Leptosphaeria maculans         | no data found                                     |
| Q2A0P1        | Fusarium oxysporum             | f. sp. Lycopersici 4287                           |
| Q2A0P1        | Fusarium oxysporum             | f. sp. lycopersici Fol007                         |
| Q2A0P1        | Verticillium dahliae           | JR2                                               |
| Q2A0P1        | Verticillium dahliae           | JR2                                               |
| Q2A0P1        | Botrytis cinerea               | B05.10                                            |
| Q2A0P1        | Botrytis cinerea               | B05.10                                            |
| Q2A0P1        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q2A0P1        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q2A0P1        | Fusarium oxysporum             | f. sp. Lycopersici, strain Fol007                 |
| Q2A0P1        | Fusarium oxysporum             | Fo5176                                            |
| Q3BM44        | Xanthomonas campestris         | pv. vesicatoria str. 85-10                        |
| Q3BM44        | Xanthomonas campestris         | pv. vesicatoria str. 85-10                        |
| Q3BM44        | Xanthomonas campestris         | pv. vesicatoria str. 85-10                        |
| Q3BM44        | Xanthomonas campestris         | pv. vesicatoria str. 85-10                        |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B66        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B68        | Pseudomonas syringae           | pv. phaseolicola 1448A                            |
| Q48B68        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B68        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B68        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B68        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q48B68        | Pseudomonas syringae           | pv. glycinea race 4                               |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4PET8        | Ustilago maydis                | SG200                                             |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris B186                               |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris B186                               |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris B186                               |
| Q4UWF4        | Xanthomonas campestris         | pv. Campestris B186                               |
| Q4VKJ1        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ1        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ1        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Maks9                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Noks1                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emwa1                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Maks9                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Noks1                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emwa1                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Maks9                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emwa1                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Maks9                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emwa1                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Maks9                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Cala2                                             |
| Q4VKJ6        | Hyaloperonospora arabidopsidis | Emwa1                                             |
| Q4ZMD6        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZMD6        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX47        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX47        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX47        | Pseudomonas syringae           | no data found                                     |
| Q4ZX47        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q4ZX49        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX49        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q4ZX82        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX82        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX85        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZX85        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q4ZYH0        | Pseudomonas syringae           | pv. maculicola                                    |
| Q4ZYH0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q4ZYH0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q52389        | Pseudomonas syringae           | pv. tomato PT23                                   |
| Q52389        | Pseudomonas syringae           | pv. tomato PT23                                   |
| Q52432        | Pseudomonas syringae           | pv. pisi 151                                      |
| Q52432        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q52432        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q572D3        | Phytophthora infestans         | no data found                                     |
| Q572D3        | Phytophthora infestans         | HP10-45                                           |
| Q572D3        | Phytophthora infestans         | 88069                                             |
| Q572D3        | Phytophthora infestans         | 88069                                             |
| Q572D3        | Phytophthora infestans         | 88069                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emco5                                             |
| Q5G7K8        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q5G7L2        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q5G7L2        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q5G7L2        | Hyaloperonospora arabidopsidis | Emoy2                                             |
| Q5TIP4        | Botrytis cinerea               | no data found                                     |
| Q5TIP4        | Botrytis cinerea               | no data found                                     |
| Q5X811        | Legionella pneumophila         | 130b                                              |
| Q5X811        | Legionella pneumophila         | 130b                                              |
| Q63K35        | Burkholderia pseudomallei      | DD503                                             |
| Q63K35        | Burkholderia pseudomallei      | 10276                                             |
| Q6EES5        | Pseudomonas syringae           | pv. syringae 61                                   |
| Q6EES5        | Pseudomonas syringae           | pv. syringae 61                                   |
| Q6LAD6        | Pseudomonas syringae           | pv. tomato JL1065                                 |
| Q6LAD6        | Pseudomonas syringae           | pv. tomato JL1065                                 |
| Q6LAD6        | Pseudomonas syringae           | pv. tomato JL1065                                 |
| Q6LAD6        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q6R659        | Melampsora lini                | no data found                                     |
| Q6R659        | Melampsora lini                | no data found                                     |
| Q6R661        | Melampsora lini                | no data found                                     |
| Q6R661        | Melampsora lini                | no data found                                     |
| Q6TKR8        | Xanthomonas oryzae             | pv. Oryzicola BLS256                              |
| Q6TKR8        | Xanthomonas oryzae             | pv. Oryzicola MAI10                               |
| Q6TKR8        | Xanthomonas oryzae             | X11-5A                                            |
| Q6XVY7        | Shigella flexneri              | M90T-Sm                                           |
| Q6XVY7        | Shigella flexneri              | M90T-Sm                                           |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q79LY0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7CQD4        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q7DB85        | Escherichia coli               | EDL933                                            |
| Q7DB85        | Escherichia coli               | EDL933                                            |
| Q7DB85        | Escherichia coli               | EDL933                                            |
| Q7PC62        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q7PC62        | Pseudomonas syringae           | pv. syringae B728a                                |
| Q83B67        | Coxiella burnetii              | Nine Mile Phase II Clone 4 (RSA439)               |
| Q83B67        | Coxiella burnetii              | Nine Mile Phase II Clone 4 (RSA439)               |
| Q83FB9        | Coxiella burnetii              | NMII                                              |
| Q83FB9        | Coxiella burnetii              | Tn1832                                            |
| Q83FB9        | Coxiella burnetii              | Nine Mile (RSA493)                                |
| Q87W07        | Magnaporthe oryzae             | Ku80                                              |
| Q87W07        | Magnaporthe oryzae             | Ku80                                              |
| Q87W42        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q87W42        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q887D0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q887D0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q887D0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q887D0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q888W0        | Magnaporthe oryzae             | Ku80                                              |
| Q888W0        | Magnaporthe oryzae             | Ku80                                              |
| Q888W0        | Magnaporthe oryzae             | 70-15                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888W0        | Magnaporthe oryzae             | KV103                                             |
| Q888Y7        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q888Y7        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q888Y7        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q888Y7        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q888Y7        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q88BH0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q88BH0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q88BH0        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q88BW6        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q88BW6        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q8H6Z4        | Phytophthora infestans         | 90128                                             |
| Q8H6Z4        | Phytophthora infestans         | 90128                                             |
| Q8H6Z4        | Phytophthora infestans         | 90128                                             |
| Q8H6Z6        | Phytophthora infestans         | 90128                                             |
| Q8H6Z6        | Phytophthora infestans         | 90128                                             |
| Q8H6Z6        | Phytophthora infestans         | 90128                                             |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | race 5                                            |
| Q8NID8        | Passalora fulva                | race 5                                            |
| Q8NID8        | Passalora fulva                | race 2                                            |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8NID8        | Passalora fulva                | no data found                                     |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PC98        | Xanthomonas campestris         | pv. Campestris str. 8004                          |
| Q8PHM7        | Xanthomonas axonopodis         | pv. citri                                         |
| Q8PHM7        | Xanthomonas axonopodis         | pv. citri                                         |
| Q8PQN7        | Xanthomonas axonopodis         | pv. citri                                         |
| Q8PQN7        | Xanthomonas axonopodis         | pv. citri                                         |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRG7        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRK7        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRM3        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. Citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8PRN6        | Xanthomonas citri              | subsp. citri 306                                  |
| Q8RP09        | Pseudomonas syringae           | pv. maculicola str. ES4326                        |
| Q8RP09        | Pseudomonas syringae           | pv. maculicola str. ES4326                        |
| Q8RSY1        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q8RSY1        | Pseudomonas syringae           | no data found                                     |
| Q8RSY1        | Pseudomonas syringae           | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XPQ6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQ26        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQA2        | Ralstonia solanacearum         | no data found                                     |
| Q8XQE6        | Ralstonia solanacearum         | UW551                                             |
| Q8XQE6        | Ralstonia solanacearum         | UW551                                             |
| Q8XQE6        | Ralstonia solanacearum         | UW551                                             |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQI6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK6        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQK7        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQL0        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XQT8        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | no data found                                     |
| Q8XR46        | Ralstonia solanacearum         | GMI1000                                           |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRE0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI0        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRI4        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XRK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XSA6        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT13        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT19        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT97        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT98        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XT99        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | no data found                                     |
| Q8XTA1        | Ralstonia solanacearum         | GMI1000                                           |
| Q8XTA1        | Ralstonia solanacearum         | GMI1000                                           |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTF2        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTK9        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XTS6        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XU25        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUH6        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XUL4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVD4        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XVQ5        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XWW1        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXI2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XXL2        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYB9        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | GMI1000                                           |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYE3        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF7        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYF8        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XYN5        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN8        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8XZN9        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y0Z8        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | no data found                                     |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y125        | Ralstonia solanacearum         | GMI1000                                           |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y164        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8Y2L4        | Ralstonia solanacearum         | no data found                                     |
| Q8YAR0        | Listeria monocytogenes         | EGD-e                                             |
| Q8YAR0        | Listeria monocytogenes         | EGD-e                                             |
| Q8ZNR3        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q8ZNR3        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q8ZNR3        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q8ZNR3        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain 14028  |
| Q8ZPY9        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| Q8ZPY9        | Salmonella enterica            | subsp. enterica serovar Typhimurium strain SL1344 |
| Q94FS7        | Phytophthora cactorum          | P381                                              |
| Q94FS7        | Phytophthora cactorum          | P381                                              |
| Q9AJW7        | Shigella flexneri              | M90T-Sm                                           |
| Q9AJW7        | Shigella flexneri              | M90T-Sm                                           |
| Q9AT28        | Phytophthora parasitica        | 1828                                              |
| Q9AT28        | Phytophthora parasitica        | 1828                                              |
| Q9AT28        | Phytophthora parasitica        | 1828                                              |
| Q9C478        | Magnaporthe oryzae             | PWL2                                              |
| Q9C478        | Magnaporthe oryzae             | CP987                                             |
| Q9C478        | Magnaporthe oryzae             | 4360-R-62                                         |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzae X11-5A                                 |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9FDG1        | Xanthomonas oryzae             | pv. oryzicola BLS256                              |
| Q9JP38        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q9JP38        | Pseudomonas syringae           | pv. Tomato str DC3000                             |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBS0        | Ralstonia solanacearum         | no data found                                     |
| Q9RBW3        | Pseudomonas syringae           | pv. phaseolicola 1449B (race7)                    |
| Q9RBW3        | Pseudomonas syringae           | pv. phaseolicola 1449B (race7)                    |
| Q9SPD4        | Pythium aphanidermatum         | no data found                                     |
| Q9SPD4        | Pythium aphanidermatum         | no data found                                     |
| Q9SPD4        | Pythium aphanidermatum         | no data found                                     |
| R4TJS0        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| R4TJS0        | Xanthomonas oryzae             | pv. oryzae BXO1                                   |
| R9RX08        | Xanthomonas oryzae             | pv. oryzae PXO99                                  |
| R9RX08        | Xanthomonas oryzae             | pv. oryzae PXO99                                  |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |
| W0H866        | Pseudomonas cichorii           | JBC1                                              |

1.  Pathogen species that have same protein IDs

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
effector_data <- data.table::fread("../../data/effector_with_IDs_organism.csv") %>% 
  rename(`Protein ID` = protein_id, `Pathogen species` = pathogen_short)
```

Now, we can map the effector data with all of the strains via protein
IDs

``` r
effector_sequence_with_metadata <- effector_data %>% 
  left_join(., phi_plant_effector, by = 'Protein ID') %>% 
  select(-c('name_src', 'Pathogen species.x')) %>% 
  rename(`Pathogen species` = `Pathogen species.y`)

effector_sequence_with_metadata <- effector_sequence_with_metadata %>% 
  group_by(`Pathogen strain`, `Pathogen species`)

write.csv(effector_sequence_with_metadata , "../../data/effector_sequence_with_metadata.csv", row.names = FALSE)
```

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
