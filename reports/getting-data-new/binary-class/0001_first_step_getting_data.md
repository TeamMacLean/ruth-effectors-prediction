Getting and learning the data
=============================

Introduction
------------

In this report, it shows the process on getting effector and non
effector from the beginning (since the previous effector and noneffector
data were incorrectly retrieved).

Getting the Phi-base IDs from the Phi-base data
-----------------------------------------------

``` r
# Read the last version of Phi-base data (CSV format)
phi_base <- data.table::fread("../../../../data/phi-base-current-data/phi-base_without_column_desc.csv", header = TRUE)
```

In order to get the effector data, we need to filter all of the data
with `GeneFuction == "effector"`

``` r
# Filter al of the effector Gene Function
# Use `ignore_case = TRUE` in order to ignore whether they are lower or uppercase
phi_effector <- phi_base %>%
  dplyr::filter(str_detect(GeneFunction, regex("effector", ignore_case = TRUE)))
```

We need to make sure whether all of the gene function filtered is
effector, by find the unique function:

``` r
# Make sure there is no `NA` column so that we do not miss something
phi_effector  %>% 
  select(GeneFunction) %>% 
  is.na() %>% 
  which()
```

    ## integer(0)

``` r
# Check the unique GeneFunction
phi_effector %>% 
  select(GeneFunction) %>% 
  unique() %>% 
  knitr::kable()
```

|      | GeneFunction                                                                  |
|------|:------------------------------------------------------------------------------|
| 1    | Effector protein                                                              |
| 11   | Effector protein, polyketide synthase                                         |
| 37   | Type III effector                                                             |
| 76   | effector                                                                      |
| 131  | encode secreted effectors                                                     |
| 154  | Candidate Effector genes                                                      |
| 162  | Translocation of T6SS effector proteins essential for virulence               |
| 442  | RxLR effector                                                                 |
| 460  | Tool for translocating pathogen effectors to monocot cells                    |
| 468  | effector protein                                                              |
| 499  | TAL effectors                                                                 |
| 501  | Type III effector protein                                                     |
| 502  | RXLR effector                                                                 |
| 503  | Effector                                                                      |
| 509  | Type III Effector                                                             |
| 514  | phospholipase A2 effector Exotoxin                                            |
| 523  | Effector (plant avirulence determinant)                                       |
| 532  | Part of the largest Effector Gene Cluster                                     |
| 564  | Cyclic di-GMP Effectors                                                       |
| 566  | RxLR Effector                                                                 |
| 570  | TTSS effector protein                                                         |
| 612  | AvrE-family Type III Effector                                                 |
| 641  | Type VI secretion systems (T6SSs) effectors                                   |
| 645  | Secreted effector protein                                                     |
| 692  | type III effector                                                             |
| 693  | Type 3 Secretion System Effectors                                             |
| 705  | CRN effectors                                                                 |
| 711  | GKLR effectors                                                                |
| 714  | T6SS effector                                                                 |
| 716  | Ysa T3SS are Ysc effectors                                                    |
| 717  | T3SS2 effector                                                                |
| 719  | Dot/Icm Effector                                                              |
| 722  | T3S effector                                                                  |
| 725  | AvrE-Family Effector                                                          |
| 728  | secreted effector                                                             |
| 730  | modulates type III effector secretion                                         |
| 733  | encodes type III effector                                                     |
| 735  | YopJ Superfamily Effector                                                     |
| 736  | encodes an a-1,3-mannosyltransferase for protein N-glycosylation of effectors |
| 738  | putative RXLR effector genes                                                  |
| 745  | effector lipase                                                               |
| 746  | Virulence effector                                                            |
| 757  | Putative Effector protein                                                     |
| 848  | Effector virulence factor                                                     |
| 849  | Regulates the expression on Pgp3, the effector virulence factor               |
| 850  | Serine Protease Effectors                                                     |
| 853  | Type III type effector                                                        |
| 863  | SCR effector                                                                  |
| 867  | RXLR Effector                                                                 |
| 870  | T6SS Effector                                                                 |
| 879  | Pathogenicity island 1 effector protein                                       |
| 882  | Type III effectors (T3E) and T3E-like genes                                   |
| 1439 | T3SS secreted effector                                                        |
| 1446 | T3SS3 effector                                                                |
| 1458 | T6SS5 Effector                                                                |
| 1496 | CRN effector                                                                  |
| 1609 | T2SS effector                                                                 |
| 1680 | Type III secretion system effector                                            |
| 1736 | Candidate Secreted Effector Proteins                                          |
| 1742 | Putative effector protein                                                     |
| 1744 | Acetyltransferase effector                                                    |
| 1755 | iT3SS: translocating effectors in the host cells                              |
| 1770 | LysM effector                                                                 |
| 1783 | H2-T6SS-dependent phospholipase effector                                      |
| 1784 | transcription activation-like (TAL) effector gene                             |
| 1812 | T3SS effector protein                                                         |
| 1821 | Type IVb effector                                                             |
| 1823 | T6SS Effector protein                                                         |
| 1832 | Type III Effector Tyrosine Phosphatase                                        |
| 1833 | TIR effector                                                                  |
| 1835 | TAL effector                                                                  |
| 1840 | E3 ubiquitin ligase effector                                                  |
| 1844 | Effector type III secreted protein                                            |
| 1848 | Type III secreted effector                                                    |
| 1853 | Effector protein, phosphothreonine lyase                                      |
| 1866 | LysM effector protein                                                         |
| 1868 | T3SS Effector protein                                                         |
| 1870 | Putative type VI secretion system, effector protein                           |
| 1940 | Transcription activator-like (TAL) effector protein                           |
| 1947 | Type III secretion system 2 effector                                          |
| 2069 | SPI-2 effector protein                                                        |
| 2071 | Necrotrophic effector protein                                                 |

Note that, although the gene function above contains the word
“effector”, but not all are actually effector. Then now we can make a
list of the description above, which is not effector but it contains
string `effector`.

``` r
# Create a list of GeneFunction that are not effectors
excluded_GeneFunction <- c("encode secreted effectors", 
                           "Translocation of T6SS effector proteins essential for virulence",
                           "Tool for translocating pathogen effectors to monocot cells", 
                           "encodes an a-1,3-mannosyltransferase for protein N-glycosylation of effectors", 
                           "modulates type III effector secretion", 
                           "encodes type III effector", 
                           "iT3SS: translocating effectors in the host cells")
```

Next step is to excluded any data with the gene function as in the list
above:

``` r
phi_effector <- phi_effector %>% 
  dplyr::filter(., !GeneFunction %in% excluded_GeneFunction)

phi_effector %>%  
  nrow()
```

    ## [1] 2219

Now, since we will focus more on plant effector, we need also to exclude
all of the effector data with non-plant host. Now we can see what are
the host list in our data, which apparently, not only plants:

``` r
# Make sure there is no `NA` column so that we do not miss something
phi_effector  %>% 
  select(Hostdescription) %>% 
  is.na() %>% 
  which()
```

    ## integer(0)

Since there is no missing data, then now we can proceed to the next
step.

``` r
phi_effector %>% 
  select(Hostdescription) %>% 
  unique() %>% 
  knitr::kable()
```

|      | Hostdescription      |
|------|:---------------------|
| 1    | Eudicots             |
| 4    | Monocots             |
| 502  | Rodents              |
| 564  | Cellular slime molds |
| 635  | Birds                |
| 710  | Flies                |
| 711  | Rabbits & hares      |
| 714  | Moths                |
| 1698 | Bony fishes          |
| 1769 | Primates             |
| 1819 | Nematodes            |

Short definition of host description:

1.  *Eudicots*: or Eudicotidae or eudicotyledons are a clade of
    flowering plants (Wikipedia)
2.  *Monocots*: commonly referred to as monocots, are flowering plants
    the seeds of which typically contain only one embryonic leaf, or
    cotyledon.

For the data that only have hosts as Eudicots or Monocots

``` r
# Restrict the host only for plants host
host <- c("Eudicots", "Monocots")

# Get the effector data only for the host "Eudicots" and "Monocots"
phi_effector_plant <- phi_effector %>% 
  dplyr::filter(., Hostdescription %in% host)

# Get all of effector data with host which are not "Eudicots" and "Monocots"
phi_effector_not_plant <- phi_effector %>% 
  dplyr::filter(., !(Hostdescription %in% host))

# Get HostID for non plant effector
phi_HostID_non_plant <- phi_effector_not_plant[["HostID"]] %>% unique()

# Get HostID for plant effector
phi_HostID_plant <- phi_effector_plant[["HostID"]] %>% unique()

print(phi_HostID_non_plant)
```

    ##  [1] "10090" "44689" "10141" "8839"  "7227"  "9986"  "7137"  "7998"  "9031" 
    ## [10] "9606"  "6239"  "10036"

``` r
print(phi_HostID_plant)
```

    ##  [1] "4081"   "3708"   "4513"   "4530"   "4113"   "3702"   "4100"   "3847"  
    ##  [9] "4006"   "4097"   "3747"   "4043"   "856898" "73574"  "3885"   "3888"  
    ## [17] "2711"   "3983"   "4577"   "4072"   "3659"   "4565"   "3712"   "4073"  
    ## [25] "3656"   "4111"   "3726"   "4236"   "29760"  "3880"   "4679"   "2708"  
    ## [33] "62141"  "37656"  "23211"  "3634"   "4146"   "3635"   "3750"   "22663" 
    ## [41] "42229"  "4571"

### Making sure that the effector data are plant effectors

In order to make sure that our out plant effector data have plant as the
host (“Viridiplantae”), then we can use the `lib(taxize)` which is a
library in R to check the taxonomy details from NCBI given taxonomy IDs.

``` r
# Check for the list ID of plant effector
ncbi_list_plant <- taxize::classification(phi_HostID_plant, db = "ncbi")

df_phi_HostID_plant <- data.frame(HostID = phi_HostID_plant)

df_phi_HostID_plant_check <- df_phi_HostID_plant %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    `Is Plant` = "Viridiplantae" %in% ncbi_list_plant[[HostID]]$name
  )

df_phi_HostID_plant_check %>% 
  knitr::kable()
```

| HostID | Is Plant |
|:-------|:---------|
| 4081   | TRUE     |
| 3708   | TRUE     |
| 4513   | TRUE     |
| 4530   | TRUE     |
| 4113   | TRUE     |
| 3702   | TRUE     |
| 4100   | TRUE     |
| 3847   | TRUE     |
| 4006   | TRUE     |
| 4097   | TRUE     |
| 3747   | TRUE     |
| 4043   | TRUE     |
| 856898 | TRUE     |
| 73574  | TRUE     |
| 3885   | TRUE     |
| 3888   | TRUE     |
| 2711   | TRUE     |
| 3983   | TRUE     |
| 4577   | TRUE     |
| 4072   | TRUE     |
| 3659   | TRUE     |
| 4565   | TRUE     |
| 3712   | TRUE     |
| 4073   | TRUE     |
| 3656   | TRUE     |
| 4111   | TRUE     |
| 3726   | TRUE     |
| 4236   | TRUE     |
| 29760  | TRUE     |
| 3880   | TRUE     |
| 4679   | TRUE     |
| 2708   | TRUE     |
| 62141  | TRUE     |
| 37656  | TRUE     |
| 23211  | TRUE     |
| 3634   | TRUE     |
| 4146   | TRUE     |
| 3635   | TRUE     |
| 3750   | TRUE     |
| 22663  | TRUE     |
| 42229  | TRUE     |
| 4571   | TRUE     |

The results above shows that all of the data that are already
categorized based the host “Eudicots” and “Monocots” are effectors on
plants.

``` r
# Check for the list ID of non plant effector
ncbi_list_non_plant <- taxize::classification(phi_HostID_non_plant, db = "ncbi")

df_phi_HostID_non_plant <- data.frame(HostID = phi_HostID_non_plant)

df_phi_HostID_non_plant_check <- df_phi_HostID_non_plant %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    `Is Plant` = "Viridiplantae" %in% ncbi_list_non_plant[[HostID]]$name
  )

df_phi_HostID_non_plant_check %>% 
  knitr::kable()
```

| HostID | Is Plant |
|:-------|:---------|
| 10090  | FALSE    |
| 44689  | FALSE    |
| 10141  | FALSE    |
| 8839   | FALSE    |
| 7227   | FALSE    |
| 9986   | FALSE    |
| 7137   | FALSE    |
| 7998   | FALSE    |
| 9031   | FALSE    |
| 9606   | FALSE    |
| 6239   | FALSE    |
| 10036  | FALSE    |

The checking process above, using the package `taxize` is very useful to
make sure that all the data we picked have plant as host.

### Selecting the relevant columns of the data

``` r
# Make list of all of the column that are needed
select_columns <- c("PHIMolConnID",
"ProteinIDsource",
"ProteinID",
"PathogenID",   
"Pathogenspecies",  
"PathogenstrainID",
"Pathogenstrain",   
"Hostdescription",
"HostID",
"Hostspecies",
"Hoststrain",
"GeneFunction",
"MutantPhenotype") 

# Select from the data phi_effector_plant
phi_effector_plant <- phi_effector_plant %>% 
  dplyr::select(select_columns)
```

Getting know the data
---------------------

After getting all of the data, we need to understand the data itself.
Since the Phi-base CSV data does not contain any sequence, then after
getting all of the IDs of effector data, we need to retrieve the
sequences from other resources. Now we have to take a look where are the
protein sequence sources.

### Protein source data

``` r
# Check the list of Protein sequence sources
phi_effector_plant %>% 
  group_by(ProteinIDsource) %>% 
  summarise(count = n()) %>% 
  knitr::kable()
```

| ProteinIDsource |  count|
|:----------------|------:|
|                 |    277|
| Uniprot         |   1714|

The proteinID sources are mainly from Uniprot, where around 16% of the
rest the data are unknown.

``` r
phi_effector_plant %>% 
  filter(ProteinIDsource == "")
```

    ##     PHIMolConnID ProteinIDsource ProteinID PathogenID
    ## 1       PHI:2577                                13164
    ## 2       PHI:2578                                13164
    ## 3       PHI:2579                                13164
    ## 4       PHI:2591                                 6303
    ## 5       PHI:2598                                13164
    ## 6       PHI:2599                                13164
    ## 7       PHI:2600                                13164
    ## 8       PHI:2612                                 6303
    ## 9       PHI:2804                                93612
    ## 10      PHI:2804                                93612
    ## 11      PHI:2942                                 4792
    ## 12      PHI:3217                               318829
    ## 13      PHI:3217                               318829
    ## 14      PHI:4527                               203908
    ## 15      PHI:4528                               203908
    ## 16      PHI:4529                               203908
    ## 17      PHI:4530                               203908
    ## 18      PHI:4531                               203908
    ## 19      PHI:4532                               203908
    ## 20      PHI:4533                               203908
    ## 21      PHI:4534                               203908
    ## 22      PHI:4535                               203908
    ## 23      PHI:4536                               203908
    ## 24      PHI:4537                               203908
    ## 25      PHI:4538                               203908
    ## 26      PHI:4539                               203908
    ## 27      PHI:4540                               203908
    ## 28      PHI:4541                               203908
    ## 29      PHI:4542                               203908
    ## 30      PHI:4543                               203908
    ## 31      PHI:4544                               203908
    ## 32      PHI:4545                               203908
    ## 33      PHI:4546                               203908
    ## 34      PHI:4547                               203908
    ## 35      PHI:4548                               203908
    ## 36      PHI:4549                               203908
    ## 37      PHI:4550                               203908
    ## 38      PHI:3938                                 4784
    ## 39      PHI:3939                                 4784
    ## 40      PHI:3940                                 4784
    ## 41      PHI:3938                                 4784
    ## 42      PHI:3939                                 4784
    ## 43      PHI:3940                                 4784
    ## 44      PHI:3969                                 4779
    ## 45      PHI:3970                                 4779
    ## 46      PHI:4153                                  317
    ## 47      PHI:4682                                 4787
    ## 48      PHI:4760                               272952
    ## 49      PHI:4838                                 4784
    ## 50      PHI:4838                                 4784
    ## 51      PHI:4838                                 4784
    ## 52      PHI:4838                                 4784
    ## 53      PHI:4838                                 4784
    ## 54      PHI:4838                                 4784
    ## 55      PHI:2435                              1873960
    ## 56      PHI:2435                              1873960
    ## 57      PHI:2436                              1873960
    ## 58      PHI:5032                                29920
    ## 59      PHI:5032                                29920
    ## 60      PHI:5070                                  305
    ## 61      PHI:5146                                  305
    ## 62      PHI:5147                                  305
    ## 63      PHI:5148                                  305
    ## 64      PHI:5149                                  305
    ## 65      PHI:5150                                  305
    ## 66      PHI:5151                                  305
    ## 67      PHI:5152                                  305
    ## 68      PHI:5153                                  305
    ## 69      PHI:5154                                  305
    ## 70      PHI:5155                                  305
    ## 71      PHI:5156                                  305
    ## 72      PHI:5157                                  305
    ## 73      PHI:5161                                  305
    ## 74      PHI:5167                                  305
    ## 75      PHI:5168                                  305
    ## 76      PHI:5170                                  305
    ## 77      PHI:5176                                  305
    ## 78      PHI:5180                                  305
    ## 79      PHI:5182                                  305
    ## 80      PHI:5146                                  305
    ## 81      PHI:5147                                  305
    ## 82      PHI:5148                                  305
    ## 83      PHI:5149                                  305
    ## 84      PHI:5150                                  305
    ## 85      PHI:5151                                  305
    ## 86      PHI:5152                                  305
    ## 87      PHI:5153                                  305
    ## 88      PHI:5154                                  305
    ## 89      PHI:5155                                  305
    ## 90      PHI:5156                                  305
    ## 91      PHI:5157                                  305
    ## 92      PHI:5161                                  305
    ## 93      PHI:5167                                  305
    ## 94      PHI:5168                                  305
    ## 95      PHI:5170                                  305
    ## 96      PHI:5176                                  305
    ## 97      PHI:5180                                  305
    ## 98      PHI:5182                                  305
    ## 99      PHI:5146                                  305
    ## 100     PHI:5147                                  305
    ## 101     PHI:5148                                  305
    ## 102     PHI:5149                                  305
    ## 103     PHI:5150                                  305
    ## 104     PHI:5151                                  305
    ## 105     PHI:5152                                  305
    ## 106     PHI:5153                                  305
    ## 107     PHI:5154                                  305
    ## 108     PHI:5155                                  305
    ## 109     PHI:5156                                  305
    ## 110     PHI:5157                                  305
    ## 111     PHI:5161                                  305
    ## 112     PHI:5167                                  305
    ## 113     PHI:5168                                  305
    ## 114     PHI:5170                                  305
    ## 115     PHI:5176                                  305
    ## 116     PHI:5180                                  305
    ## 117     PHI:5182                                  305
    ## 118     PHI:5146                                  305
    ## 119     PHI:5147                                  305
    ## 120     PHI:5148                                  305
    ## 121     PHI:5149                                  305
    ## 122     PHI:5150                                  305
    ## 123     PHI:5151                                  305
    ## 124     PHI:5152                                  305
    ## 125     PHI:5153                                  305
    ## 126     PHI:5154                                  305
    ## 127     PHI:5155                                  305
    ## 128     PHI:5156                                  305
    ## 129     PHI:5157                                  305
    ## 130     PHI:5161                                  305
    ## 131     PHI:5167                                  305
    ## 132     PHI:5168                                  305
    ## 133     PHI:5170                                  305
    ## 134     PHI:5176                                  305
    ## 135     PHI:5180                                  305
    ## 136     PHI:5182                                  305
    ## 137     PHI:5146                                  305
    ## 138     PHI:5147                                  305
    ## 139     PHI:5148                                  305
    ## 140     PHI:5149                                  305
    ## 141     PHI:5150                                  305
    ## 142     PHI:5151                                  305
    ## 143     PHI:5152                                  305
    ## 144     PHI:5153                                  305
    ## 145     PHI:5154                                  305
    ## 146     PHI:5155                                  305
    ## 147     PHI:5156                                  305
    ## 148     PHI:5157                                  305
    ## 149     PHI:5161                                  305
    ## 150     PHI:5167                                  305
    ## 151     PHI:5168                                  305
    ## 152     PHI:5170                                  305
    ## 153     PHI:5176                                  305
    ## 154     PHI:5180                                  305
    ## 155     PHI:5182                                  305
    ## 156     PHI:5146                                  305
    ## 157     PHI:5147                                  305
    ## 158     PHI:5148                                  305
    ## 159     PHI:5149                                  305
    ## 160     PHI:5150                                  305
    ## 161     PHI:5151                                  305
    ## 162     PHI:5152                                  305
    ## 163     PHI:5153                                  305
    ## 164     PHI:5154                                  305
    ## 165     PHI:5155                                  305
    ## 166     PHI:5156                                  305
    ## 167     PHI:5157                                  305
    ## 168     PHI:5161                                  305
    ## 169     PHI:5167                                  305
    ## 170     PHI:5168                                  305
    ## 171     PHI:5170                                  305
    ## 172     PHI:5176                                  305
    ## 173     PHI:5180                                  305
    ## 174     PHI:5182                                  305
    ## 175     PHI:5146                                  305
    ## 176     PHI:5147                                  305
    ## 177     PHI:5148                                  305
    ## 178     PHI:5149                                  305
    ## 179     PHI:5150                                  305
    ## 180     PHI:5151                                  305
    ## 181     PHI:5152                                  305
    ## 182     PHI:5153                                  305
    ## 183     PHI:5154                                  305
    ## 184     PHI:5155                                  305
    ## 185     PHI:5156                                  305
    ## 186     PHI:5157                                  305
    ## 187     PHI:5161                                  305
    ## 188     PHI:5167                                  305
    ## 189     PHI:5168                                  305
    ## 190     PHI:5170                                  305
    ## 191     PHI:5176                                  305
    ## 192     PHI:5180                                  305
    ## 193     PHI:5182                                  305
    ## 194     PHI:5146                                  305
    ## 195     PHI:5147                                  305
    ## 196     PHI:5148                                  305
    ## 197     PHI:5149                                  305
    ## 198     PHI:5150                                  305
    ## 199     PHI:5151                                  305
    ## 200     PHI:5152                                  305
    ## 201     PHI:5153                                  305
    ## 202     PHI:5154                                  305
    ## 203     PHI:5155                                  305
    ## 204     PHI:5156                                  305
    ## 205     PHI:5157                                  305
    ## 206     PHI:5161                                  305
    ## 207     PHI:5167                                  305
    ## 208     PHI:5168                                  305
    ## 209     PHI:5170                                  305
    ## 210     PHI:5176                                  305
    ## 211     PHI:5180                                  305
    ## 212     PHI:5182                                  305
    ## 213     PHI:5146                                  305
    ## 214     PHI:5147                                  305
    ## 215     PHI:5148                                  305
    ## 216     PHI:5149                                  305
    ## 217     PHI:5150                                  305
    ## 218     PHI:5151                                  305
    ## 219     PHI:5152                                  305
    ## 220     PHI:5153                                  305
    ## 221     PHI:5154                                  305
    ## 222     PHI:5155                                  305
    ## 223     PHI:5156                                  305
    ## 224     PHI:5157                                  305
    ## 225     PHI:5161                                  305
    ## 226     PHI:5167                                  305
    ## 227     PHI:5168                                  305
    ## 228     PHI:5170                                  305
    ## 229     PHI:5176                                  305
    ## 230     PHI:5180                                  305
    ## 231     PHI:5182                                  305
    ## 232     PHI:5355                                 5507
    ## 233     PHI:5356                                 5507
    ## 234     PHI:5357                                 5507
    ## 235     PHI:5358                                 5507
    ## 236     PHI:5359                                 5507
    ## 237     PHI:5360                                 5507
    ## 238     PHI:5361                                 5507
    ## 239     PHI:5362                                 5507
    ## 240     PHI:5363                                 5507
    ## 241     PHI:5364                                 5507
    ## 242     PHI:5488                               318829
    ## 243     PHI:5489                               318829
    ## 244     PHI:5490                               318829
    ## 245     PHI:6082                               318829
    ## 246     PHI:6083                               318829
    ## 247     PHI:6160                                27350
    ## 248     PHI:6161                                27350
    ## 249     PHI:6162                                27350
    ## 250     PHI:6163                                27350
    ## 251     PHI:6164                                27350
    ## 252     PHI:6165                                27350
    ## 253     PHI:6166                                27350
    ## 254     PHI:6168                                27350
    ## 255     PHI:6169                                27350
    ## 256     PHI:6170                                27350
    ## 257     PHI:6171                                27350
    ## 258     PHI:6172                                27350
    ## 259     PHI:6173                                27350
    ## 260     PHI:6174                                27350
    ## 261     PHI:6175                                27350
    ## 262     PHI:7313                                 4784
    ## 263     PHI:7736                                53413
    ## 264     PHI:7989                                 5297
    ## 265     PHI:7989                                 5297
    ## 266     PHI:7989                                 5297
    ## 267     PHI:7989                                 5297
    ## 268     PHI:8167                                  347
    ## 269     PHI:8170                                  347
    ## 270     PHI:8172                                  347
    ## 271     PHI:8173                                  347
    ## 272     PHI:8174                                  347
    ## 273     PHI:8175                                  347
    ## 274     PHI:8422                                  317
    ## 275     PHI:8422                                  317
    ## 276     PHI:8422                                  317
    ## 277     PHI:8422                                  317
    ##                    Pathogenspecies PathogenstrainID
    ## 1                   Myzus persicae               NA
    ## 2                   Myzus persicae               NA
    ## 3                   Myzus persicae               NA
    ## 4             Meloidogyne javanica               NA
    ## 5                   Myzus persicae               NA
    ## 6                   Myzus persicae               NA
    ## 7                   Myzus persicae               NA
    ## 8             Meloidogyne javanica               NA
    ## 9             Setosphaeria turcica           671987
    ## 10            Setosphaeria turcica           671987
    ## 11         Phytophthora parasitica               NA
    ## 12              Magnaporthe oryzae           242507
    ## 13              Magnaporthe oryzae           242507
    ## 14      Melampsora larici-populina           747676
    ## 15      Melampsora larici-populina           747676
    ## 16      Melampsora larici-populina           747676
    ## 17      Melampsora larici-populina           747676
    ## 18      Melampsora larici-populina           747676
    ## 19      Melampsora larici-populina           747676
    ## 20      Melampsora larici-populina           747676
    ## 21      Melampsora larici-populina           747676
    ## 22      Melampsora larici-populina           747676
    ## 23      Melampsora larici-populina           747676
    ## 24      Melampsora larici-populina           747676
    ## 25      Melampsora larici-populina           747676
    ## 26      Melampsora larici-populina           747676
    ## 27      Melampsora larici-populina           747676
    ## 28      Melampsora larici-populina           747676
    ## 29      Melampsora larici-populina           747676
    ## 30      Melampsora larici-populina           747676
    ## 31      Melampsora larici-populina           747676
    ## 32      Melampsora larici-populina           747676
    ## 33      Melampsora larici-populina           747676
    ## 34      Melampsora larici-populina           747676
    ## 35      Melampsora larici-populina           747676
    ## 36      Melampsora larici-populina           747676
    ## 37      Melampsora larici-populina           747676
    ## 38            Phytophthora capsici               NA
    ## 39            Phytophthora capsici               NA
    ## 40            Phytophthora capsici               NA
    ## 41            Phytophthora capsici               NA
    ## 42            Phytophthora capsici               NA
    ## 43            Phytophthora capsici               NA
    ## 44                 Bremia lactucae               NA
    ## 45                 Bremia lactucae               NA
    ## 46            Pseudomonas syringae           629260
    ## 47          Phytophthora infestans               NA
    ## 48  Hyaloperonospora arabidopsidis           559515
    ## 49            Phytophthora capsici               NA
    ## 50            Phytophthora capsici               NA
    ## 51            Phytophthora capsici               NA
    ## 52            Phytophthora capsici               NA
    ## 53            Phytophthora capsici               NA
    ## 54            Phytophthora capsici               NA
    ## 55      Pseudocercospora fijiensis           383855
    ## 56      Pseudocercospora fijiensis           383855
    ## 57      Pseudocercospora fijiensis           383855
    ## 58           Phytophthora cactorum               NA
    ## 59           Phytophthora cactorum               NA
    ## 60          Ralstonia solanacearum           267608
    ## 61          Ralstonia solanacearum           267608
    ## 62          Ralstonia solanacearum           267608
    ## 63          Ralstonia solanacearum           267608
    ## 64          Ralstonia solanacearum           267608
    ## 65          Ralstonia solanacearum           267608
    ## 66          Ralstonia solanacearum           267608
    ## 67          Ralstonia solanacearum           267608
    ## 68          Ralstonia solanacearum           267608
    ## 69          Ralstonia solanacearum           267608
    ## 70          Ralstonia solanacearum           267608
    ## 71          Ralstonia solanacearum           267608
    ## 72          Ralstonia solanacearum           267608
    ## 73          Ralstonia solanacearum           267608
    ## 74          Ralstonia solanacearum           267608
    ## 75          Ralstonia solanacearum           267608
    ## 76          Ralstonia solanacearum           267608
    ## 77          Ralstonia solanacearum           267608
    ## 78          Ralstonia solanacearum           267608
    ## 79          Ralstonia solanacearum           267608
    ## 80          Ralstonia solanacearum           267608
    ## 81          Ralstonia solanacearum           267608
    ## 82          Ralstonia solanacearum           267608
    ## 83          Ralstonia solanacearum           267608
    ## 84          Ralstonia solanacearum           267608
    ## 85          Ralstonia solanacearum           267608
    ## 86          Ralstonia solanacearum           267608
    ## 87          Ralstonia solanacearum           267608
    ## 88          Ralstonia solanacearum           267608
    ## 89          Ralstonia solanacearum           267608
    ## 90          Ralstonia solanacearum           267608
    ## 91          Ralstonia solanacearum           267608
    ## 92          Ralstonia solanacearum           267608
    ## 93          Ralstonia solanacearum           267608
    ## 94          Ralstonia solanacearum           267608
    ## 95          Ralstonia solanacearum           267608
    ## 96          Ralstonia solanacearum           267608
    ## 97          Ralstonia solanacearum           267608
    ## 98          Ralstonia solanacearum           267608
    ## 99          Ralstonia solanacearum           267608
    ## 100         Ralstonia solanacearum           267608
    ## 101         Ralstonia solanacearum           267608
    ## 102         Ralstonia solanacearum           267608
    ## 103         Ralstonia solanacearum           267608
    ## 104         Ralstonia solanacearum           267608
    ## 105         Ralstonia solanacearum           267608
    ## 106         Ralstonia solanacearum           267608
    ## 107         Ralstonia solanacearum           267608
    ## 108         Ralstonia solanacearum           267608
    ## 109         Ralstonia solanacearum           267608
    ## 110         Ralstonia solanacearum           267608
    ## 111         Ralstonia solanacearum           267608
    ## 112         Ralstonia solanacearum           267608
    ## 113         Ralstonia solanacearum           267608
    ## 114         Ralstonia solanacearum           267608
    ## 115         Ralstonia solanacearum           267608
    ## 116         Ralstonia solanacearum           267608
    ## 117         Ralstonia solanacearum           267608
    ## 118         Ralstonia solanacearum           267608
    ## 119         Ralstonia solanacearum           267608
    ## 120         Ralstonia solanacearum           267608
    ## 121         Ralstonia solanacearum           267608
    ## 122         Ralstonia solanacearum           267608
    ## 123         Ralstonia solanacearum           267608
    ## 124         Ralstonia solanacearum           267608
    ## 125         Ralstonia solanacearum           267608
    ## 126         Ralstonia solanacearum           267608
    ## 127         Ralstonia solanacearum           267608
    ## 128         Ralstonia solanacearum           267608
    ## 129         Ralstonia solanacearum           267608
    ## 130         Ralstonia solanacearum           267608
    ## 131         Ralstonia solanacearum           267608
    ## 132         Ralstonia solanacearum           267608
    ## 133         Ralstonia solanacearum           267608
    ## 134         Ralstonia solanacearum           267608
    ## 135         Ralstonia solanacearum           267608
    ## 136         Ralstonia solanacearum           267608
    ## 137         Ralstonia solanacearum           267608
    ## 138         Ralstonia solanacearum           267608
    ## 139         Ralstonia solanacearum           267608
    ## 140         Ralstonia solanacearum           267608
    ## 141         Ralstonia solanacearum           267608
    ## 142         Ralstonia solanacearum           267608
    ## 143         Ralstonia solanacearum           267608
    ## 144         Ralstonia solanacearum           267608
    ## 145         Ralstonia solanacearum           267608
    ## 146         Ralstonia solanacearum           267608
    ## 147         Ralstonia solanacearum           267608
    ## 148         Ralstonia solanacearum           267608
    ## 149         Ralstonia solanacearum           267608
    ## 150         Ralstonia solanacearum           267608
    ## 151         Ralstonia solanacearum           267608
    ## 152         Ralstonia solanacearum           267608
    ## 153         Ralstonia solanacearum           267608
    ## 154         Ralstonia solanacearum           267608
    ## 155         Ralstonia solanacearum           267608
    ## 156         Ralstonia solanacearum           267608
    ## 157         Ralstonia solanacearum           267608
    ## 158         Ralstonia solanacearum           267608
    ## 159         Ralstonia solanacearum           267608
    ## 160         Ralstonia solanacearum           267608
    ## 161         Ralstonia solanacearum           267608
    ## 162         Ralstonia solanacearum           267608
    ## 163         Ralstonia solanacearum           267608
    ## 164         Ralstonia solanacearum           267608
    ## 165         Ralstonia solanacearum           267608
    ## 166         Ralstonia solanacearum           267608
    ## 167         Ralstonia solanacearum           267608
    ## 168         Ralstonia solanacearum           267608
    ## 169         Ralstonia solanacearum           267608
    ## 170         Ralstonia solanacearum           267608
    ## 171         Ralstonia solanacearum           267608
    ## 172         Ralstonia solanacearum           267608
    ## 173         Ralstonia solanacearum           267608
    ## 174         Ralstonia solanacearum           267608
    ## 175         Ralstonia solanacearum           267608
    ## 176         Ralstonia solanacearum           267608
    ## 177         Ralstonia solanacearum           267608
    ## 178         Ralstonia solanacearum           267608
    ## 179         Ralstonia solanacearum           267608
    ## 180         Ralstonia solanacearum           267608
    ## 181         Ralstonia solanacearum           267608
    ## 182         Ralstonia solanacearum           267608
    ## 183         Ralstonia solanacearum           267608
    ## 184         Ralstonia solanacearum           267608
    ## 185         Ralstonia solanacearum           267608
    ## 186         Ralstonia solanacearum           267608
    ## 187         Ralstonia solanacearum           267608
    ## 188         Ralstonia solanacearum           267608
    ## 189         Ralstonia solanacearum           267608
    ## 190         Ralstonia solanacearum           267608
    ## 191         Ralstonia solanacearum           267608
    ## 192         Ralstonia solanacearum           267608
    ## 193         Ralstonia solanacearum           267608
    ## 194         Ralstonia solanacearum           267608
    ## 195         Ralstonia solanacearum           267608
    ## 196         Ralstonia solanacearum           267608
    ## 197         Ralstonia solanacearum           267608
    ## 198         Ralstonia solanacearum           267608
    ## 199         Ralstonia solanacearum           267608
    ## 200         Ralstonia solanacearum           267608
    ## 201         Ralstonia solanacearum           267608
    ## 202         Ralstonia solanacearum           267608
    ## 203         Ralstonia solanacearum           267608
    ## 204         Ralstonia solanacearum           267608
    ## 205         Ralstonia solanacearum           267608
    ## 206         Ralstonia solanacearum           267608
    ## 207         Ralstonia solanacearum           267608
    ## 208         Ralstonia solanacearum           267608
    ## 209         Ralstonia solanacearum           267608
    ## 210         Ralstonia solanacearum           267608
    ## 211         Ralstonia solanacearum           267608
    ## 212         Ralstonia solanacearum           267608
    ## 213         Ralstonia solanacearum           267608
    ## 214         Ralstonia solanacearum           267608
    ## 215         Ralstonia solanacearum           267608
    ## 216         Ralstonia solanacearum           267608
    ## 217         Ralstonia solanacearum           267608
    ## 218         Ralstonia solanacearum           267608
    ## 219         Ralstonia solanacearum           267608
    ## 220         Ralstonia solanacearum           267608
    ## 221         Ralstonia solanacearum           267608
    ## 222         Ralstonia solanacearum           267608
    ## 223         Ralstonia solanacearum           267608
    ## 224         Ralstonia solanacearum           267608
    ## 225         Ralstonia solanacearum           267608
    ## 226         Ralstonia solanacearum           267608
    ## 227         Ralstonia solanacearum           267608
    ## 228         Ralstonia solanacearum           267608
    ## 229         Ralstonia solanacearum           267608
    ## 230         Ralstonia solanacearum           267608
    ## 231         Ralstonia solanacearum           267608
    ## 232             Fusarium oxysporum           396571
    ## 233             Fusarium oxysporum           396571
    ## 234             Fusarium oxysporum           396571
    ## 235             Fusarium oxysporum           396571
    ## 236             Fusarium oxysporum           396571
    ## 237             Fusarium oxysporum           396571
    ## 238             Fusarium oxysporum           396571
    ## 239             Fusarium oxysporum           396571
    ## 240             Fusarium oxysporum           396571
    ## 241             Fusarium oxysporum           396571
    ## 242             Magnaporthe oryzae               NA
    ## 243             Magnaporthe oryzae               NA
    ## 244             Magnaporthe oryzae               NA
    ## 245             Magnaporthe oryzae               NA
    ## 246             Magnaporthe oryzae               NA
    ## 247           Puccinia striiformis           168172
    ## 248           Puccinia striiformis           168172
    ## 249           Puccinia striiformis           168172
    ## 250           Puccinia striiformis           168172
    ## 251           Puccinia striiformis           168172
    ## 252           Puccinia striiformis           168172
    ## 253           Puccinia striiformis           168172
    ## 254           Puccinia striiformis           168172
    ## 255           Puccinia striiformis           168172
    ## 256           Puccinia striiformis           168172
    ## 257           Puccinia striiformis           168172
    ## 258           Puccinia striiformis           168172
    ## 259           Puccinia striiformis           168172
    ## 260           Puccinia striiformis           168172
    ## 261           Puccinia striiformis           168172
    ## 262           Phytophthora capsici           763924
    ## 263         Xanthomonas axonopodis           487838
    ## 264              Puccinia graminis            56615
    ## 265              Puccinia graminis            56615
    ## 266              Puccinia graminis            56615
    ## 267              Puccinia graminis            56615
    ## 268             Xanthomonas oryzae          1009853
    ## 269             Xanthomonas oryzae          1009853
    ## 270             Xanthomonas oryzae          1009853
    ## 271             Xanthomonas oryzae          1009853
    ## 272             Xanthomonas oryzae          1009853
    ## 273             Xanthomonas oryzae          1009853
    ## 274           Pseudomonas syringae               NA
    ## 275           Pseudomonas syringae           264730
    ## 276           Pseudomonas syringae           223283
    ## 277           Pseudomonas syringae           264730
    ##                   Pathogenstrain Hostdescription HostID
    ## 1              RRes (genotype O)        Eudicots   4100
    ## 2              RRes (genotype O)        Eudicots   4100
    ## 3              RRes (genotype O)        Eudicots   4100
    ## 4                            VW4        Eudicots   4081
    ## 5              RRes (genotype O)        Eudicots   4100
    ## 6              RRes (genotype O)        Eudicots   4100
    ## 7              RRes (genotype O)        Eudicots   4100
    ## 8                            VW4        Eudicots   4081
    ## 9                          Et28A        Monocots   4577
    ## 10                         Et28A        Monocots   4577
    ## 11                           310        Eudicots   3702
    ## 12                         70-15        Monocots   4530
    ## 13                         70-15        Eudicots   4100
    ## 14                        98AG31        Eudicots   4100
    ## 15                        98AG31        Eudicots   4100
    ## 16                        98AG31        Eudicots   4100
    ## 17                        98AG31        Eudicots   4100
    ## 18                        98AG31        Eudicots   4100
    ## 19                        98AG31        Eudicots   4100
    ## 20                        98AG31        Eudicots   4100
    ## 21                        98AG31        Eudicots   4100
    ## 22                        98AG31        Eudicots   4100
    ## 23                        98AG31        Eudicots   4100
    ## 24                        98AG31        Eudicots   4100
    ## 25                        98AG31        Eudicots   4100
    ## 26                        98AG31        Eudicots   4100
    ## 27                        98AG31        Eudicots   4100
    ## 28                        98AG31        Eudicots   4100
    ## 29                        98AG31        Eudicots   4100
    ## 30                        98AG31        Eudicots   4100
    ## 31                        98AG31        Eudicots   4100
    ## 32                        98AG31        Eudicots   4100
    ## 33                        98AG31        Eudicots   4100
    ## 34                        98AG31        Eudicots   4100
    ## 35                        98AG31        Eudicots   4100
    ## 36                        98AG31        Eudicots   4100
    ## 37                        98AG31        Eudicots   4100
    ## 38                                      Eudicots   4100
    ## 39                                      Eudicots   4100
    ## 40                                      Eudicots   4100
    ## 41                                      Eudicots   4097
    ## 42                                      Eudicots   4097
    ## 43                                      Eudicots   4097
    ## 44                                      Eudicots   4236
    ## 45                                      Eudicots   4236
    ## 46        pv. lachrymans M301315        Eudicots   4072
    ## 47                         88069        Eudicots   3702
    ## 48                         Emoy2        Eudicots   3702
    ## 49                         LT263        Eudicots   4100
    ## 50                         LT263        Eudicots   3702
    ## 51                         LT263        Eudicots   4100
    ## 52                         LT263        Eudicots   3702
    ## 53                         LT263        Eudicots   4100
    ## 54                         LT263        Eudicots   3702
    ## 55                       CIRAD86        Eudicots   4081
    ## 56                       CIRAD86        Eudicots   4100
    ## 57                       CIRAD86        Eudicots   4081
    ## 58                         10300        Eudicots   4081
    ## 59                         10300        Eudicots   4097
    ## 60                       GMI1000        Eudicots   3880
    ## 61                                      Eudicots   4111
    ## 62                                      Eudicots   4111
    ## 63                                      Eudicots   4111
    ## 64                                      Eudicots   4111
    ## 65                                      Eudicots   4111
    ## 66                                      Eudicots   4111
    ## 67                                      Eudicots   4111
    ## 68                                      Eudicots   4111
    ## 69                                      Eudicots   4111
    ## 70                                      Eudicots   4111
    ## 71                                      Eudicots   4111
    ## 72                                      Eudicots   4111
    ## 73                                      Eudicots   4111
    ## 74                                      Eudicots   4111
    ## 75                                      Eudicots   4111
    ## 76                                      Eudicots   4111
    ## 77                                      Eudicots   4111
    ## 78                                      Eudicots   4111
    ## 79                                      Eudicots   4111
    ## 80                                      Eudicots   4111
    ## 81                                      Eudicots   4111
    ## 82                                      Eudicots   4111
    ## 83                                      Eudicots   4111
    ## 84                                      Eudicots   4111
    ## 85                                      Eudicots   4111
    ## 86                                      Eudicots   4111
    ## 87                                      Eudicots   4111
    ## 88                                      Eudicots   4111
    ## 89                                      Eudicots   4111
    ## 90                                      Eudicots   4111
    ## 91                                      Eudicots   4111
    ## 92                                      Eudicots   4111
    ## 93                                      Eudicots   4111
    ## 94                                      Eudicots   4111
    ## 95                                      Eudicots   4111
    ## 96                                      Eudicots   4111
    ## 97                                      Eudicots   4111
    ## 98                                      Eudicots   4111
    ## 99                                      Eudicots   4111
    ## 100                                     Eudicots   4111
    ## 101                                     Eudicots   4111
    ## 102                                     Eudicots   4111
    ## 103                                     Eudicots   4111
    ## 104                                     Eudicots   4111
    ## 105                                     Eudicots   4111
    ## 106                                     Eudicots   4111
    ## 107                                     Eudicots   4111
    ## 108                                     Eudicots   4111
    ## 109                                     Eudicots   4111
    ## 110                                     Eudicots   4111
    ## 111                                     Eudicots   4111
    ## 112                                     Eudicots   4111
    ## 113                                     Eudicots   4111
    ## 114                                     Eudicots   4111
    ## 115                                     Eudicots   4111
    ## 116                                     Eudicots   4111
    ## 117                                     Eudicots   4111
    ## 118                                     Eudicots   4111
    ## 119                                     Eudicots   4111
    ## 120                                     Eudicots   4111
    ## 121                                     Eudicots   4111
    ## 122                                     Eudicots   4111
    ## 123                                     Eudicots   4111
    ## 124                                     Eudicots   4111
    ## 125                                     Eudicots   4111
    ## 126                                     Eudicots   4111
    ## 127                                     Eudicots   4111
    ## 128                                     Eudicots   4111
    ## 129                                     Eudicots   4111
    ## 130                                     Eudicots   4111
    ## 131                                     Eudicots   4111
    ## 132                                     Eudicots   4111
    ## 133                                     Eudicots   4111
    ## 134                                     Eudicots   4111
    ## 135                                     Eudicots   4111
    ## 136                                     Eudicots   4111
    ## 137                                     Eudicots   4111
    ## 138                                     Eudicots   4111
    ## 139                                     Eudicots   4111
    ## 140                                     Eudicots   4111
    ## 141                                     Eudicots   4111
    ## 142                                     Eudicots   4111
    ## 143                                     Eudicots   4111
    ## 144                                     Eudicots   4111
    ## 145                                     Eudicots   4111
    ## 146                                     Eudicots   4111
    ## 147                                     Eudicots   4111
    ## 148                                     Eudicots   4111
    ## 149                                     Eudicots   4111
    ## 150                                     Eudicots   4111
    ## 151                                     Eudicots   4111
    ## 152                                     Eudicots   4111
    ## 153                                     Eudicots   4111
    ## 154                                     Eudicots   4111
    ## 155                                     Eudicots   4111
    ## 156                                     Eudicots   4072
    ## 157                                     Eudicots   4072
    ## 158                                     Eudicots   4072
    ## 159                                     Eudicots   4072
    ## 160                                     Eudicots   4072
    ## 161                                     Eudicots   4072
    ## 162                                     Eudicots   4072
    ## 163                                     Eudicots   4072
    ## 164                                     Eudicots   4072
    ## 165                                     Eudicots   4072
    ## 166                                     Eudicots   4072
    ## 167                                     Eudicots   4072
    ## 168                                     Eudicots   4072
    ## 169                                     Eudicots   4072
    ## 170                                     Eudicots   4072
    ## 171                                     Eudicots   4072
    ## 172                                     Eudicots   4072
    ## 173                                     Eudicots   4072
    ## 174                                     Eudicots   4072
    ## 175                                     Eudicots   4072
    ## 176                                     Eudicots   4072
    ## 177                                     Eudicots   4072
    ## 178                                     Eudicots   4072
    ## 179                                     Eudicots   4072
    ## 180                                     Eudicots   4072
    ## 181                                     Eudicots   4072
    ## 182                                     Eudicots   4072
    ## 183                                     Eudicots   4072
    ## 184                                     Eudicots   4072
    ## 185                                     Eudicots   4072
    ## 186                                     Eudicots   4072
    ## 187                                     Eudicots   4072
    ## 188                                     Eudicots   4072
    ## 189                                     Eudicots   4072
    ## 190                                     Eudicots   4072
    ## 191                                     Eudicots   4072
    ## 192                                     Eudicots   4072
    ## 193                                     Eudicots   4072
    ## 194                                     Eudicots   4072
    ## 195                                     Eudicots   4072
    ## 196                                     Eudicots   4072
    ## 197                                     Eudicots   4072
    ## 198                                     Eudicots   4072
    ## 199                                     Eudicots   4072
    ## 200                                     Eudicots   4072
    ## 201                                     Eudicots   4072
    ## 202                                     Eudicots   4072
    ## 203                                     Eudicots   4072
    ## 204                                     Eudicots   4072
    ## 205                                     Eudicots   4072
    ## 206                                     Eudicots   4072
    ## 207                                     Eudicots   4072
    ## 208                                     Eudicots   4072
    ## 209                                     Eudicots   4072
    ## 210                                     Eudicots   4072
    ## 211                                     Eudicots   4072
    ## 212                                     Eudicots   4072
    ## 213                                     Eudicots   4081
    ## 214                                     Eudicots   4081
    ## 215                                     Eudicots   4081
    ## 216                                     Eudicots   4081
    ## 217                                     Eudicots   4081
    ## 218                                     Eudicots   4081
    ## 219                                     Eudicots   4081
    ## 220                                     Eudicots   4081
    ## 221                                     Eudicots   4081
    ## 222                                     Eudicots   4081
    ## 223                                     Eudicots   4081
    ## 224                                     Eudicots   4081
    ## 225                                     Eudicots   4081
    ## 226                                     Eudicots   4081
    ## 227                                     Eudicots   4081
    ## 228                                     Eudicots   4081
    ## 229                                     Eudicots   4081
    ## 230                                     Eudicots   4081
    ## 231                                     Eudicots   4081
    ## 232                 f. sp. cepae        Monocots   4679
    ## 233                 f. sp. cepae        Monocots   4679
    ## 234                 f. sp. cepae        Monocots   4679
    ## 235                 f. sp. cepae        Monocots   4679
    ## 236                 f. sp. cepae        Monocots   4679
    ## 237                 f. sp. cepae        Monocots   4679
    ## 238                 f. sp. cepae        Monocots   4679
    ## 239                 f. sp. cepae        Monocots   4679
    ## 240                 f. sp. cepae        Monocots   4679
    ## 241                 f. sp. cepae        Monocots   4679
    ## 242                   08/05/4091        Monocots   4530
    ## 243                   08/05/4091        Monocots   4530
    ## 244                   08/05/4091        Monocots   4530
    ## 245                        INA72        Eudicots   4100
    ## 246                        INA72        Eudicots   4100
    ## 247               f. sp. tritici        Eudicots   4100
    ## 248               f. sp. tritici        Eudicots   4100
    ## 249               f. sp. tritici        Eudicots   4100
    ## 250               f. sp. tritici        Eudicots   4100
    ## 251               f. sp. tritici        Eudicots   4100
    ## 252               f. sp. tritici        Eudicots   4100
    ## 253               f. sp. tritici        Eudicots   4100
    ## 254               f. sp. tritici        Eudicots   4100
    ## 255               f. sp. tritici        Eudicots   4100
    ## 256               f. sp. tritici        Eudicots   4100
    ## 257               f. sp. tritici        Eudicots   4100
    ## 258               f. sp. tritici        Eudicots   4100
    ## 259               f. sp. tritici        Eudicots   4100
    ## 260               f. sp. tritici        Eudicots   4100
    ## 261               f. sp. tritici        Eudicots   4100
    ## 262                       LT1534        Eudicots   4100
    ## 263 pv. Punicae str. ITCCBD 0003        Eudicots  22663
    ## 264               f. sp. tritici        Eudicots   4100
    ## 265               f. sp. tritici        Eudicots   4097
    ## 266               f. sp. tritici        Monocots   4565
    ## 267               f. sp. tritici        Monocots   4565
    ## 268                       X11-5A        Monocots   4530
    ## 269                       X11-5A        Monocots   4530
    ## 270                       X11-5A        Monocots   4530
    ## 271                       X11-5A        Monocots   4530
    ## 272                       X11-5A        Monocots   4530
    ## 273                       X11-5A        Monocots   4530
    ## 274     pv. syringae strain 7B40        Eudicots   3885
    ## 275       pv. phaseolicola 1448A        Eudicots   3885
    ## 276        pv. Tomato str DC3000        Eudicots   4100
    ## 277       pv. phaseolicola 1448A        Eudicots   4081
    ##                                      Hostspecies
    ## 1   Nicotiana benthamiana (no common name found)
    ## 2   Nicotiana benthamiana (no common name found)
    ## 3   Nicotiana benthamiana (no common name found)
    ## 4         Solanum lycopersicum (related: tomato)
    ## 5   Nicotiana benthamiana (no common name found)
    ## 6   Nicotiana benthamiana (no common name found)
    ## 7   Nicotiana benthamiana (no common name found)
    ## 8         Solanum lycopersicum (related: tomato)
    ## 9                      Zea mays (related: maize)
    ## 10                     Zea mays (related: maize)
    ## 11   Arabidopsis thaliana (related: thale cress)
    ## 12                  Oryza sativa (related: rice)
    ## 13  Nicotiana benthamiana (no common name found)
    ## 14  Nicotiana benthamiana (no common name found)
    ## 15  Nicotiana benthamiana (no common name found)
    ## 16  Nicotiana benthamiana (no common name found)
    ## 17  Nicotiana benthamiana (no common name found)
    ## 18  Nicotiana benthamiana (no common name found)
    ## 19  Nicotiana benthamiana (no common name found)
    ## 20  Nicotiana benthamiana (no common name found)
    ## 21  Nicotiana benthamiana (no common name found)
    ## 22  Nicotiana benthamiana (no common name found)
    ## 23  Nicotiana benthamiana (no common name found)
    ## 24  Nicotiana benthamiana (no common name found)
    ## 25  Nicotiana benthamiana (no common name found)
    ## 26  Nicotiana benthamiana (no common name found)
    ## 27  Nicotiana benthamiana (no common name found)
    ## 28  Nicotiana benthamiana (no common name found)
    ## 29  Nicotiana benthamiana (no common name found)
    ## 30  Nicotiana benthamiana (no common name found)
    ## 31  Nicotiana benthamiana (no common name found)
    ## 32  Nicotiana benthamiana (no common name found)
    ## 33  Nicotiana benthamiana (no common name found)
    ## 34  Nicotiana benthamiana (no common name found)
    ## 35  Nicotiana benthamiana (no common name found)
    ## 36  Nicotiana benthamiana (no common name found)
    ## 37  Nicotiana benthamiana (no common name found)
    ## 38  Nicotiana benthamiana (no common name found)
    ## 39  Nicotiana benthamiana (no common name found)
    ## 40  Nicotiana benthamiana (no common name found)
    ## 41   Nicotiana tabacum (related: common tobacco)
    ## 42   Nicotiana tabacum (related: common tobacco)
    ## 43   Nicotiana tabacum (related: common tobacco)
    ## 44      Lactuca sativa (related: garden lettuce)
    ## 45      Lactuca sativa (related: garden lettuce)
    ## 46        Capsicum annuum (no common name found)
    ## 47   Arabidopsis thaliana (related: thale cress)
    ## 48   Arabidopsis thaliana (related: thale cress)
    ## 49  Nicotiana benthamiana (no common name found)
    ## 50   Arabidopsis thaliana (related: thale cress)
    ## 51  Nicotiana benthamiana (no common name found)
    ## 52   Arabidopsis thaliana (related: thale cress)
    ## 53  Nicotiana benthamiana (no common name found)
    ## 54   Arabidopsis thaliana (related: thale cress)
    ## 55        Solanum lycopersicum (related: tomato)
    ## 56  Nicotiana benthamiana (no common name found)
    ## 57        Solanum lycopersicum (related: tomato)
    ## 58        Solanum lycopersicum (related: tomato)
    ## 59   Nicotiana tabacum (related: common tobacco)
    ## 60   Medicago truncatula (related: barrel medic)
    ## 61         Solanum melongena (related: eggplant)
    ## 62         Solanum melongena (related: eggplant)
    ## 63         Solanum melongena (related: eggplant)
    ## 64         Solanum melongena (related: eggplant)
    ## 65         Solanum melongena (related: eggplant)
    ## 66         Solanum melongena (related: eggplant)
    ## 67         Solanum melongena (related: eggplant)
    ## 68         Solanum melongena (related: eggplant)
    ## 69         Solanum melongena (related: eggplant)
    ## 70         Solanum melongena (related: eggplant)
    ## 71         Solanum melongena (related: eggplant)
    ## 72         Solanum melongena (related: eggplant)
    ## 73         Solanum melongena (related: eggplant)
    ## 74         Solanum melongena (related: eggplant)
    ## 75         Solanum melongena (related: eggplant)
    ## 76         Solanum melongena (related: eggplant)
    ## 77         Solanum melongena (related: eggplant)
    ## 78         Solanum melongena (related: eggplant)
    ## 79         Solanum melongena (related: eggplant)
    ## 80         Solanum melongena (related: eggplant)
    ## 81         Solanum melongena (related: eggplant)
    ## 82         Solanum melongena (related: eggplant)
    ## 83         Solanum melongena (related: eggplant)
    ## 84         Solanum melongena (related: eggplant)
    ## 85         Solanum melongena (related: eggplant)
    ## 86         Solanum melongena (related: eggplant)
    ## 87         Solanum melongena (related: eggplant)
    ## 88         Solanum melongena (related: eggplant)
    ## 89         Solanum melongena (related: eggplant)
    ## 90         Solanum melongena (related: eggplant)
    ## 91         Solanum melongena (related: eggplant)
    ## 92         Solanum melongena (related: eggplant)
    ## 93         Solanum melongena (related: eggplant)
    ## 94         Solanum melongena (related: eggplant)
    ## 95         Solanum melongena (related: eggplant)
    ## 96         Solanum melongena (related: eggplant)
    ## 97         Solanum melongena (related: eggplant)
    ## 98         Solanum melongena (related: eggplant)
    ## 99         Solanum melongena (related: eggplant)
    ## 100        Solanum melongena (related: eggplant)
    ## 101        Solanum melongena (related: eggplant)
    ## 102        Solanum melongena (related: eggplant)
    ## 103        Solanum melongena (related: eggplant)
    ## 104        Solanum melongena (related: eggplant)
    ## 105        Solanum melongena (related: eggplant)
    ## 106        Solanum melongena (related: eggplant)
    ## 107        Solanum melongena (related: eggplant)
    ## 108        Solanum melongena (related: eggplant)
    ## 109        Solanum melongena (related: eggplant)
    ## 110        Solanum melongena (related: eggplant)
    ## 111        Solanum melongena (related: eggplant)
    ## 112        Solanum melongena (related: eggplant)
    ## 113        Solanum melongena (related: eggplant)
    ## 114        Solanum melongena (related: eggplant)
    ## 115        Solanum melongena (related: eggplant)
    ## 116        Solanum melongena (related: eggplant)
    ## 117        Solanum melongena (related: eggplant)
    ## 118        Solanum melongena (related: eggplant)
    ## 119        Solanum melongena (related: eggplant)
    ## 120        Solanum melongena (related: eggplant)
    ## 121        Solanum melongena (related: eggplant)
    ## 122        Solanum melongena (related: eggplant)
    ## 123        Solanum melongena (related: eggplant)
    ## 124        Solanum melongena (related: eggplant)
    ## 125        Solanum melongena (related: eggplant)
    ## 126        Solanum melongena (related: eggplant)
    ## 127        Solanum melongena (related: eggplant)
    ## 128        Solanum melongena (related: eggplant)
    ## 129        Solanum melongena (related: eggplant)
    ## 130        Solanum melongena (related: eggplant)
    ## 131        Solanum melongena (related: eggplant)
    ## 132        Solanum melongena (related: eggplant)
    ## 133        Solanum melongena (related: eggplant)
    ## 134        Solanum melongena (related: eggplant)
    ## 135        Solanum melongena (related: eggplant)
    ## 136        Solanum melongena (related: eggplant)
    ## 137        Solanum melongena (related: eggplant)
    ## 138        Solanum melongena (related: eggplant)
    ## 139        Solanum melongena (related: eggplant)
    ## 140        Solanum melongena (related: eggplant)
    ## 141        Solanum melongena (related: eggplant)
    ## 142        Solanum melongena (related: eggplant)
    ## 143        Solanum melongena (related: eggplant)
    ## 144        Solanum melongena (related: eggplant)
    ## 145        Solanum melongena (related: eggplant)
    ## 146        Solanum melongena (related: eggplant)
    ## 147        Solanum melongena (related: eggplant)
    ## 148        Solanum melongena (related: eggplant)
    ## 149        Solanum melongena (related: eggplant)
    ## 150        Solanum melongena (related: eggplant)
    ## 151        Solanum melongena (related: eggplant)
    ## 152        Solanum melongena (related: eggplant)
    ## 153        Solanum melongena (related: eggplant)
    ## 154        Solanum melongena (related: eggplant)
    ## 155        Solanum melongena (related: eggplant)
    ## 156       Capsicum annuum (no common name found)
    ## 157       Capsicum annuum (no common name found)
    ## 158       Capsicum annuum (no common name found)
    ## 159       Capsicum annuum (no common name found)
    ## 160       Capsicum annuum (no common name found)
    ## 161       Capsicum annuum (no common name found)
    ## 162       Capsicum annuum (no common name found)
    ## 163       Capsicum annuum (no common name found)
    ## 164       Capsicum annuum (no common name found)
    ## 165       Capsicum annuum (no common name found)
    ## 166       Capsicum annuum (no common name found)
    ## 167       Capsicum annuum (no common name found)
    ## 168       Capsicum annuum (no common name found)
    ## 169       Capsicum annuum (no common name found)
    ## 170       Capsicum annuum (no common name found)
    ## 171       Capsicum annuum (no common name found)
    ## 172       Capsicum annuum (no common name found)
    ## 173       Capsicum annuum (no common name found)
    ## 174       Capsicum annuum (no common name found)
    ## 175       Capsicum annuum (no common name found)
    ## 176       Capsicum annuum (no common name found)
    ## 177       Capsicum annuum (no common name found)
    ## 178       Capsicum annuum (no common name found)
    ## 179       Capsicum annuum (no common name found)
    ## 180       Capsicum annuum (no common name found)
    ## 181       Capsicum annuum (no common name found)
    ## 182       Capsicum annuum (no common name found)
    ## 183       Capsicum annuum (no common name found)
    ## 184       Capsicum annuum (no common name found)
    ## 185       Capsicum annuum (no common name found)
    ## 186       Capsicum annuum (no common name found)
    ## 187       Capsicum annuum (no common name found)
    ## 188       Capsicum annuum (no common name found)
    ## 189       Capsicum annuum (no common name found)
    ## 190       Capsicum annuum (no common name found)
    ## 191       Capsicum annuum (no common name found)
    ## 192       Capsicum annuum (no common name found)
    ## 193       Capsicum annuum (no common name found)
    ## 194       Capsicum annuum (no common name found)
    ## 195       Capsicum annuum (no common name found)
    ## 196       Capsicum annuum (no common name found)
    ## 197       Capsicum annuum (no common name found)
    ## 198       Capsicum annuum (no common name found)
    ## 199       Capsicum annuum (no common name found)
    ## 200       Capsicum annuum (no common name found)
    ## 201       Capsicum annuum (no common name found)
    ## 202       Capsicum annuum (no common name found)
    ## 203       Capsicum annuum (no common name found)
    ## 204       Capsicum annuum (no common name found)
    ## 205       Capsicum annuum (no common name found)
    ## 206       Capsicum annuum (no common name found)
    ## 207       Capsicum annuum (no common name found)
    ## 208       Capsicum annuum (no common name found)
    ## 209       Capsicum annuum (no common name found)
    ## 210       Capsicum annuum (no common name found)
    ## 211       Capsicum annuum (no common name found)
    ## 212       Capsicum annuum (no common name found)
    ## 213       Solanum lycopersicum (related: tomato)
    ## 214       Solanum lycopersicum (related: tomato)
    ## 215       Solanum lycopersicum (related: tomato)
    ## 216       Solanum lycopersicum (related: tomato)
    ## 217       Solanum lycopersicum (related: tomato)
    ## 218       Solanum lycopersicum (related: tomato)
    ## 219       Solanum lycopersicum (related: tomato)
    ## 220       Solanum lycopersicum (related: tomato)
    ## 221       Solanum lycopersicum (related: tomato)
    ## 222       Solanum lycopersicum (related: tomato)
    ## 223       Solanum lycopersicum (related: tomato)
    ## 224       Solanum lycopersicum (related: tomato)
    ## 225       Solanum lycopersicum (related: tomato)
    ## 226       Solanum lycopersicum (related: tomato)
    ## 227       Solanum lycopersicum (related: tomato)
    ## 228       Solanum lycopersicum (related: tomato)
    ## 229       Solanum lycopersicum (related: tomato)
    ## 230       Solanum lycopersicum (related: tomato)
    ## 231       Solanum lycopersicum (related: tomato)
    ## 232                 Allium cepa (related: onion)
    ## 233                 Allium cepa (related: onion)
    ## 234                 Allium cepa (related: onion)
    ## 235                 Allium cepa (related: onion)
    ## 236                 Allium cepa (related: onion)
    ## 237                 Allium cepa (related: onion)
    ## 238                 Allium cepa (related: onion)
    ## 239                 Allium cepa (related: onion)
    ## 240                 Allium cepa (related: onion)
    ## 241                 Allium cepa (related: onion)
    ## 242                 Oryza sativa (related: rice)
    ## 243                 Oryza sativa (related: rice)
    ## 244                 Oryza sativa (related: rice)
    ## 245 Nicotiana benthamiana (no common name found)
    ## 246 Nicotiana benthamiana (no common name found)
    ## 247 Nicotiana benthamiana (no common name found)
    ## 248 Nicotiana benthamiana (no common name found)
    ## 249 Nicotiana benthamiana (no common name found)
    ## 250 Nicotiana benthamiana (no common name found)
    ## 251 Nicotiana benthamiana (no common name found)
    ## 252 Nicotiana benthamiana (no common name found)
    ## 253 Nicotiana benthamiana (no common name found)
    ## 254 Nicotiana benthamiana (no common name found)
    ## 255 Nicotiana benthamiana (no common name found)
    ## 256 Nicotiana benthamiana (no common name found)
    ## 257 Nicotiana benthamiana (no common name found)
    ## 258 Nicotiana benthamiana (no common name found)
    ## 259 Nicotiana benthamiana (no common name found)
    ## 260 Nicotiana benthamiana (no common name found)
    ## 261 Nicotiana benthamiana (no common name found)
    ## 262 Nicotiana benthamiana (no common name found)
    ## 263       Punica granatum (related: pomegranate)
    ## 264 Nicotiana benthamiana (no common name found)
    ## 265  Nicotiana tabacum (related: common tobacco)
    ## 266     Triticum aestivum (related: bread wheat)
    ## 267     Triticum aestivum (related: bread wheat)
    ## 268                 Oryza sativa (related: rice)
    ## 269                 Oryza sativa (related: rice)
    ## 270                 Oryza sativa (related: rice)
    ## 271                 Oryza sativa (related: rice)
    ## 272                 Oryza sativa (related: rice)
    ## 273                 Oryza sativa (related: rice)
    ## 274    Phaseolus vulgaris (related: kidney bean)
    ## 275    Phaseolus vulgaris (related: kidney bean)
    ## 276 Nicotiana benthamiana (no common name found)
    ## 277       Solanum lycopersicum (related: tomato)
    ##                              Hoststrain
    ## 1                                      
    ## 2                                      
    ## 3                                      
    ## 4                            Moneymaker
    ## 5                                      
    ## 6                                      
    ## 7                                      
    ## 8                            Moneymaker
    ## 9   subsp. mays W64A-N (related: maize)
    ## 10    subsp. mays Pa91 (related: maize)
    ## 11                                Col-0
    ## 12                                CO-39
    ## 13                                     
    ## 14                                     
    ## 15                                     
    ## 16                                     
    ## 17                                     
    ## 18                                     
    ## 19                                     
    ## 20                                     
    ## 21                                     
    ## 22                                     
    ## 23                                     
    ## 24                                     
    ## 25                                     
    ## 26                                     
    ## 27                                     
    ## 28                                     
    ## 29                                     
    ## 30                                     
    ## 31                                     
    ## 32                                     
    ## 33                                     
    ## 34                                     
    ## 35                                     
    ## 36                                     
    ## 37                                     
    ## 38                              H2B-RFP
    ## 39                              H2B-RFP
    ## 40                              H2B-RFP
    ## 41                              H2B-RFP
    ## 42                              H2B-RFP
    ## 43                              H2B-RFP
    ## 44                                     
    ## 45                                     
    ## 46                     Early Cal Wonder
    ## 47                                     
    ## 48                                     
    ## 49                                     
    ## 50                                     
    ## 51                                     
    ## 52                                     
    ## 53                                     
    ## 54                                     
    ## 55                          Money maker
    ## 56                                     
    ## 57                          Money maker
    ## 58                           Maofen-802
    ## 59                            Samsun NN
    ## 60                                  A17
    ## 61              Dingras Multiple Purple
    ## 62              Dingras Multiple Purple
    ## 63              Dingras Multiple Purple
    ## 64              Dingras Multiple Purple
    ## 65              Dingras Multiple Purple
    ## 66              Dingras Multiple Purple
    ## 67              Dingras Multiple Purple
    ## 68              Dingras Multiple Purple
    ## 69              Dingras Multiple Purple
    ## 70              Dingras Multiple Purple
    ## 71              Dingras Multiple Purple
    ## 72              Dingras Multiple Purple
    ## 73              Dingras Multiple Purple
    ## 74              Dingras Multiple Purple
    ## 75              Dingras Multiple Purple
    ## 76              Dingras Multiple Purple
    ## 77              Dingras Multiple Purple
    ## 78              Dingras Multiple Purple
    ## 79              Dingras Multiple Purple
    ## 80                                  SM6
    ## 81                                  SM6
    ## 82                                  SM6
    ## 83                                  SM6
    ## 84                                  SM6
    ## 85                                  SM6
    ## 86                                  SM6
    ## 87                                  SM6
    ## 88                                  SM6
    ## 89                                  SM6
    ## 90                                  SM6
    ## 91                                  SM6
    ## 92                                  SM6
    ## 93                                  SM6
    ## 94                                  SM6
    ## 95                                  SM6
    ## 96                                  SM6
    ## 97                                  SM6
    ## 98                                  SM6
    ## 99                               Ceylan
    ## 100                              Ceylan
    ## 101                              Ceylan
    ## 102                              Ceylan
    ## 103                              Ceylan
    ## 104                              Ceylan
    ## 105                              Ceylan
    ## 106                              Ceylan
    ## 107                              Ceylan
    ## 108                              Ceylan
    ## 109                              Ceylan
    ## 110                              Ceylan
    ## 111                              Ceylan
    ## 112                              Ceylan
    ## 113                              Ceylan
    ## 114                              Ceylan
    ## 115                              Ceylan
    ## 116                              Ceylan
    ## 117                              Ceylan
    ## 118                               Surya
    ## 119                               Surya
    ## 120                               Surya
    ## 121                               Surya
    ## 122                               Surya
    ## 123                               Surya
    ## 124                               Surya
    ## 125                               Surya
    ## 126                               Surya
    ## 127                               Surya
    ## 128                               Surya
    ## 129                               Surya
    ## 130                               Surya
    ## 131                               Surya
    ## 132                               Surya
    ## 133                               Surya
    ## 134                               Surya
    ## 135                               Surya
    ## 136                               Surya
    ## 137                             AG91-25
    ## 138                             AG91-25
    ## 139                             AG91-25
    ## 140                             AG91-25
    ## 141                             AG91-25
    ## 142                             AG91-25
    ## 143                             AG91-25
    ## 144                             AG91-25
    ## 145                             AG91-25
    ## 146                             AG91-25
    ## 147                             AG91-25
    ## 148                             AG91-25
    ## 149                             AG91-25
    ## 150                             AG91-25
    ## 151                             AG91-25
    ## 152                             AG91-25
    ## 153                             AG91-25
    ## 154                             AG91-25
    ## 155                             AG91-25
    ## 156                               PM687
    ## 157                               PM687
    ## 158                               PM687
    ## 159                               PM687
    ## 160                               PM687
    ## 161                               PM687
    ## 162                               PM687
    ## 163                               PM687
    ## 164                               PM687
    ## 165                               PM687
    ## 166                               PM687
    ## 167                               PM687
    ## 168                               PM687
    ## 169                               PM687
    ## 170                               PM687
    ## 171                               PM687
    ## 172                               PM687
    ## 173                               PM687
    ## 174                               PM687
    ## 175                                 CA8
    ## 176                                 CA8
    ## 177                                 CA8
    ## 178                                 CA8
    ## 179                                 CA8
    ## 180                                 CA8
    ## 181                                 CA8
    ## 182                                 CA8
    ## 183                                 CA8
    ## 184                                 CA8
    ## 185                                 CA8
    ## 186                                 CA8
    ## 187                                 CA8
    ## 188                                 CA8
    ## 189                                 CA8
    ## 190                                 CA8
    ## 191                                 CA8
    ## 192                                 CA8
    ## 193                                 CA8
    ## 194                               PM659
    ## 195                               PM659
    ## 196                               PM659
    ## 197                               PM659
    ## 198                               PM659
    ## 199                               PM659
    ## 200                               PM659
    ## 201                               PM659
    ## 202                               PM659
    ## 203                               PM659
    ## 204                               PM659
    ## 205                               PM659
    ## 206                               PM659
    ## 207                               PM659
    ## 208                               PM659
    ## 209                               PM659
    ## 210                               PM659
    ## 211                               PM659
    ## 212                               PM659
    ## 213                          Hawaii7996
    ## 214                          Hawaii7996
    ## 215                          Hawaii7996
    ## 216                          Hawaii7996
    ## 217                          Hawaii7996
    ## 218                          Hawaii7996
    ## 219                          Hawaii7996
    ## 220                          Hawaii7996
    ## 221                          Hawaii7996
    ## 222                          Hawaii7996
    ## 223                          Hawaii7996
    ## 224                          Hawaii7996
    ## 225                          Hawaii7996
    ## 226                          Hawaii7996
    ## 227                          Hawaii7996
    ## 228                          Hawaii7996
    ## 229                          Hawaii7996
    ## 230                          Hawaii7996
    ## 231                          Hawaii7996
    ## 232                            Napoleon
    ## 233                            Napoleon
    ## 234                            Napoleon
    ## 235                            Napoleon
    ## 236                            Napoleon
    ## 237                            Napoleon
    ## 238                            Napoleon
    ## 239                            Napoleon
    ## 240                            Napoleon
    ## 241                            Napoleon
    ## 242                                CO39
    ## 243                                M201
    ## 244                       Yashiro-mochi
    ## 245                                    
    ## 246                                    
    ## 247                                    
    ## 248                                    
    ## 249                                    
    ## 250                                    
    ## 251                                    
    ## 252                                    
    ## 253                                    
    ## 254                                    
    ## 255                                    
    ## 256                                    
    ## 257                                    
    ## 258                                    
    ## 259                                    
    ## 260                                    
    ## 261                                    
    ## 262                                    
    ## 263                                    
    ## 264                                    
    ## 265                                    
    ## 266                                Gabo
    ## 267                             Fielder
    ## 268                             Azucena
    ## 269                             Azucena
    ## 270                             Azucena
    ## 271                             Azucena
    ## 272                             Azucena
    ## 273                             Azucena
    ## 274                     Canadian Wonder
    ## 275                     Canadian Wonder
    ## 276                                    
    ## 277                                    
    ##                                    GeneFunction
    ## 1                      Candidate Effector genes
    ## 2                      Candidate Effector genes
    ## 3                      Candidate Effector genes
    ## 4                              Effector protein
    ## 5                      Candidate Effector genes
    ## 6                      Candidate Effector genes
    ## 7                      Candidate Effector genes
    ## 8                              Effector protein
    ## 9                              Effector protein
    ## 10                             Effector protein
    ## 11                                RXLR effector
    ## 12      Effector (plant avirulence determinant)
    ## 13      Effector (plant avirulence determinant)
    ## 14                                     Effector
    ## 15                                     Effector
    ## 16                                     Effector
    ## 17                                     Effector
    ## 18                                     Effector
    ## 19                                     Effector
    ## 20                                     Effector
    ## 21                                     Effector
    ## 22                                     Effector
    ## 23                                     Effector
    ## 24                                     Effector
    ## 25                                     Effector
    ## 26                                     Effector
    ## 27                                     Effector
    ## 28                                     Effector
    ## 29                                     Effector
    ## 30                                     Effector
    ## 31                                     Effector
    ## 32                                     Effector
    ## 33                                     Effector
    ## 34                                     Effector
    ## 35                                     Effector
    ## 36                                     Effector
    ## 37                                     Effector
    ## 38                                CRN effectors
    ## 39                                CRN effectors
    ## 40                                CRN effectors
    ## 41                                CRN effectors
    ## 42                                CRN effectors
    ## 43                                CRN effectors
    ## 44                               GKLR effectors
    ## 45                               GKLR effectors
    ## 46                    YopJ Superfamily Effector
    ## 47                             Effector protein
    ## 48                             Effector protein
    ## 49      Effector (plant avirulence determinant)
    ## 50      Effector (plant avirulence determinant)
    ## 51      Effector (plant avirulence determinant)
    ## 52      Effector (plant avirulence determinant)
    ## 53      Effector (plant avirulence determinant)
    ## 54      Effector (plant avirulence determinant)
    ## 55                             Effector protein
    ## 56                             Effector protein
    ## 57                             Effector protein
    ## 58                                 SCR effector
    ## 59                                 SCR effector
    ## 60                            Type III Effector
    ## 61  Type III effectors (T3E) and T3E-like genes
    ## 62  Type III effectors (T3E) and T3E-like genes
    ## 63  Type III effectors (T3E) and T3E-like genes
    ## 64  Type III effectors (T3E) and T3E-like genes
    ## 65  Type III effectors (T3E) and T3E-like genes
    ## 66  Type III effectors (T3E) and T3E-like genes
    ## 67  Type III effectors (T3E) and T3E-like genes
    ## 68  Type III effectors (T3E) and T3E-like genes
    ## 69  Type III effectors (T3E) and T3E-like genes
    ## 70  Type III effectors (T3E) and T3E-like genes
    ## 71  Type III effectors (T3E) and T3E-like genes
    ## 72  Type III effectors (T3E) and T3E-like genes
    ## 73  Type III effectors (T3E) and T3E-like genes
    ## 74  Type III effectors (T3E) and T3E-like genes
    ## 75  Type III effectors (T3E) and T3E-like genes
    ## 76  Type III effectors (T3E) and T3E-like genes
    ## 77  Type III effectors (T3E) and T3E-like genes
    ## 78  Type III effectors (T3E) and T3E-like genes
    ## 79  Type III effectors (T3E) and T3E-like genes
    ## 80  Type III effectors (T3E) and T3E-like genes
    ## 81  Type III effectors (T3E) and T3E-like genes
    ## 82  Type III effectors (T3E) and T3E-like genes
    ## 83  Type III effectors (T3E) and T3E-like genes
    ## 84  Type III effectors (T3E) and T3E-like genes
    ## 85  Type III effectors (T3E) and T3E-like genes
    ## 86  Type III effectors (T3E) and T3E-like genes
    ## 87  Type III effectors (T3E) and T3E-like genes
    ## 88  Type III effectors (T3E) and T3E-like genes
    ## 89  Type III effectors (T3E) and T3E-like genes
    ## 90  Type III effectors (T3E) and T3E-like genes
    ## 91  Type III effectors (T3E) and T3E-like genes
    ## 92  Type III effectors (T3E) and T3E-like genes
    ## 93  Type III effectors (T3E) and T3E-like genes
    ## 94  Type III effectors (T3E) and T3E-like genes
    ## 95  Type III effectors (T3E) and T3E-like genes
    ## 96  Type III effectors (T3E) and T3E-like genes
    ## 97  Type III effectors (T3E) and T3E-like genes
    ## 98  Type III effectors (T3E) and T3E-like genes
    ## 99  Type III effectors (T3E) and T3E-like genes
    ## 100 Type III effectors (T3E) and T3E-like genes
    ## 101 Type III effectors (T3E) and T3E-like genes
    ## 102 Type III effectors (T3E) and T3E-like genes
    ## 103 Type III effectors (T3E) and T3E-like genes
    ## 104 Type III effectors (T3E) and T3E-like genes
    ## 105 Type III effectors (T3E) and T3E-like genes
    ## 106 Type III effectors (T3E) and T3E-like genes
    ## 107 Type III effectors (T3E) and T3E-like genes
    ## 108 Type III effectors (T3E) and T3E-like genes
    ## 109 Type III effectors (T3E) and T3E-like genes
    ## 110 Type III effectors (T3E) and T3E-like genes
    ## 111 Type III effectors (T3E) and T3E-like genes
    ## 112 Type III effectors (T3E) and T3E-like genes
    ## 113 Type III effectors (T3E) and T3E-like genes
    ## 114 Type III effectors (T3E) and T3E-like genes
    ## 115 Type III effectors (T3E) and T3E-like genes
    ## 116 Type III effectors (T3E) and T3E-like genes
    ## 117 Type III effectors (T3E) and T3E-like genes
    ## 118 Type III effectors (T3E) and T3E-like genes
    ## 119 Type III effectors (T3E) and T3E-like genes
    ## 120 Type III effectors (T3E) and T3E-like genes
    ## 121 Type III effectors (T3E) and T3E-like genes
    ## 122 Type III effectors (T3E) and T3E-like genes
    ## 123 Type III effectors (T3E) and T3E-like genes
    ## 124 Type III effectors (T3E) and T3E-like genes
    ## 125 Type III effectors (T3E) and T3E-like genes
    ## 126 Type III effectors (T3E) and T3E-like genes
    ## 127 Type III effectors (T3E) and T3E-like genes
    ## 128 Type III effectors (T3E) and T3E-like genes
    ## 129 Type III effectors (T3E) and T3E-like genes
    ## 130 Type III effectors (T3E) and T3E-like genes
    ## 131 Type III effectors (T3E) and T3E-like genes
    ## 132 Type III effectors (T3E) and T3E-like genes
    ## 133 Type III effectors (T3E) and T3E-like genes
    ## 134 Type III effectors (T3E) and T3E-like genes
    ## 135 Type III effectors (T3E) and T3E-like genes
    ## 136 Type III effectors (T3E) and T3E-like genes
    ## 137 Type III effectors (T3E) and T3E-like genes
    ## 138 Type III effectors (T3E) and T3E-like genes
    ## 139 Type III effectors (T3E) and T3E-like genes
    ## 140 Type III effectors (T3E) and T3E-like genes
    ## 141 Type III effectors (T3E) and T3E-like genes
    ## 142 Type III effectors (T3E) and T3E-like genes
    ## 143 Type III effectors (T3E) and T3E-like genes
    ## 144 Type III effectors (T3E) and T3E-like genes
    ## 145 Type III effectors (T3E) and T3E-like genes
    ## 146 Type III effectors (T3E) and T3E-like genes
    ## 147 Type III effectors (T3E) and T3E-like genes
    ## 148 Type III effectors (T3E) and T3E-like genes
    ## 149 Type III effectors (T3E) and T3E-like genes
    ## 150 Type III effectors (T3E) and T3E-like genes
    ## 151 Type III effectors (T3E) and T3E-like genes
    ## 152 Type III effectors (T3E) and T3E-like genes
    ## 153 Type III effectors (T3E) and T3E-like genes
    ## 154 Type III effectors (T3E) and T3E-like genes
    ## 155 Type III effectors (T3E) and T3E-like genes
    ## 156 Type III effectors (T3E) and T3E-like genes
    ## 157 Type III effectors (T3E) and T3E-like genes
    ## 158 Type III effectors (T3E) and T3E-like genes
    ## 159 Type III effectors (T3E) and T3E-like genes
    ## 160 Type III effectors (T3E) and T3E-like genes
    ## 161 Type III effectors (T3E) and T3E-like genes
    ## 162 Type III effectors (T3E) and T3E-like genes
    ## 163 Type III effectors (T3E) and T3E-like genes
    ## 164 Type III effectors (T3E) and T3E-like genes
    ## 165 Type III effectors (T3E) and T3E-like genes
    ## 166 Type III effectors (T3E) and T3E-like genes
    ## 167 Type III effectors (T3E) and T3E-like genes
    ## 168 Type III effectors (T3E) and T3E-like genes
    ## 169 Type III effectors (T3E) and T3E-like genes
    ## 170 Type III effectors (T3E) and T3E-like genes
    ## 171 Type III effectors (T3E) and T3E-like genes
    ## 172 Type III effectors (T3E) and T3E-like genes
    ## 173 Type III effectors (T3E) and T3E-like genes
    ## 174 Type III effectors (T3E) and T3E-like genes
    ## 175 Type III effectors (T3E) and T3E-like genes
    ## 176 Type III effectors (T3E) and T3E-like genes
    ## 177 Type III effectors (T3E) and T3E-like genes
    ## 178 Type III effectors (T3E) and T3E-like genes
    ## 179 Type III effectors (T3E) and T3E-like genes
    ## 180 Type III effectors (T3E) and T3E-like genes
    ## 181 Type III effectors (T3E) and T3E-like genes
    ## 182 Type III effectors (T3E) and T3E-like genes
    ## 183 Type III effectors (T3E) and T3E-like genes
    ## 184 Type III effectors (T3E) and T3E-like genes
    ## 185 Type III effectors (T3E) and T3E-like genes
    ## 186 Type III effectors (T3E) and T3E-like genes
    ## 187 Type III effectors (T3E) and T3E-like genes
    ## 188 Type III effectors (T3E) and T3E-like genes
    ## 189 Type III effectors (T3E) and T3E-like genes
    ## 190 Type III effectors (T3E) and T3E-like genes
    ## 191 Type III effectors (T3E) and T3E-like genes
    ## 192 Type III effectors (T3E) and T3E-like genes
    ## 193 Type III effectors (T3E) and T3E-like genes
    ## 194 Type III effectors (T3E) and T3E-like genes
    ## 195 Type III effectors (T3E) and T3E-like genes
    ## 196 Type III effectors (T3E) and T3E-like genes
    ## 197 Type III effectors (T3E) and T3E-like genes
    ## 198 Type III effectors (T3E) and T3E-like genes
    ## 199 Type III effectors (T3E) and T3E-like genes
    ## 200 Type III effectors (T3E) and T3E-like genes
    ## 201 Type III effectors (T3E) and T3E-like genes
    ## 202 Type III effectors (T3E) and T3E-like genes
    ## 203 Type III effectors (T3E) and T3E-like genes
    ## 204 Type III effectors (T3E) and T3E-like genes
    ## 205 Type III effectors (T3E) and T3E-like genes
    ## 206 Type III effectors (T3E) and T3E-like genes
    ## 207 Type III effectors (T3E) and T3E-like genes
    ## 208 Type III effectors (T3E) and T3E-like genes
    ## 209 Type III effectors (T3E) and T3E-like genes
    ## 210 Type III effectors (T3E) and T3E-like genes
    ## 211 Type III effectors (T3E) and T3E-like genes
    ## 212 Type III effectors (T3E) and T3E-like genes
    ## 213 Type III effectors (T3E) and T3E-like genes
    ## 214 Type III effectors (T3E) and T3E-like genes
    ## 215 Type III effectors (T3E) and T3E-like genes
    ## 216 Type III effectors (T3E) and T3E-like genes
    ## 217 Type III effectors (T3E) and T3E-like genes
    ## 218 Type III effectors (T3E) and T3E-like genes
    ## 219 Type III effectors (T3E) and T3E-like genes
    ## 220 Type III effectors (T3E) and T3E-like genes
    ## 221 Type III effectors (T3E) and T3E-like genes
    ## 222 Type III effectors (T3E) and T3E-like genes
    ## 223 Type III effectors (T3E) and T3E-like genes
    ## 224 Type III effectors (T3E) and T3E-like genes
    ## 225 Type III effectors (T3E) and T3E-like genes
    ## 226 Type III effectors (T3E) and T3E-like genes
    ## 227 Type III effectors (T3E) and T3E-like genes
    ## 228 Type III effectors (T3E) and T3E-like genes
    ## 229 Type III effectors (T3E) and T3E-like genes
    ## 230 Type III effectors (T3E) and T3E-like genes
    ## 231 Type III effectors (T3E) and T3E-like genes
    ## 232                                    Effector
    ## 233                                    Effector
    ## 234                                    Effector
    ## 235                                    Effector
    ## 236                                    Effector
    ## 237                                    Effector
    ## 238                                    Effector
    ## 239                                    Effector
    ## 240                                    Effector
    ## 241                                    Effector
    ## 242                                    Effector
    ## 243                                    Effector
    ## 244                                    Effector
    ## 245                                    Effector
    ## 246                                    Effector
    ## 247                                    Effector
    ## 248                                    Effector
    ## 249                                    Effector
    ## 250                                    Effector
    ## 251                                    Effector
    ## 252                                    Effector
    ## 253                                    Effector
    ## 254                                    Effector
    ## 255                                    Effector
    ## 256                                    Effector
    ## 257                                    Effector
    ## 258                                    Effector
    ## 259                                    Effector
    ## 260                                    Effector
    ## 261                                    Effector
    ## 262                            Effector protein
    ## 263                            Effector protein
    ## 264                            Effector protein
    ## 265                            Effector protein
    ## 266                            Effector protein
    ## 267                            Effector protein
    ## 268                            Effector protein
    ## 269                            Effector protein
    ## 270                            Effector protein
    ## 271                            Effector protein
    ## 272                            Effector protein
    ## 273                            Effector protein
    ## 274                            Effector protein
    ## 275                            Effector protein
    ## 276                            Effector protein
    ## 277                            Effector protein
    ##                             MutantPhenotype
    ## 1   effector (plant avirulence determinant)
    ## 2   effector (plant avirulence determinant)
    ## 3   effector (plant avirulence determinant)
    ## 4   effector (plant avirulence determinant)
    ## 5   effector (plant avirulence determinant)
    ## 6   effector (plant avirulence determinant)
    ## 7   effector (plant avirulence determinant)
    ## 8   effector (plant avirulence determinant)
    ## 9                         reduced virulence
    ## 10                        reduced virulence
    ## 11  effector (plant avirulence determinant)
    ## 12  effector (plant avirulence determinant)
    ## 13  effector (plant avirulence determinant)
    ## 14  effector (plant avirulence determinant)
    ## 15  effector (plant avirulence determinant)
    ## 16  effector (plant avirulence determinant)
    ## 17  effector (plant avirulence determinant)
    ## 18  effector (plant avirulence determinant)
    ## 19  effector (plant avirulence determinant)
    ## 20  effector (plant avirulence determinant)
    ## 21  effector (plant avirulence determinant)
    ## 22  effector (plant avirulence determinant)
    ## 23  effector (plant avirulence determinant)
    ## 24  effector (plant avirulence determinant)
    ## 25  effector (plant avirulence determinant)
    ## 26  effector (plant avirulence determinant)
    ## 27  effector (plant avirulence determinant)
    ## 28  effector (plant avirulence determinant)
    ## 29  effector (plant avirulence determinant)
    ## 30  effector (plant avirulence determinant)
    ## 31  effector (plant avirulence determinant)
    ## 32  effector (plant avirulence determinant)
    ## 33  effector (plant avirulence determinant)
    ## 34  effector (plant avirulence determinant)
    ## 35  effector (plant avirulence determinant)
    ## 36  effector (plant avirulence determinant)
    ## 37  effector (plant avirulence determinant)
    ## 38  effector (plant avirulence determinant)
    ## 39  effector (plant avirulence determinant)
    ## 40  effector (plant avirulence determinant)
    ## 41  effector (plant avirulence determinant)
    ## 42  effector (plant avirulence determinant)
    ## 43  effector (plant avirulence determinant)
    ## 44  effector (plant avirulence determinant)
    ## 45  effector (plant avirulence determinant)
    ## 46  effector (plant avirulence determinant)
    ## 47  effector (plant avirulence determinant)
    ## 48  effector (plant avirulence determinant)
    ## 49  effector (plant avirulence determinant)
    ## 50  effector (plant avirulence determinant)
    ## 51  effector (plant avirulence determinant)
    ## 52  effector (plant avirulence determinant)
    ## 53  effector (plant avirulence determinant)
    ## 54  effector (plant avirulence determinant)
    ## 55  effector (plant avirulence determinant)
    ## 56  effector (plant avirulence determinant)
    ## 57  effector (plant avirulence determinant)
    ## 58  effector (plant avirulence determinant)
    ## 59  effector (plant avirulence determinant)
    ## 60  effector (plant avirulence determinant)
    ## 61  effector (plant avirulence determinant)
    ## 62  effector (plant avirulence determinant)
    ## 63  effector (plant avirulence determinant)
    ## 64  effector (plant avirulence determinant)
    ## 65  effector (plant avirulence determinant)
    ## 66  effector (plant avirulence determinant)
    ## 67  effector (plant avirulence determinant)
    ## 68  effector (plant avirulence determinant)
    ## 69  effector (plant avirulence determinant)
    ## 70  effector (plant avirulence determinant)
    ## 71  effector (plant avirulence determinant)
    ## 72  effector (plant avirulence determinant)
    ## 73  effector (plant avirulence determinant)
    ## 74  effector (plant avirulence determinant)
    ## 75  effector (plant avirulence determinant)
    ## 76  effector (plant avirulence determinant)
    ## 77  effector (plant avirulence determinant)
    ## 78  effector (plant avirulence determinant)
    ## 79  effector (plant avirulence determinant)
    ## 80  effector (plant avirulence determinant)
    ## 81  effector (plant avirulence determinant)
    ## 82  effector (plant avirulence determinant)
    ## 83  effector (plant avirulence determinant)
    ## 84  effector (plant avirulence determinant)
    ## 85  effector (plant avirulence determinant)
    ## 86  effector (plant avirulence determinant)
    ## 87  effector (plant avirulence determinant)
    ## 88  effector (plant avirulence determinant)
    ## 89  effector (plant avirulence determinant)
    ## 90  effector (plant avirulence determinant)
    ## 91  effector (plant avirulence determinant)
    ## 92  effector (plant avirulence determinant)
    ## 93  effector (plant avirulence determinant)
    ## 94  effector (plant avirulence determinant)
    ## 95  effector (plant avirulence determinant)
    ## 96  effector (plant avirulence determinant)
    ## 97  effector (plant avirulence determinant)
    ## 98  effector (plant avirulence determinant)
    ## 99  effector (plant avirulence determinant)
    ## 100 effector (plant avirulence determinant)
    ## 101 effector (plant avirulence determinant)
    ## 102 effector (plant avirulence determinant)
    ## 103 effector (plant avirulence determinant)
    ## 104 effector (plant avirulence determinant)
    ## 105 effector (plant avirulence determinant)
    ## 106 effector (plant avirulence determinant)
    ## 107 effector (plant avirulence determinant)
    ## 108 effector (plant avirulence determinant)
    ## 109 effector (plant avirulence determinant)
    ## 110 effector (plant avirulence determinant)
    ## 111 effector (plant avirulence determinant)
    ## 112 effector (plant avirulence determinant)
    ## 113 effector (plant avirulence determinant)
    ## 114 effector (plant avirulence determinant)
    ## 115 effector (plant avirulence determinant)
    ## 116 effector (plant avirulence determinant)
    ## 117 effector (plant avirulence determinant)
    ## 118 effector (plant avirulence determinant)
    ## 119 effector (plant avirulence determinant)
    ## 120 effector (plant avirulence determinant)
    ## 121 effector (plant avirulence determinant)
    ## 122 effector (plant avirulence determinant)
    ## 123 effector (plant avirulence determinant)
    ## 124 effector (plant avirulence determinant)
    ## 125 effector (plant avirulence determinant)
    ## 126 effector (plant avirulence determinant)
    ## 127 effector (plant avirulence determinant)
    ## 128 effector (plant avirulence determinant)
    ## 129 effector (plant avirulence determinant)
    ## 130 effector (plant avirulence determinant)
    ## 131 effector (plant avirulence determinant)
    ## 132 effector (plant avirulence determinant)
    ## 133 effector (plant avirulence determinant)
    ## 134 effector (plant avirulence determinant)
    ## 135 effector (plant avirulence determinant)
    ## 136 effector (plant avirulence determinant)
    ## 137 effector (plant avirulence determinant)
    ## 138 effector (plant avirulence determinant)
    ## 139 effector (plant avirulence determinant)
    ## 140 effector (plant avirulence determinant)
    ## 141 effector (plant avirulence determinant)
    ## 142 effector (plant avirulence determinant)
    ## 143 effector (plant avirulence determinant)
    ## 144 effector (plant avirulence determinant)
    ## 145 effector (plant avirulence determinant)
    ## 146 effector (plant avirulence determinant)
    ## 147 effector (plant avirulence determinant)
    ## 148 effector (plant avirulence determinant)
    ## 149 effector (plant avirulence determinant)
    ## 150 effector (plant avirulence determinant)
    ## 151 effector (plant avirulence determinant)
    ## 152 effector (plant avirulence determinant)
    ## 153 effector (plant avirulence determinant)
    ## 154 effector (plant avirulence determinant)
    ## 155 effector (plant avirulence determinant)
    ## 156 effector (plant avirulence determinant)
    ## 157 effector (plant avirulence determinant)
    ## 158 effector (plant avirulence determinant)
    ## 159 effector (plant avirulence determinant)
    ## 160 effector (plant avirulence determinant)
    ## 161 effector (plant avirulence determinant)
    ## 162 effector (plant avirulence determinant)
    ## 163 effector (plant avirulence determinant)
    ## 164 effector (plant avirulence determinant)
    ## 165 effector (plant avirulence determinant)
    ## 166 effector (plant avirulence determinant)
    ## 167 effector (plant avirulence determinant)
    ## 168 effector (plant avirulence determinant)
    ## 169 effector (plant avirulence determinant)
    ## 170 effector (plant avirulence determinant)
    ## 171 effector (plant avirulence determinant)
    ## 172 effector (plant avirulence determinant)
    ## 173 effector (plant avirulence determinant)
    ## 174 effector (plant avirulence determinant)
    ## 175 effector (plant avirulence determinant)
    ## 176 effector (plant avirulence determinant)
    ## 177 effector (plant avirulence determinant)
    ## 178 effector (plant avirulence determinant)
    ## 179 effector (plant avirulence determinant)
    ## 180 effector (plant avirulence determinant)
    ## 181 effector (plant avirulence determinant)
    ## 182 effector (plant avirulence determinant)
    ## 183 effector (plant avirulence determinant)
    ## 184 effector (plant avirulence determinant)
    ## 185 effector (plant avirulence determinant)
    ## 186 effector (plant avirulence determinant)
    ## 187 effector (plant avirulence determinant)
    ## 188 effector (plant avirulence determinant)
    ## 189 effector (plant avirulence determinant)
    ## 190 effector (plant avirulence determinant)
    ## 191 effector (plant avirulence determinant)
    ## 192 effector (plant avirulence determinant)
    ## 193 effector (plant avirulence determinant)
    ## 194 effector (plant avirulence determinant)
    ## 195 effector (plant avirulence determinant)
    ## 196 effector (plant avirulence determinant)
    ## 197 effector (plant avirulence determinant)
    ## 198 effector (plant avirulence determinant)
    ## 199 effector (plant avirulence determinant)
    ## 200 effector (plant avirulence determinant)
    ## 201 effector (plant avirulence determinant)
    ## 202 effector (plant avirulence determinant)
    ## 203 effector (plant avirulence determinant)
    ## 204 effector (plant avirulence determinant)
    ## 205 effector (plant avirulence determinant)
    ## 206 effector (plant avirulence determinant)
    ## 207 effector (plant avirulence determinant)
    ## 208 effector (plant avirulence determinant)
    ## 209 effector (plant avirulence determinant)
    ## 210 effector (plant avirulence determinant)
    ## 211 effector (plant avirulence determinant)
    ## 212 effector (plant avirulence determinant)
    ## 213 effector (plant avirulence determinant)
    ## 214 effector (plant avirulence determinant)
    ## 215 effector (plant avirulence determinant)
    ## 216 effector (plant avirulence determinant)
    ## 217 effector (plant avirulence determinant)
    ## 218 effector (plant avirulence determinant)
    ## 219 effector (plant avirulence determinant)
    ## 220 effector (plant avirulence determinant)
    ## 221 effector (plant avirulence determinant)
    ## 222 effector (plant avirulence determinant)
    ## 223 effector (plant avirulence determinant)
    ## 224 effector (plant avirulence determinant)
    ## 225 effector (plant avirulence determinant)
    ## 226 effector (plant avirulence determinant)
    ## 227 effector (plant avirulence determinant)
    ## 228 effector (plant avirulence determinant)
    ## 229 effector (plant avirulence determinant)
    ## 230 effector (plant avirulence determinant)
    ## 231 effector (plant avirulence determinant)
    ## 232 effector (plant avirulence determinant)
    ## 233 effector (plant avirulence determinant)
    ## 234 effector (plant avirulence determinant)
    ## 235 effector (plant avirulence determinant)
    ## 236 effector (plant avirulence determinant)
    ## 237 effector (plant avirulence determinant)
    ## 238 effector (plant avirulence determinant)
    ## 239 effector (plant avirulence determinant)
    ## 240 effector (plant avirulence determinant)
    ## 241 effector (plant avirulence determinant)
    ## 242 effector (plant avirulence determinant)
    ## 243 effector (plant avirulence determinant)
    ## 244 effector (plant avirulence determinant)
    ## 245 effector (plant avirulence determinant)
    ## 246 effector (plant avirulence determinant)
    ## 247 effector (plant avirulence determinant)
    ## 248 effector (plant avirulence determinant)
    ## 249 effector (plant avirulence determinant)
    ## 250 effector (plant avirulence determinant)
    ## 251 effector (plant avirulence determinant)
    ## 252 effector (plant avirulence determinant)
    ## 253 effector (plant avirulence determinant)
    ## 254 effector (plant avirulence determinant)
    ## 255 effector (plant avirulence determinant)
    ## 256 effector (plant avirulence determinant)
    ## 257 effector (plant avirulence determinant)
    ## 258 effector (plant avirulence determinant)
    ## 259 effector (plant avirulence determinant)
    ## 260 effector (plant avirulence determinant)
    ## 261 effector (plant avirulence determinant)
    ## 262 effector (plant avirulence determinant)
    ## 263 effector (plant avirulence determinant)
    ## 264 effector (plant avirulence determinant)
    ## 265 effector (plant avirulence determinant)
    ## 266 effector (plant avirulence determinant)
    ## 267 effector (plant avirulence determinant)
    ## 268 effector (plant avirulence determinant)
    ## 269 effector (plant avirulence determinant)
    ## 270 effector (plant avirulence determinant)
    ## 271 effector (plant avirulence determinant)
    ## 272 effector (plant avirulence determinant)
    ## 273 effector (plant avirulence determinant)
    ## 274 effector (plant avirulence determinant)
    ## 275 effector (plant avirulence determinant)
    ## 276 effector (plant avirulence determinant)
    ## 277 effector (plant avirulence determinant)

    According to the email sent by Martin Urban, the data without any ProteinIDsource is due to the older PHI-base records, therefore there is no link to any database. It means that we can not get any effector data from those data without any ProteinID and ProteinIDSource. 

### Check the duplicate Phi base ID

Since we have the Phi-base fasta data, then we can get the sequences
from Phi-base fasta format data with Phi-base ID. (However, not all of
the sequences are provided)

``` r
# Check the number of unique Phi-base IDs
phibase_id_uniq <- phi_effector_plant %>% 
  dplyr::select(PHIMolConnID) %>% 
  unique() %>% 
  unlist()

# Check the number of unique Phi-base IDs
phibase_id_uniq %>% 
  length()
```

    ## [1] 638

Getting the sequence from Phi-base Fasta
----------------------------------------

Now, by using the unique those Phi-base unique IDs, we can retrieve the
sequences from Phi-base fasta data. The fasta data have been transformed
to dataframe using the function in
`0001-getting-fastaPhiBase-to-df.Rmd`. Now we can just load them.

``` r
# Load the Phi-base fasta data from RDS objects
phi_base_fasta <- readRDS("../../../../data/getting-data-new/binary-class-data/phi_base_fasta.RDS") 

# There are 7 columns or variables in Phi base fasta:
# ProteinID
# PHIMolConnID: Phi-base accession for each database entry to aid curation
# Gene: Gene name
# PathogenID: Pathogen NCBI species Taxonomy ID 
# Pathogenspecies:  Pathogen species
# MutantPhenotype: Phenotype of mutant 
# Sequence: Amino acid sequence

# Print the first ten of the data to give the idea how the data looks like
phi_base_fasta %>% 
  head(10) %>% 
  knitr::kable()
```

| ProteinID  | PHIMolConnID         | Gene                | PathogenID | Pathogenspecies               | MutantPhenotype                                 | Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|:-----------|:---------------------|:--------------------|:-----------|:------------------------------|:------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| A0A016PMX1 | PHI:6978             | Fgtom1              | 5518       | Fusarium\_graminearum         | reduced\_virulence                              | MKVTFLIVALVAAAEAVKVPCTGTTSLREGAAKKDLLIGSGAINPAYLDDAQFRAVLSQQFNSLSPENELKWNFFHTAKDTYDWHKLDRLVEFAEANDMAVKGHGLLSSCCNPDYVLNITSPDALRSEITKHFEAVMHRYRGKMDRWDVVSEALKTNGSGLASNQFYEILGPDWVEEAFRIARAADPGAKLFLNENLVESMPNKRQELYEMVPKLVSRGVPIDRIALQTHVTLEPLVPGVIRDMVNSYKALGLEVSIAELDVHTYNATLQAEIYGDVVKEALDAGITDISFWGFTDKHLYTWLPGSKPLIFNETYYPKEAFYSTHEALANFFQRS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| A0A016PTV5 | PHI:5393             | Glx                 | 5518       | Fusarium\_graminearum         | reduced\_virulence                              | MLVLRFAALAAASLPLTSASAIFRRDLTAPESPRDGWSYLGCYIDSTSKRALDGPVHYDETGLTAETCVAHCVGLGYAFAGMEYSKECFCGSKRPTTKTDEADCNMPCTGDDEQPCGAGDRLTVFGKASAGGTSVPSGSPTETVSATASASGTASPVASDLPGTWSYAGCYTDPPGRALQAAGTSKSMTPQKCIATCIADGFKIAGVEYAEECFCGNALNNAASKVKESDCNMPCSGDSSQMCGAGSRLSLYSDGDFEVNPIPVAKKDGLPGDWEYKGCVFDNNNPYLLQWLYEDAGSYATSNMTIETCLNRCQKFGYSAGGVEYGRQCVCGDLKAVENRGDVWKDDSFCSMACPGDRNSTCGAGNHINYYEWTGASLNTFHYASGPKAGKYDHFSTSPIIPLISSVGINDKIVFVEKHGTSDDDTEGSFEFDYTTNIYRELALKTDVFCSASFTLPDKAGRIINIGGWSAESVYGIRFFTPDSPQGVDNGTNVWEEDYTQLRLFDPRWYPTALVLSNGSILAMGGESGSDAPIVPTAEVLPHPAGVTESTYVDYLERAENIGRTNSYPHMAILPSGNIFFTQFNESRLLSQVDFQSIKKLPDMPGQINNPLTGRNYPLQGTLMVLPHKAPYSDPVEILICGGTTHEPGNDALDNCVLMAPDVEGAEWAIERMPSKRVMPNMVALPDGRYLILGGAQVGRGGFGLADNANLNAVMYDPEEPLGQRMTVLANTTIARLYHSEAVLLSDGKVLVSGSDPQDQGKHPQEKRIEYFWPDYLLSGATQPNFTISDRDWTYGESYTFTLTSDLEEGASKLRVSLMASVGATHGVSMGQRTLFPEFSCSGKTCSVTAPPNAFVSPPSWYQMFVLDGPTPSHAIWVRIGGDPGKLGDWPKLPGFTPPGV                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| A0A016PUA3 | PHI:3659             | Nps1                | 5518       | Fusarium\_graminearum         | unaffected\_pathogenicity\_\_reduced\_virulence | MDTRFSVQEPPVTQLRPLHTYQQSLQSSNDQLLHWTSSFEPRPDQDTAIKAFVKIVSRIIILDEGESFCIQDSQRGGYILANSDAETVKFVQYEVDQDIPTQFSIGDGKDVQFHIEDNGQTRLTVKSTTVPEAALQALGSMFDDILGQENQTGVSSWSQPSVLNFLPRTSPMDIARQDETQSNHFQSDSSKTGFLHCWFEQRALESPEAIAIDYLTDLDSGSRVQFTYRQVSNIATTLAARILELAAQSSQSLKTVAVLMGPCPELYISYLAALKTGVAFCPIAVEAPKERKEALMADLKPSAVLTTSSIFSSDTWLGTRDVLTEVVDVTPYLASPDVEPPQLSSQSTTVDNIAYVLYTSGTTGLPKGVAVSHLSASCTISSLSTHYNFTLPPAGSKPVRWFQGAAPTFDISLFEIFWTLSSGSTLCCAPRYLTLQNIDKVVTTLEADITNITPSFASLLDPKSIKGLMVGGETLNARLLQDFAQWNPPGDNEATHVPRGIFNGYGPTETAIYCIAQAHVLDNQRGSVIGTPLATCGALIVDQSSSQGLQPVPMGAVGELVITGPQVSKLGYLNRPEETALAFVDDATWGRAYKTGDRARVVFDPNGSPVIEFLGRISDDQVKLSGRRVELGEIESVLASKVDGVRETMACVWKQGAESGSEKVVSLVVVEPRGGVSFDTVHQECLEAAKLHLPDYMRPFRILQVDTLPRSASGKADRKAALRIVNDLLAKDGGVSCQTPSQGGNSHTLQPLEDSDDAKLEKELLYIVRGVLGDGSTEGLIVTPVTVLADAGIDSLRAMRILRDVRKRWPVSERPSSHSSRARLQPSLAVLLDANATVRSVFFPSTEDIETSISSSTNESDTRTKLEAFGAKHLPEAVEKLGLKESDIEMVLPATSTQSQLAVSFAMDKRNYISHTVLRLKSGVSTTALKDAVSDVLSRQAIYRSAILPCDDNLSPFIQATLTVGAWQGLVGESSVVMHKASDGRPSSDAQLWIELAEEELNLETHRLYHIQVVESEDTDNGGLLVISAAHCICDGPSLEVLMSDIARQYAGMEPLPRQGIYEAVFDWASNISADTDQLWQKALEGWETEPLGAISGNNVKPSAARARQHAMVQHSSNISWNVLEAKSRALGASPLTILQASWAMLLHILSEADTDDVTFGSVLSGHDEFIHAPTFSVVPYRVALPESQSVRQLIHNLMGYSRFAQGHRHTSFGVFKTLPYNTALALQAYPAPEYGSDGAALWTEVSNPAIRYDFAMFAEVFPTNPNSADRNGRFDDVFFKVTYRDETFSEISASCIVKQFAALTEVILRSVPDDLVQTLPARMERALLSTEGTIPTKPDVEALPADLREAYERTQVLHAQFEDQAAATPDLLALSFYSSLDAPPVELSYAELDARANGLSNILREQDIGIIPICMERSVELYVSILAIIKAGSAWCPIDTTSPVQRRTSLIARAQSKFLLTNTESLPLVEPCLEQGALEDVQIILVDKYTDCKTSIRAKPRDSIASSKVSGKDLAYLLWTSGTTGEPKGVMIQHSAAAQAMRDLQVQVEHEDSEQVRTLQLSAYSFDVFVQDLFYTWGLAGTVISGTRELVLGTFTEFIWKSQPTHAHLTPSFGASIAVNELKGSTLQYVTFIGEKLTEDVAEAWAAPDITTRAYNTYGPAENAVVSTMRRFYGKSRDNAKAANVGFPLNPCTAYVMREVGTAEGTNRWELVPRYGVGELALGGAQVGKGYLANEAKTTKAFIQGGPGIDERIYLTGDLVRLNDHGFEFLGRNDDLVKITGIRIELSEISAACASLKDDDAAIEHVETLYLQRPGAPAENNNKVVVTFVSVKAANVDTGKIRTQVFKRAKDMLPSYMVPGHVVVLDTTMPRTASNKVDRKALQAIYAESDLNVLAGRDDVAGSQGLGPASTRQQWTAEQLNILKTIASNFNIAIENLSPEDSLAGLGLSSLQVTKLAWLLRRELQCQVGVLDLMRCESLGELVDITLSRMPKTTQSQPSDTKPIETSWLAPIKAALTANIKGALRPSDTTSILPATPMQESLIVETMLEPQAYWAHREFDLSHLDEIDSRQLKRAWIAAAKRFDILRTIFVPLTQLDVECNDNTVGWAKELGVKSTILQLIRDKPTVRWTLISSEQNQSLGEVARVLQTELSPTTTIEPPWAVTFVEESNKMMLSMHHALYDAVASDMFLEAVSKFYNMEPLEDLGKVAQFDTGLDLGLLPSPAKRDEATALWNRRLEELSKTVSGGMLNAPLPDLTQSRQKQVQKILLAKQDIPASLLDAALGVSLPTLLQSAFGCVLASYLELDSIILGHTVSQRALHPDLENVVGPAIATLPLIIRSNATSAEELWKEMARDSAELFKTTHNLHPVDVKKMLNRGSGSSNAPFPGLFVYHPASDDSNDGSVAHIFNEVGQALSLNVEHPLALNVFEGNKSIELTGDGRRISQPQLELLLDQILDQARVMVEASQVPLSQLQNRMAGRLLSTSGTVTTIENTPANPTDSVALHASQHPDWVAAEEVIFQSSDDDEEIVTKTVTYAQLDTLTNAIATKLVQHEANLQPDDPVAMCLGRDIKSLAVTLAIFRAGFIYLPVDEDLPSARKQLLVRDAGAKLVITTDELVGDLGLDANDDAPVITLPDGQDDLDVIYAWPVAEHLFEAGDGGYLLYTSGSTGRPKGVRVSNSNLCHFISAFSNRLIEHSPATAKLGSAGKYLNLTSRAFDPHLTQLFVPWHLGHRVVIGKDRTSMMANLKELINELSITHFGSVPSVLTQLRLKPEDVPSVRVITTGGEKASNELLNTWAQDTKAGDDRAVLFNFYGPTEVTIGCLGHAVNHDSNARNLGLPLQGLETVLLYTSTSSDDLVVARRGQPGELCIAGPQVALGYLNRPVEDAKSFQTVSILGSTKRVYRTGDMMRMMKDGTLEFLGRADQQTKIRGQRLELDEVVGFLKQVAGDAGDLDFAAAVASSGNGSTQQQQLLGFVAKKTSGEEVELLYDYDRESGLLLDTIEQECQAKLPAFMVPTMVWVTKIPYLAASGKVDTKLLGRLASDFVARQRDEEYEDAIGQSLSTGDDLSPEESLVVAAIEESVGVSVKATSSSGIHRLGIDSLSAMTHILTASCTVSDIARLADNVASSADSTPIETPTSTDGQHDMLVNTKAHTVTDLGPLLTSLEANNIEAVLPCLPLQSALIARSLMWLSAYDSDSEGYDVPYVAQFNYRLDSNTDIARWKNAAEQVLISEAMLRTCFIQRDEDGQIFQVVLRSPPVSPLQDTRTADIVAKMHVQPPIRLLVTENSDGATVALKIHHALFDGVAIDSIKQKLEQAYNDPSSASITSIHTLNALSRISHHCNPTGEEFEVTKRSWQEELSGVQPCRIKASDDDNTSGSMVRATLCLGYTTSELNAKLKSGSDATISTSSAFQLATSLCLAQLTQQSSITYGFVMSLRPLLDDVADNVDSLVAPCLNTLIRTLALHDMGESLPELANRVHKGHTNVCQGTMPLVSVDKVQRWSGSEDKLFDSLLSINIIPARDDEEPKPGHMTALPTSSSSDMALAIDVDLHADGKIVLTLSSAGALNEKQLGDVSQLFEKALYSCSDSNVKISDIVPTLLATKAPSSITQKGSNSGIPELNTTDTDYKKALADVQETISRFLRLDLSEVSSKPESTSLYQLGLDSITVLPFVKQINKSNKTKLSSHTVIKARTIKGVAELLRDARTKQTSTMNGVPTTSVAKPEINALSEGDIYDRTLARLAKDLMFVATPLQEGMLSASLSTGGNAYTYVHAIQLSDTALDADTHNLDQFLAAVKDTVMACEILRTRFMFTDDDESPWVGIVSPTEQSDLVNWEVKSSPRGRVHLRIHHALYDATSIQAIWRTLEENYLRRLNGEVSQSPSLSKHLFRPFAKTVASSQKSSIAFWTGEVRDYTYTPIDFPTDELQASSAFHFSLEEQDLSLLHSICKDLGVTVKAALQLAWVKVLCESLYHQADVVYGEVVSTDADDDVAIVGPTINTLPMRVKLGGVSDAITVSQALNLVQKQIDYARGTNSMASLRKVQTQWRSSQQGQTSAGLFQSLFVFDGIVGSSDAANASLPLFRPAEASSGDSTGPAYDDYPLIVSFHIKDNKLNGKLRAKISATEVNAIGEQLSASLKHVLHGSDSPAIDTKHLRVAEKVEVTKPKSDTVKYTDLNGLTPTADAILNVVEKVIGTRRGGKKIGYDTKLISIGLDSILAIRLSKLLKDRVGLSVSVFEIMKGASVRDIASRTPKGGKALVKKVTSSTVDEGLKAPMAKALGIPESLIKSVIPALSGQRSHLEQWVHNGKKFFEPPWVYRVDDSFTEDKVSSAWAELVRTHDVLRTTFVCENAGELFQVTLNEQWLPLQNFTTIRESAKTIEALIHDHVSKENVKSSDMKSPAARLSFLETSDGKAVVLRLHHALYDGWSIKMVEKDMNQTLTLGKMEGSHPSLQNIVRQMRDIRDPDAETAYWKEHLSNAEETVFGTVDTATPSPFGPQFKATFPTNIPQAAVDLLSQKSQSQTSTAIILAYAKTLENITQTTRPTFGLNHASRSLSSPDGTQTLDLTDASIPTLTVTPLSVDLVGSLKSQLGFIQDHLAQLACFAQADGLSKISPRFNSYLNIIHRSDTTQSNTSDRGALQRYRLPEPLASSYFTTTEPSSTVSTVDQLSTSHLCAHRLFFNIIVSQSDGIKITLSADKAIWGGDDKMIENLVRCFRDRLEEVVKELD                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| A0A016PW58 | PHI:3658             | Nps2                | 5518       | Fusarium\_graminearum         | unaffected\_pathogenicity\_\_reduced\_virulence | MDITGLSIMNAQASRQPGPHLLHQLVKPPNQNVALDYMGSNQRVNITYHQLHEAATSLASRITKTSGSAQGQFVVPLLIHQSPSLYISLLAILKAGGAFCPLNIDAPPERVKFILDDVAATVVLVSKELASAIPNGISAAVIIVDEEEDQSSTLQSLSTEVSSRVPGPEDLAYVMYTSGSTGTPKGVGISHDAATQALIAHDRHIPSFSRFLQFAAPTFDVSVFEIFFPFFRGATLVSVRRIEMLDDLPGVLRTMEVDACELTPTVAGSLLRTRSNAPELKVLLTIGEMLNAPVVEEFGGDENRPSMLWAMYGPTEATIHCTLQAEFSSDSSTGNIGVPLDTVSCFIIEVPDSDSEQSEIRVLPQGEVGELAVGGCQLATGYINRPEQTNSVFIDSPFGRIYRTGDKARLLPNGKLECFGRLSDGQVKLRGQRLELGEVEQAVLRTSGCHSAVAAVARSILVVFCAVDAGVTEDAVLKHCGDWLPQYMVPGEVVLMSEFPRLPSGKVDRKRLKAEYEEHKEAMLEDIADSEPVDEFESSLLVVVSRVMNFKVNKSTSLAAIGMDSLSAIKLASSLRNAGYSIDTTDLLTAKTVFDIIFAARRQIQNQTVSSLPFSTNLSLDFGQILQQNATLADMSGLIEEIIPCTPLQAAMLAETSHNSTAYCNQVELAIPLSYSAYQISESFAQLSQKNPILRTGFAIVDRRFVTLVYGELRPEQIKVVDQVQGDFSLSSPEDFCSPMRLQIQRDSVDEKSRVLLQIHHSLYDGWSMDVLLSDWSKLLLQEPVSEHSSFREVVKFYQQLQTSDDARMFWTENLAGWKQTPLPKLRDKLVHSGEVLSFRRPMSLSRRRVTAEVQRNGFSPQVLFQASLALLWTSVTGARDITIGSVTSGRTIPVIDIEQIIGPCIAALPVRIDFDKISIGLELLKNLHSSNRKIMQYCTVSLSEVKKLVGLQLGESLYDVLFVYQESLASSERTQCMVKEATHLDRLETPILFEVEPTEDGFTLQVTYHEAIVPPATVQHMVDQFEALAHSILERPTQEIKCALREIKCTPSVDNINASPPQRVPDLARLFEDVVRKHPEENALLFYHSLKVANNVVWSFRELNNEANQIAHYLQSCGIQVGQVVAVIMEKSPALYASILAIIKCGCGYLPILPSTPLARTREILLQAEIKYCLVDSSPDQLASMLELSTITVNTNLFNEFSTANLDNEVDGSRLAYVIYTSGTTGTPKGVAVQQQSIVANIEHLEATYPKPSSSQGRLLQACSQAFDVSVFEIFYTWCAGMCLCAGTNDTILEDIERSIRDLEITHLSMTPTVAALVEPSNVPSVQFLVTAGEPMTQAVHNKWCHQLWQGYGPSETTNICTVKKMATDDHIEHLGHVFPNTSAVVLSPATLDTVPLNWVGEFCFGGAQVAQGYLNMPELTAQKFIHHPQYGKLYRSGDMGRMLPDGSLVILGRIDDQVKLRGQRIEIGEINSTVTMAGFATSAATVLVEHEESTIKQLALFFVPQHDPTEFRVLEIDNEVQQSLAAHMQSRLPGYMVPSYLVPISSMPMTSSGKVDKRRLHDCFDKLDRHYLEKVSRSSGDNPDEGDWSRMDLVISEVIKESVAASAGGFGRWTPFTVLGVDSISAIDLARALNAKLGARVAVSDILRNPTIAQLAKHLEGKPLYEETAFEGTQREFFPAAFTAAIKEVFLGESKAIKDILPCTPLQEAMLSRGQRGYYNKVLLRLKAEPGAIRSYWEVMSKRHDILRTCFATTTDSKHAIAQVVLEDWEIPWRTFDISEPSFDGAIEEHLKSLPDPVDSRTPPVSLALLRYRGSAFLSFICHHALYDGVAMERLLKEVEALAGGEDLLPPVSYKEFLEISTNLPNDTEEFWQQHLRGYKALSIFTQSSSSEIDQSTCTTSLDMPLANLQGRLRDFGTTLLSVCQASWATVLAMTYRQPDVCFGNVMSGRTLDIDGLERLVAPCFNTIPIRVALPTTSSNIDMVKHLQKLNTEMLTYQFTPLRLIQRSINRTGKHIFDTLLLLQKPLQDIDQTVWELEADSGDMDIPLVCEVVPCPGLNSLVINLHRDMSIVTEDVASAMADAFKVILKAILTAPHSTPMTAEDLPDSLRSILQQLKPQYDKKDNTGNIPDGEEEWSEVELDVRQVLAKLSGVSEQQIKRRTTIFQLGLDSINAVQVASILRQRGFIVSASDVIECPSCSKIAAKLLENRSRTKSEDLKRYDIGRFSHQVYSEIAGRLPQTATIEAVLPCTPLQSAMLASFIQSGGENYLNAMEYIVTDEISLESLTKAWQLLHERHPMLRTGFVPVQHPDATFAMVRYAPGSMKTPVSIAESEGDEVSDLLNLKGNTSERVLTALHQPPWTVVLAQTPQQTSMKFVAHHVLYDAHALQMMLYDLSRLVKSERLPPVSRIEPAVSAILINSLDEQGSEKAFWEAKASSTVVNKFPLMTPLRVESRCMLADTTVASLSFTKLKRATQASNVTIQAVIQAAWTRVLASYLGENSVVFGVALSGRTTDETKDAPFPCLNTVPVVGNNVTSNAELVSYMMEYNQRLHKHQFSPLGKVQRWLGHPTGPVFDTLIAYQKMADAGSSSLPWKLVKDEARVEYSVSLEIEPTENDHVRLCITYYNDILPREQAHLLMKQFDSSLNHIACNHLATEDEAFNLSPDLYSVLPPSHPVLDSPVQFLHQFVELGAAVHPNKLALEFVSAFDGDTCLKKQWDYRQLNIMGNRVANMLQEKLTPGSIVAIHFDKCPEAYFSILGILKAGCSFVALDPSAPKARKQFIVEDSRAPCLLTRSLEDLDFEAKTAILEVKVESLSVLEEEELIFQPAISPSDTCYCLYTSGTTGTPKGCEITHDNAVQAMMAFQELFKGHWDADSRWLQFAALHFDVSVLEQYWSWSVGMAVVAAPKDLILDDLTASINKLEITHIDLTPSLARLTHPDEIPSLCRGVFITGGEQLKQEILDVWGPKAVIYNAYGPTEATIGVTMFQRVPVNGRPSNIGKQFPNVGSFIFKQNTNTPVLRGAVGELCVSGRLVGKGYLNRPQLTEERFPTLEEFGERVYRTGDLVRVLHDGCFDFLGRADDQVKLRGQRLEIAEINHIIRTDVTEVHDAATIVARHGTSGKDVLVSFIVSEHLTTGPLRVVSDDEGLATKAKEACRAKLPGYMVPTYILLLPYIPLSSNNKAEIKDLKKLFSELAPEQLMELSHAATAPVSRGAQDILVLLYDALAQFSNISKDDISPTTSIFDVGVDSITALKLASLLKSRGLHAVSPAMLLKNPVIGDLANCLAKAASSQRQKLAREIKQSIQAYAHRHRGLVYSSLNIGPADIEYIAPCSPLQEGIISRSLTSTKPGAYFNTFQLKLHQSTTTTKFQQAWKDLVFSESILRTVFVPSTDGFLQVALRNPLFPWESTAFGSNELAEYYFAEQKENWIQRNKSSITQPLLLTYVETPTSRLFTVHIFHALYDGNSFDLMMDRVAANYAGTSVQKAPSFFEALTSGPLTRHDNCKGFWEKHLEGWVPSSITAHERSTRGSVVVAERDMPISNFEAMRSSHNVTLQAVIMALWTSVLQNLVESQITIGVVVSGRAVDLPGVENTIGPLFNTVPFFRQAVQHEDWKSLVRRCHDFNASVLDFQHVPLKNIQKWCSHNKALFDTLFTYQIDEAKTDDNELPFEIQNSEVTPDYPLALEAVYTKTGKLRFTLVAQGHVVSQSILDNLLNEIERFADLAAESPQSEVPVPQFKIPIVDDFHSGNAEKDSQNSFEWTSEAQAIQNEISVLVGINPAEIAQDVSILELGLDSIDVIKLAAKLKRKSINLAPSQIMRQQTIAKMITELASITNDSSCPPRDNFLSRIGYRLREHLEASEVDISNVESVLPPTHLQESMVARMIHSGFEAYFNHDVLRVSDHVDTTQLIQAWKDLIHQTPVLRTGFYQVESQDFDMTFCQVVSKSIDIDFEATRVEDLSELHQITDAAKSAAKNGRGQKKLFQLKLVVIGQERYMVLSIAHALYDGWSLSLLFQDLQALLEGRLITRPPVEQFIARVMESTTSKAKDFWMQYLQDAPSSTILTKVQLPTVEEKVQRFESVSKVLCQACWAVTLARQIRSTDVTFGTVLSGRDFDGADSLVFPTMNTVALRCILHGSAAEFLRYLEENMTDIRDFQHYPLRKAQSAAKVDGQDLFNTLFILQRSPVSSDPADRPLLTSVEASSSTEYPLCVEAEAVSDSLVWRLALQPQCSWNGGPQSLLETLDNVMSFLLKSKDPEILSFSERGVSICGTPPVALPESIIHEDASDNYSSRDEKIEWNQNEIGIREVLFQVSNVPILSIKLSDNLYHLGLDSISAIKVSSLLRKAGINLRPQDLIKSSSISEMAQKADTKLKKPLQTLETVEDWLPPADIDVNKLLADNGINKDEAEVLPALPMQVYMLTAWENSDGSVFFPEFPCRIKTSASLGEIDQAWGKLVSETPLLRTCFASTQSSTIPFIQIILKELGIPLSSLQPGERSNCCIRPLVEVNIEQEDKDTWLLRLKLHHALYDGVSLPALLQRLSELLNGSGTMENKGLSQWKQFTMRHTTDEARIARREFWTSYLKGSSSSPIIANSDTDVKMRTSHLNRSAISDISSIQALATQSGVSFQSLFLSAYARALAKQNNVSDTVFGLYLANRAAGENLPQTYPTLNLVPLRVSSPINRPLAAVAADIQRDIHLITSESRAEVGLWEIAQWTGIRITSFVNFLTLPGDTDPTGNSITVLPETNTGVVVKDNLPDRPRTPYLESIFRNDIPMAIDIEASVDGKNLAIGVFGSLQQISSEEALTLVANIAEILDGGI                                                                                                                                                                                                                                                                                                                                                                                                  |
| A0A016PX00 | PHI:4602             | Fdb2                | 5518       | Fusarium\_graminearum         | reduced\_virulence                              | MYGAGQGAQTGISTPRSSASLRPLTLSHGSLETSFLIPTNLHFHASQLKDKFVATLPEATDELAQDDEPSSVPDLVARYLDLIAHEVDEGEDDDAGSYEEVLKLVLNEFERNFLRGNEAHAVAANLPGIESKKLDFIRSYYTARSVCNRPIKPHASALFRAAEDGDAEIYTIFGGQGNIEEYFEELREIFKTYSSFVGDLITRSAELLQTLSKDPKAEKMFPKGLDIMNWLHHEDSTPDVDYLISAPVSFPLIGLVQLAHYEVTCKVLGVHPGMLRERITGSTGHSQGIVMAAATAAADSWDSWQDISSSVLTMLFWIGTRSQQAFPTTSMTPTMLRESSEHGEGSPTPMLSIRDLPQAEVQKHIDATNQYLPEDRHISISLINSPRNLVVTGPPTSLYGLNTQLRKVKAPVGLDQNRIPFTERKVRFANRFLPITAPFHSKYLSEATAMIDEDLKDINIDSSDLGIAVFDTNTGKDLREEVKGNIIPALVRLITRDPVNWEKATVFPDATHILDFGPGGVSGLGVLTSRNKEGTGVRVILAGTVDGGMNDLGYKAELFDRDEENAVKYAIDWVKEFGPKLVNNKSGRTYLDTKMSRLLGLPPVLVAGMTPCTVPWDFVAATMNAGYHIELAGGGYFVGPMMTDAITKIEKAIPAGRGISINLIYVNPRAMAWQIPLIQKLRSEGVPIEGLTIGAGVPSIEVAQEYIETLGLKHISFKPGSLDAIQSVINIAKANPHFPVMLQWTGGRGGGHHSFEDFHHPILQMYGRIRRQENIILVAGSGFGGADDTYPYITGDWSKKYGYPPMPFDGCLFGSRMMVAKEAHTSKAAKEAIVAAPGLEDSEWEQTYKGPAGGVLTVRSEMGEPIHKLATRGVRFWAEMDQKIFSLPKEKRVAELKKNKDYYIKKLNDDCQKVWFGRNKEGKAVDLDDMTYAEVLRRLVELLYVKHQSRWIDRSYTVLVGDYIHRLEERLTATPGQASLLQSYSELSDPFEVIDRILAAYPDAETQIINAQDVQYFLIMCMRPTQKPVTFIPVFDDNFEFYFKKDSLWQSEDLDAVVDQDVGRTCILQGPAAAKFSTEIDQPIKSILDDIHETHVKYLTRDQYNDNAASIPYIEYFGGKLVDPEIPLDVEGLTVSYDTYKNTYRLSSSAAATLPELDSWLSLLAGPNRNWRHALLMADVVAQGQKFQTNPIRRIFAPSRGLFVEISYPNDPSKTTIVVREQPRHNQYVDVIEVKLVGENKIQVNMIKETTALGKPAALPLEFTYHPEAGYSPIREVMEGRNDRIKEFYWKAWFGEEALDLDADVTAKFDGGKTTITGEAINDFVHAVGNTGEAFVDRPGKTVYAPMDFAIVIGWKAITKPIFPRTIDGDLLKLVHLSNQFRMIPGAEPLKKGDEVSTTAQINAVINQEAGKMVEVCGTLVRDGKPVMEVTSQFLYRGTYTDYENTFQRKVETPIQVHLATTRDVAILRSKDWFSVDELPQNIELLGQTLTFRLQSLVRYKNKTVFSSIETRGQVLLELPTKEIIQVGSVDYETGESHGNPVIDYLERNGQPIDQPINFENAIPLSGRSPLALRAPASNENYARVSGDYNPIHVSRVFSSYANLPGTITHGMYSSASVRSLVETWAAENNIGRVRSFHASLVGMVLPNDDLEVKLQHSGMVSGRKIIKVEVRNKETEDKVLEGEAEVEQPATAYVFTGQGSQEQGMGMELYASSPVAKEVWDRADKYLLDAYGFSITNIVKNNPKELTIHFGGPRGKAIRQNYIAMTFETVAADGSIKSEKIFKEIDETTTSYTYRSPTGLLSATQFTQPALTLMEKASFEDMKAKGLVPNNSTFAGHSLGEYSALAALAEVMPIESLVSVVFYRGLTMQVAVERDAAGRSNYSMCAVNPSRISKTFNEAALQFVVDNIAEETGWLLEIVNYNIANQQYVCAGDLRALDTLAGVTNFIKIQKIDIEEVKDNIEEVKGHLREIIRGCAEKTLAKPTPLELERGFATIPLRGIDVPFHSTFLRSGVKPFRSFLLKKINKTSIDPSKLVGKYIPNVTAKPFALTKEYFEDVYKLTNSPKIGAVLANWDKYNQDGVAANGVEGSESSSGEYEGSGRAA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| A0A016Q4W5 | PHI:3657             | Nps6                | 5518       | Fusarium\_graminearum         | reduced\_virulence                              | MTELDLSAEPNMPQTKTVRHNVPVLRKHHEQDAQELEILLLAWSLLLYRHNHGNHVEFSWGLTKIESSTCRTFTLNTAKLQWDGSNIVASELEVFKTYTQQQLQSEVPSKRDQYKLFFNDEPATGDLVNYINEDGDISVSWGNVQIQATLEDDALCLRPTWREPLGAEFLANHLAQAFVEVLNTTLADPDATLSSVLPLGKLDKSVIWNWNKDLPPPVPECIHTLISEQVRLRPDAQAICSWDGDLTYAEMDNLSTLLAQHLINLGVKNGDIVPLCFEKSRWTTVGVMGVIKAGAAFVLMDPSQPIQRRQVMAQQVKATHILTSRDQAKYGPEIAPEAKHVIVDTETLDSLAKTIEDPSRELPQVPPESLLYIIFTSGSTGTPKGVMLSHETYTASALARSTGIGYSSISRSLDFTSYAFDVSIDSILCTLIRGGCLCIPTDMDRVNDLSGAIRRLKVNMVNITPSVARILDPDIIPSLNSLGIGGEACSAGDIAIWGQHTRIVIGYGPAECTIGCTVNPSAAGKPYVSMGPGTGACIWLVDPDDHNKLVPVGAVGELLIEGPIVGQGYLGDPEKTKEAFISDPDFLLAGADGIPGRQGRLYKSGDLVRYDPDGENGFIFVGRKDTQIKLRGQRVELGEIEHHIKNLLPSGAEVVAEIIAPRNQNKESMLVAFVADREAKDEGDARQIDFPPRFREALEVLNDKLSKVVPVYMVPTQYITLSKIPYLVSGKTDRKSLRALGAEISANMQASAAANESSEIREPQSEAELFLRDSWCRLLGLENKQVSTTHNFFTSGGDSVLAMKLVPIVRDWGYTLSVADIFNYPVLSDMANSMQKGDSSKGVDMQIPEFSLLKDDMDREALCAEAAHNSGCDVSAIEDIYPCAPMQEIHMAFYTRSKENYVAQRIADIPASSSIDKLKSAWNAVYKESPILRTRIVEFKQHGFMQVVVNEALQWEEVDSSLEDFIEKDKKEPMSPGAPLSRFAIVTDKALDKRYFVWTAHHAIYDGWSTDLIVEHARAAYKGQEVSRPAQFKHFIRYLAEDSRESSKAYWKTQLAGATGPQFPSLPSRSYIPDPTSLTERFIKLDKAAKSDITIATVIRAAWALLASQYSMRDDVVFGETFMGRTIPLPGAELIEGPILATVPVHIRLDRTTTVQDFLRAVQEQSVKRAAHEHLGIQHIRRLSDDAQIACEVTMGLVVQPQDPDPTETENDALPSFRGGDAALEALHFNSYPLMLAVSMQKTGFRLLASFDSELLSHVQVERVLSQFEVAINQLRGDMSRSLNELTCLGEEELGQIWEANKHAPVPPKDISRFLSTGDKYPTVQYVPWVVQPGNEKLLMPLGSTGELLLEGIASAGDDDDVVDAPEWLKEGAMGYPGRQGKLVRTGDLVKYADDMSLVFIGRKDAMTSVDGRVVDLNATNIELNRLLPSNTKAVSRLVVPKGSNSQTPVVVALVQETSAKDTQLLKLGLDISDASMPLSQAVSVELATAIIGLNKAMVETLPPYAIPSICIPLDNSADIDSIDTIVENITLSLVVELRKSFASLKKTIADTTVLTTKERVLRASWSKFLGIEEEKLALDDNFFRLGGDSIVAMRMVSALRQDGYRLSVASIFQNMQLRGMADSLAEISTESVEAAKEYTAFSMLDTNDVDSFLAKVIRPQLADPKWKIKDVLHATGPQAGDVKQSVFAPRSSVQYNMLYLDQSIDTARLVDSFQYLVSQHAILRTVFVENEGQTLQVVLDDLKVPVTEQNVEGSIDAAAKQLAVSDINSSIDLTHGSSFIRLFVLRGESENALVIRISHAQYDGVSLPELLRQLELRYRGLEIPSSEPFETYIQHLSAAKPQNVEFWRKTLQGSSYTEIAPAADPQKKTAFMTKDVDISGASPNTTHAMLLTAGWAKVLSQHLNVSDVTFGGIVSGRDVDVAGIDTIMGPCYQYQPVRVKFESNWTATDLLDSVRSQSLEGSQRATLSFQEVLKECTNWPAGTPFYGSFTNHLNKEFFDSIPFAGTKCRVDYSIPHPEPATPPRVVSFLEDGRTQIGIEADEERREFWEARLGELAKVIEGFVKNPQALI                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| A0A023H5D8 | PHI:6442             | Eepr                | 615        | Serratia\_marcescens          | reduced\_virulence                              | MDNNHQKFDSQSIANRVRELFLHYGIGKRQHARELSRILDLSFSHAHRKLKGQSPWTLEQINSVAAALGETPAAIADLSAEHETTEPNMARDAIFFVAGVAMPCVGHIGDELPAGRPAEFVALRVEGQWHIYRADEAPAGPRYGVELIEIRPGYGDDERLSIAVLDDSHQAADELAKYLGDCGFNAVAFYDVDSFCQALQQSLFDGYVVDWLIGEETADRCIATIRASDNPDAPVLVLTGELGTDRRESEIAQAMREYDVLGPYEKPVRLHVIEAALQRCFNL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| A0A023NA98 | PHI:3354             | Rtxa1               | 672        | Vibrio\_vulnificus            | reduced\_virulence                              | MGKPFWRSVEYFFTGNYSADDGNNSIVAIGFGGEIHAYGGDDHVTVGSIGAKVYTGSGNDTVVGGSAYLRVEDTTGHLSVKGAAGYADINKSGDGNVSFAGAAGGVSIDHLGNHGDVNYGGAAAYNGITRKGLSGNVTFKGAGGYNALWHETNQGNLSFAGAGAGNKLDRTWFNRYQGSRGDVTFDGAGAANSISSRVETGNITFRGAGADNHLVRKGKVGDITLQGAGASNRIERTRQAEDVYAQTRGNIRFEGVGGYNSLYSDVAHGDIHFSGGGAYNTITRKGSGSSFDAQGMEYAKAEEIVLTAAQMHGLSIDNGNKFHAVTAVKSEREPNTYLFAIADGTYTKINKVRLYNDPETGKLKYYSEAWFKRGNHLAELARSDVSSAGGFEVNPINGGYTLSNIAVEHQQSLTVHAVEKDLTEYEWVTYANGVLIDAKDVALSEAKMGGHAISTDGTTVDVQAVKSNRKPNTYVYAKVLGPYTKIVVVELANDPKTGALKYQARSWYKEGNHTANLANEDISSANGYHSMGKSGYSLSDLHYSVNAVRSTSETVADIDEYTDQTLFKPATDSGESSGDVRFNGAGGGNVIKSNVTRGNVYFNGGGIANVILHSSQFGHTEFNGGGAANVIVKSGEEGDLTFRGAGLANVLVHQSKQGKMDVYAGGAVNVLVRIGDGQYLAHLLAYGNISVHKGNGNSRVVMLGGYNTHTQIGSGNGLWLAAGGFNVMTQVGKGDVASVLAGGANVLTKVGDGDLTAGMLGGANVITHISSDNETSNTTAVALGGANILTKKGKGNTVAVMGGGANVLTHVGDGTTTGVMVGGANILTKVGNGDTTGIMLGVGNVLTHVGDGQTLGVMGAAGNIFTKVGDGTSIAVMIGAGNIFTHVGEGNAWALMGGLGNVFTKVGNGDALALMVAEANVFTHIGDGMSVALMLAKGNVATKVGNGTTLAAMVGNANIFTHLGSGSTFAAMIGQANIMTKVGNDLTAALMVGKANIYTHVGDGTSLGIFAGEVNVMTKVGNGTTLAAMFGKANIMTHVGDGLTGVLALGEANIVTKVGDDFMGVVAAAKANVVTHVGDATTAAVLAGKGNILTKVGEGTTVGLLISDIGNVMTHVGDGTTIGIAKGKANIITKVGDGLGVNVAWGQANVFTQVGDGDRYNFAKGEANIITKVGDGQEVSVVQGKANIITHVGNGDDYTGAWGKANVITKVGDGRNVVLAKGEANIVTQVGDGDSFNALWSKGNIVTKVGDGMQVTAAKGKANITTTVGDGLSVTAAYGDANINTKVGDGVSVNVAWGKYNINTKVGDGLNVAVMKGKANANIHVGDGLNINASYAQNNVAIKVGNGDFYSLAVASSNTSSNKLSALFDNIKQTLLGVGGSQAINYLVQGDEASSSGTQKGRGAIATPEITKLDGFQMEAIEEVGSDLGDSLTGSVTKVDTPDLNKVQNALDVDGSSDQTQAPNLIVNGDFEQGDQGWKSTHGVEASHSGNVYGVNGEGHGARVTELDTYTNTSLYQDLTDLTEGEVIAVSFDFAKRAGLSNNEGIEVLWNGEVVFSSSGDASAWQQKTLKLTAHAGSNRIEFKGTGHNDGLGYILDNVVAKSESSLQAKAVSEHATQNQASQNALSDKERAEADRQRLEQEKQKQLDAVAGSQSQLESTDQQALENNGQAQRDAVKEESEAVTAELTKLAQGLDVLDGQATHTGESGDQWRNDFAGGLLDGVQRQLDDAKQLANDKIAAAKQTHADNQNKVKDAVAKSEAGVAKGEQNRAGAEQDIADAKADAEKRKADALAKGKDAQQAESDAHHAVNNAQSRGERDVQLAENKANQAQADAQGAKQNEGDRPDRQGVAGSGLSGNAHSVEGAGETGSHVNADSSTNADGRFSEGLSEQEQEALEGATNAVNRLQINAGIRGKNSGSTITSMFTETNSDSIVVPTTASQDVVRKEIRISGVNLEGLGEASHDSAESLVAARAEKVANLYRWLDTENDVATDKYVPVPGFERVDADVSDEVKQRMIQSMSGYIEHTDNQVPKDQAEALATLFVESTLDYDWDKRVEFLTKLESYGYSFEAPHAEKSIVSFWSGKNFKQYRDVLDNAQTDGKKVVYDIDVKGNAFAMDLNKHLMRWGGLFLDPDNAEQNQLKSSIDAATFSNTGFWSSVYATGAQNDVYVIAEGGVRLGNYFWNVDLPALRQLQREGLVGEIRLLDKPVSEYKDLPADQIGRRLTDAGVAVKVRFDALSHERQAELLADNPDGYKADTLVELDVKLSAIDSMLRESLPFYSLRTERNLLVQEGEEGFEVRSWPGTDGKSKTILLDNPEDAAQQKSIERFILANFDNFEQMPDELFLVDNKVLSHHDGRTRILAQKEDGAWTYNTNVELMSVTELLDAAHVSGKVRGESYQQVIDALTEYHASTAEHADYELTSVEKLLNLRKQVEGYVLGHPDSGRVQAMNSLLNQVNSRLEAVSVLVVSEQSIKAHDSFSHLYDQLDNANLKESKHLYLDGNGDFVTKGKGNLANIDKLGGSDAVLEKVKAAVTHEYGQVVADTIFAGLSANDLAKDGKGIDIAGLNKVHQAIEQHMSPVSATMYIWKPSKHSTLGHAALQIGQGRTQLEGQAAADFNKQNYVSWWPLGSKSPNIGNILNVATKDQPDLKLRWSDFSQPAHQSDTLEHDMAAEENDSFGLKKGEAKLKRFIEELNAAKGIDASFKMASEGYASLLLGNPDMLASTGIPAHVFQPFLDQWNDTSYDMMDVANRFAEELQKQAKIEVNPEQIEQQISEVVKEFAQDELDKIQAFKVAQADQGRVFRINLEGLDVAAMQAEWHRLSNDPDARYQLLTKNCSSTVAKVLKAGGADKLIGHTWLPKFGVWTPTELFNFGQALQEAQLEIAAKKQSHQVTDVLDALSGNEKHKENVAIENDGTPPRDKESLSPLTRFLNNELYGEKDARRKIGEITQTLLDHAVENGESQKVTLKGEAGRLTGYYHQGAASSEGETSATSGKVVLFLHGSGSSAEEQASEIRNHYQKQGIDMLAVNLRGYGESDGGPSEKGLYQDARTMFNYLVNDKGIDPSNIIIHGYSMGGPIAADLARYAAQNGQAVSGLLLDRPMPSMTKAITAHEVANPAGIVGAIAKAVNGQFSVEKNLKGLPKETPILLLTDNEGLGEEGEKLRAKLAIAGYNVSGEQTFYGHEASNRLMGQYADQIVSGLFNAEQAAVEVKDIRATEDLSVVKTVASDTELGTNTDAPHKNYQSRDLVLEPIVQPETIELGMPDSDQKILAEVAERENVIIGVRPVDEKSKSLIDSKLYSSKGLFVKAKSSDWGPMSGFIPVDQAFAKASARRDLDKFNGYAEQSIESGNAVSADLYLNQVRIDELVSKYQSLTALEFDAESGMYKTTATNGDQTVTFFLNKVTVDSKDLWQVHYIKDGKLAPFKVIGDPVSKQPMTADYDLLTVMYSYSDLGPQDKLKQPLTWEQWKESVTYEELTPKYKELYNSEVLYNKKDGASLGVVSDRLKALKDVINTSLGRTDGLEMVHHGADDANPYAVMADNFPATFFVPKSFFMEDGLGEGKGSIQTYFNVNEQGAVVIRDPQEFSNFQQVAINVSYRASLNDKWNVGLDDPLFTPKSKLSHDFLNAKEEVIKKLSGEVETNVRTTQLLTDNEGLGNEGEKLRTKPTASGFFESSDVQPEVMKGLLQNVGDKIFDLKAGGKELDMFSFHFSQSADLLNKLITLIPEVGNVLITNGNDVQSKDFLEGVCFALSARYMMEERVHGLGGGKAYMEWLKDTVQAYNDNITNKKNDIGSVEQKLLNQYRRQNLGLAIKDLLSMQYSQLMDTSTSAARDAANKEYSGKLRANGLTGPNINDALNYGADGYESVMDKLRNVNKSTYMTFMSQGHAMSVVVHKKGNHKVWSFYDPNFGTKSFAQYDDFRGFMDNFHKGLLTQYKFQDSEEAGQSFYVRFKKFEEGDISSYDGLWKNAREGEQSYVLRALKEQGKTFSMGKNITGKLVDFNDDVITLDVTSKNGRKVLVEVAVSDISQAANLVKTNISQVFSDPLASKLSIQSHAESATITVLEVSGQEAISEVVEGAKIPGQKDAWTGATSKADNQNVNDWERVVVTPAVDGGETRFDGQIIVQMENDAVAAKAAANLAGKHPESSVVVQLDSDGNYRVVYGDPSKLDGKIRWQLVGHGRDHSESNNTRLSGYSADELAVKLASFQQMFNQAEKISSKPDHISIVGCSLVSDDKQKGFGHQFINAMDANGLRVDVSVRSAKVYINEMGRKLYFDGKDSWVNKAINSKVLLSWNGQGDVVAKDERIRNGIAEGDIDLSRIGISDVDEPARGAIGDNKDVFDAPEKRKAETETSSSSANNKLSYSGNIQVNVGDGEFTAVNWGTSNVGIKVGTGGFKSLAFGDNNVMVHIGNGESKHSFDIGGYQALEGAQMFIGNRNVSFNKGRSNDLIVMMDKSIPTPPLVNPFDGAARISGVLQSIATSGEGQDWLAAQEQQWTLSGAKKFVKDMSGLDQSSSVDYTSLVELDSQNERSSRGLKHDAEAALNKQYNQWLSGNGDSDTSKLSRADKLRQANEKLAFNFAVGGQGADIQVTTGNWNFMFGDNIQSILDTNLGSLFGLMTQQFSATGQAKTTFTYTPEDLPRQLKNKLLGQLAGVGAETTLADIFGVDYTASGQIVSRNGEAVDGVSILKEMLEVIGEFSGDQLQAFVDPAKLLDSLKSGINMGADGIKSFAETHGLKEKAPEEEEDNSSVSVNGASVNSAQGATVADGSTETAETPDRAFGFNSLNLPNLFATIFSQDKQKEMKSLVENLKENLTADLLNMKEKTFDFLRNSGHLQGDGDINISLGNYNFNWGGDGKDLGAYLGDNNNFWGGRGDDVFYATGTSNIFTGGEGNDMGVLMGRENMMFGGDGNDTAVVAGRINHVFLGAGDDQSFVFGEGGEIDTGSGRDYVVTSGNFNRVDTGDDQDYSVTIGNNNQVELGAGNDFANVFGNYNRINASAGNDVVKLMGYHAVLNGGEGEDHLIAAAISKFSQFNGGEGRDLMVLGGYQNTFKGGTDVDSFVVSGDVIDNLVEDIRSEDNIVFNGIDWQKLWFERSGYDLKLSILRDPASDSDQAKFEHIGSVTFSDYFNGNRAQVIIAMGEKDATGEREYTTLSESAIDALVQAMSGFDPQAGDNGFIDNLDSKSRVAITTAWADVVHKKGITV |
| A0A023UJQ9 | PHI:5538\_\_PHI:5557 | Avr5(cfce1)\_\_Avr5 | 5499       | Passalora\_fulva              | effector\_(plant\_avirulence\_determinant)      | MKSPIVITILATALGALGSYDALPINCRDTTNYCFNGNGRHEVCSYCNQAKEEPLKLGRRGGQRDCGVAGSQCNDVDHQQCDARCCSKIGSPTFYGVRCPYPY                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| A0A023Y9U3 | PHI:3462             | Rpff                | 40324      | Stenotrophomonas\_maltophilia | reduced\_virulence\_\_unaffected\_pathogenicity | MSAVRPIITRPSQHPTLRITEEPERDVYWIHMHANLVNQPGRPCFASRLVDDIVDYQRELGDRLSASHTLSPHVVLASDSDVFNLGGDLELFCRLIREGDRARLLD                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |

``` r
# Map the unique Phi-base IDs of effector data to the Phi-base fatsa data
# `phibase_id_uniq`: list of all unique Phi-base IDs of plant effector data
effector_exists_from_phibase_fasta <- phi_base_fasta %>% 
  rowwise() %>% 
  dplyr::filter(
    dplyr::intersect(
      stringr::str_split(PHIMolConnID, "__") %>% unlist(), 
      phibase_id_uniq
    ) %>% 
      length() != 0
  )

# Get all of the unique ProteinIDs in the data
effector_proteinIDs_available_in_phibase_fasta <- effector_exists_from_phibase_fasta %>% 
  select(ProteinID) %>% 
  unique()

effector_proteinIDs_available_in_phibase_fasta %>% 
  nrow()
```

    ## [1] 379

``` r
# write_csv(effector_proteinIDs_available_in_phibase_fatsa, "effector_proteinIDs_available_in_phibase_fatsa.csv", col_names = FALSE)
```

In order to make sure that the mapping is correct, then we can make a
dataframe containing of the list of mapped Phi-base IDs in every row of
list Phi-base IDs in fasta data:

``` r
# Non-filtered version ----
phi_base_fasta_with_effector_IDs_mapping <- phi_base_fasta %>%
  rowwise() %>%
  dplyr::mutate(
    # Internediary variable with vector of intersection
    intersectIDs = dplyr::intersect(
      stringr::str_split(PHIMolConnID, "__") %>% unlist(),
      phibase_id_uniq
    ) %>% list(),
    # Amount of elements in intersection
    intersectCount = intersectIDs %>% unlist() %>% length()
  ) %>%
  # Replicate original format for more than one ID
  dplyr::mutate(
    PHIMolConnIDIntersect = ifelse(intersectCount > 0,
                          intersectIDs %>% str_c(collapse = ","),
                          NA)
  ) %>%
  # Reorder variables
  dplyr::select(ProteinID, PHIMolConnID, intersectIDs, PHIMolConnIDIntersect, intersectCount) #%>%
  # Get rid of intermediary variable
  # dplyr::select(-intersectIDs)

phi_base_fasta_with_effector_IDs_mapping %>% 
  head(10) %>% 
  knitr::kable()
```

| ProteinID  | PHIMolConnID         | intersectIDs              | PHIMolConnIDIntersect |  intersectCount|
|:-----------|:---------------------|:--------------------------|:----------------------|---------------:|
| A0A016PMX1 | PHI:6978             | character(0)              | NA                    |               0|
| A0A016PTV5 | PHI:5393             | character(0)              | NA                    |               0|
| A0A016PUA3 | PHI:3659             | character(0)              | NA                    |               0|
| A0A016PW58 | PHI:3658             | character(0)              | NA                    |               0|
| A0A016PX00 | PHI:4602             | character(0)              | NA                    |               0|
| A0A016Q4W5 | PHI:3657             | character(0)              | NA                    |               0|
| A0A023H5D8 | PHI:6442             | character(0)              | NA                    |               0|
| A0A023NA98 | PHI:3354             | character(0)              | NA                    |               0|
| A0A023UJQ9 | PHI:5538\_\_PHI:5557 | c(“PHI:5538”, “PHI:5557”) | PHI:5538,PHI:5557     |               2|
| A0A023Y9U3 | PHI:3462             | character(0)              | NA                    |               0|

### Analysis the sequence data which are available in Phibase FASTA

However, some of the sequences data that are available in Phibase FASTA
have multiple different Pathogen IDs where we need only one unique
Pathogen ID to retrieve the non-effector data from NCBI.

``` r
# Check the data that has more than one PathogenID
effector_exists_from_phibase_fasta_multi_PathoID <- effector_exists_from_phibase_fasta %>% 
  filter(str_detect(PathogenID, "_")) # since different Pathogen IDs are separated by '_'

# Select the ProteinIDs
IDs_effector_exists_from_phibase_fasta_multi_PathoID <- effector_exists_from_phibase_fasta_multi_PathoID %>% 
  select(ProteinID)

# Write the data into CSV data
# write_csv(IDs_effector_exists_from_phibase_fasta_multi_PathoID, "IDs_effector_exists_from_phibase_fasta_multi_PathoID.csv", col_names = FALSE)
  
# Showing the data in a table
effector_exists_from_phibase_fasta_multi_PathoID %>% 
  knitr::kable()
```

| ProteinID | PHIMolConnID                     | Gene            | PathogenID                    | Pathogenspecies                                                                            | MutantPhenotype                                                       | Sequence                                                                                                                                                                                                                                                              |
|:----------|:---------------------------------|:----------------|:------------------------------|:-------------------------------------------------------------------------------------------|:----------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| B3VBK9    | PHI:5495\_\_PHI:5543\_\_PHI:5576 | Ecp6            | 5499\_\_5507                  | Passalora\_fulva\_\_Fusarium\_oxysporum                                                    | effector\_(plant\_avirulence\_determinant)                            | MQSMILFAAALMGAAVNGFVLPRTDDPDCETKATDCGSTSNIKYTVVKGDTLTSIAKKFKSGICNIVSVNKLANPNLIELGATLIIPENCSNPDNKSCVSTPAEPTETCVPGLPGSYTIVSGDTLTNISQDFNITLDSLIAANTQIENPDAIDVGQIITVPVCPSSQCEAVGTYNIVAGDLFVDLAATYHTTIGQIKALNNNVNPSKLKVGQQIILPQDCKNVTTAVA                                  |
| Q2A0P1    | PHI:4848\_\_PHI:5283\_\_PHI:7472 | Avr2            | 5507\_\_27337\_\_40559\_\_317 | Fusarium\_oxysporum\_\_Verticillium\_dahliae\_\_Botrytis\_cinerea\_\_Pseudomonas\_syringae | effector\_(plant\_avirulence\_determinant)                            | MRFLLLIAMSMTWVCSIAGLPVEDADSSVGQLQGRGNPYCVFPGRRTSSTSFTTSFSTEPLGYARMLHRDPPYERAGNSGLNHRIYERSRVGGLRTVIDVAPPDGHQAIANYEIEVRRIPVATPNAAGDCFHTARLSTGSRGPATISWDADASYTYYLTISED                                                                                                   |
| Q888W0    | PHI:582\_\_PHI:7336              | Hopai1\_\_Hopai | 317\_\_318829                 | Pseudomonas\_syringae\_\_Magnaporthe\_oryzae                                               | unaffected\_pathogenicity\_*effector*(plant\_avirulence\_determinant) | MLALKLNTSIAQAPLKKNAEAELRHMNHAEVRAHTPTRFTLNHRAPTYEVAQSALGENHGGWTAVNKFKVTESEVFIHMERSDSRSKGDFAGDKIHLSVAPQHVASAFNAIGKILQADDSPVDKWKVTDMSCASSDLQPEKKRVTQGAQFTLYAKPDRADNTYSPEYMGKMRGMISSIERELHTAGVQQSNNRPASDVAPGHWAYASYRNEHRSERAGSSSQANELEKEPFFQLVSFPDVAASPVKSGASSRSLMPPPWTR |

``` r
# Save into RDS object
# saveRDS(effector_exists_from_phibase_fasta_multi_PathoID, "effector_exists_from_phibase_fasta_multi_PathoID.RDS")
```

Multiple PathogenIDs in one unique ProteinID can be likely to be
incosistency in Phi-base data (based on Martin Urban email), therefore
it is better to check in Uniprot data. Therefore we can get all of the
information of ProteinIDs (with multiple PathogenIDs) above from
Uniprot, so that we can get the correct PathogenID for all of the
ProteinIDs.

``` r
# Store the ProteinIDs above 
IDs_effector_exists_from_phibase_fasta_multi_PathoID_Uniprot_details <- data.table::fread("../../../../data/getting-data-new/binary-class-data/IDs_effector_exists_from_phibase_fasta_multi_PathoID_Uniprot_details.csv")

IDs_effector_exists_from_phibase_fasta_multi_PathoID_Uniprot_details %>% 
  knitr::kable()
```

| Entry  | Entry name    | Status     | Protein names               | Gene names         |  Organism ID| Organism                                                       | Taxonomic lineage (GENUS) |
|:-------|:--------------|:-----------|:----------------------------|:-------------------|------------:|:---------------------------------------------------------------|:--------------------------|
| B3VBK9 | B3VBK9\_PASFU | unreviewed | Extracellular protein 6     |                    |         5499| Passalora fulva (Tomato leaf mold) (Cladosporium fulvum)       | Passalora                 |
| Q2A0P1 | Q2A0P1\_FUSOX | unreviewed | Secreted in xylem 3 protein | SIX3               |        59765| Fusarium oxysporum f. sp. lycopersici                          | Fusarium                  |
| Q888W0 | Q888W0\_PSESM | unreviewed | Type III effector HopAI1    | hopAI1 PSPTO\_0906 |       223283| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | Pseudomonas               |

``` r
effector_exists_from_phibase_fasta_multi_PathoID_new <- IDs_effector_exists_from_phibase_fasta_multi_PathoID_Uniprot_details %>%
  dplyr::rowwise() %>%
  # Parse Organism name into structure used by PHI base
  dplyr::mutate(
    uniprotOrganism = stringr::str_split(Organism, " ") %>%
      unlist() %>%
      .[c(1, 2)] %>%
      stringr::str_c(collapse = "_")
  ) %>%
  # Select important variables
  dplyr::select(Entry, uniprotOrganism) %>%
  dplyr::rename(ProteinID = Entry) %>%
  # Join with IDs_effector_organisms
  dplyr::right_join(effector_exists_from_phibase_fasta_multi_PathoID, by = "ProteinID") %>%
  # Select proper PathogenID
  dplyr::mutate(
    # Index of the element we need to elect
    indexPathogenID = stringr::str_detect(
      Pathogenspecies %>%
        stringr::str_split("__") %>%
        unlist(),
      uniprotOrganism
    ) %>%
      which(TRUE),
    # Match index of the PathogenID
    newPathogenID = PathogenID %>%
      stringr::str_split("__") %>%
      unlist() %>%
      .[indexPathogenID]
  ) %>%
  # For display purposes
  select(ProteinID, newPathogenID, uniprotOrganism, Sequence)
  
effector_exists_from_phibase_fasta_multi_PathoID_new <- effector_exists_from_phibase_fasta_multi_PathoID_new %>% 
  `colnames<-`(c("ProteinID",
                 "PathogenID", 
                 "Pathogenspecies", 
                 "Sequence"))
```

``` r
# Get the column names of the dataframe
effector_exists_from_phibase_fasta %>% 
  colnames()
```

    ## [1] "ProteinID"       "PHIMolConnID"    "Gene"            "PathogenID"     
    ## [5] "Pathogenspecies" "MutantPhenotype" "Sequence"

#### Combine all of the data that exists in Phi-base fasta with unique PathogenID and with multiple Pathogen ID

``` r
# Select the relevant columns of the data exists in Phi-base 
effector_phibase_fasta_with_uniq_pathogenIDs <- effector_exists_from_phibase_fasta %>% 
  filter(!str_detect(PathogenID, "_")) %>% 
  select(ProteinID, PathogenID, Pathogenspecies, Sequence)

effector_phibase_fasta_all <- effector_phibase_fasta_with_uniq_pathogenIDs %>% 
  rbind(effector_exists_from_phibase_fasta_multi_PathoID_new)
  
effector_phibase_fasta_all %>% 
  nrow()
```

    ## [1] 379

``` r
# Check how many Phi-base IDs mapped to phi-base fasta
phi_IDs_exist <- phi_base_fasta_with_effector_IDs_mapping %>%
  dplyr::select(PHIMolConnIDIntersect) %>%
  unlist() %>%
  stringr::str_split(",") %>%
  unlist() %>%
  dplyr::intersect(phibase_id_uniq, .)

# Print the length
phi_IDs_exist %>% 
  length()
```

    ## [1] 493

``` r
# Set diff ----
# List of Phi-base IDs whose the sequences are not available in Phi-base fasta data
phi_IDs_not_exist <- phi_base_fasta_with_effector_IDs_mapping %>%
  dplyr::select(PHIMolConnIDIntersect) %>%
  unlist() %>%
  stringr::str_split(",") %>%
  unlist() %>%
  dplyr::setdiff(phibase_id_uniq, .)

phi_IDs_not_exist %>% 
  length()
```

    ## [1] 145

    There are 493 Phi-base IDs that are sucessfully mapped to Phi-base fasta, and results in 379 sequences (since some of the sequences have same Phi-base IDs). There are 145 (shown above) Phi-base IDs whose the sequences are not provided in Phi-base fasta. Therefore, it raises a question: Where can we get the sequences of those 145 Phi-base IDs?

We need to get the Phi-base IDs which are not in the Phi base data

``` r
# Unique Phibase IDs that are not provided in Phibase fasta (or as shown before)
phibase_IDs_for_non_existed_fasta  <- effector_exists_from_phibase_fasta %>% 
  dplyr::select(PHIMolConnID) %>% 
  unlist() %>% 
  stringr::str_split("__") %>% 
  unlist() %>% 
  dplyr::setdiff(phibase_id_uniq, .)

phibase_IDs_for_non_existed_fasta
```

    ##   [1] "PHI:2577" "PHI:2578" "PHI:2579" "PHI:2591" "PHI:2598" "PHI:2599"
    ##   [7] "PHI:2600" "PHI:2612" "PHI:2744" "PHI:2804" "PHI:2942" "PHI:2974"
    ##  [13] "PHI:3217" "PHI:3496" "PHI:4527" "PHI:4528" "PHI:4529" "PHI:4530"
    ##  [19] "PHI:4531" "PHI:4532" "PHI:4533" "PHI:4534" "PHI:4535" "PHI:4536"
    ##  [25] "PHI:4537" "PHI:4538" "PHI:4539" "PHI:4540" "PHI:4541" "PHI:4542"
    ##  [31] "PHI:4543" "PHI:4544" "PHI:4545" "PHI:4546" "PHI:4547" "PHI:4548"
    ##  [37] "PHI:4549" "PHI:4550" "PHI:3938" "PHI:3939" "PHI:3940" "PHI:3969"
    ##  [43] "PHI:3970" "PHI:4153" "PHI:4682" "PHI:4760" "PHI:4838" "PHI:2435"
    ##  [49] "PHI:2436" "PHI:5032" "PHI:5070" "PHI:5123" "PHI:5124" "PHI:5146"
    ##  [55] "PHI:5147" "PHI:5148" "PHI:5149" "PHI:5150" "PHI:5151" "PHI:5152"
    ##  [61] "PHI:5153" "PHI:5154" "PHI:5155" "PHI:5156" "PHI:5157" "PHI:5161"
    ##  [67] "PHI:5167" "PHI:5168" "PHI:5170" "PHI:5176" "PHI:5180" "PHI:5182"
    ##  [73] "PHI:5355" "PHI:5356" "PHI:5357" "PHI:5358" "PHI:5359" "PHI:5360"
    ##  [79] "PHI:5361" "PHI:5362" "PHI:5363" "PHI:5364" "PHI:5488" "PHI:5489"
    ##  [85] "PHI:5490" "PHI:6082" "PHI:6083" "PHI:6160" "PHI:6161" "PHI:6162"
    ##  [91] "PHI:6163" "PHI:6164" "PHI:6165" "PHI:6166" "PHI:6168" "PHI:6169"
    ##  [97] "PHI:6170" "PHI:6171" "PHI:6172" "PHI:6173" "PHI:6174" "PHI:6175"
    ## [103] "PHI:6389" "PHI:7313" "PHI:7736" "PHI:7989" "PHI:8068" "PHI:8096"
    ## [109] "PHI:8097" "PHI:8098" "PHI:8099" "PHI:8167" "PHI:8168" "PHI:8169"
    ## [115] "PHI:8170" "PHI:8171" "PHI:8172" "PHI:8173" "PHI:8174" "PHI:8175"
    ## [121] "PHI:8206" "PHI:8207" "PHI:8276" "PHI:8277" "PHI:8278" "PHI:8326"
    ## [127] "PHI:8327" "PHI:8328" "PHI:8329" "PHI:8364" "PHI:8380" "PHI:8381"
    ## [133] "PHI:8422" "PHI:8423" "PHI:8424" "PHI:8425" "PHI:8426" "PHI:8427"
    ## [139] "PHI:8428" "PHI:8429" "PHI:8430" "PHI:8431" "PHI:8432" "PHI:8448"
    ## [145] "PHI:8491"

Then now, we need to get the information about the Phi-base IDs that do
not exist in Phibase fasta:

``` r
# Map all of the Phibase IDs (that do not match to Phibase fasta) to the effector data  
effector_list_not_in_phibase_fasta <- phi_effector_plant %>% 
  dplyr::filter(PHIMolConnID %in% phibase_IDs_for_non_existed_fasta)

# Get the count of the ProteinsourceID for the sequence that are not in Phibase fasta
effector_list_not_in_phibase_fasta %>% 
  group_by(ProteinIDsource) %>% 
  summarise(count = n())
```

    ## # A tibble: 2 x 2
    ##   ProteinIDsource count
    ##   <chr>           <int>
    ## 1 ""                277
    ## 2 Uniprot           236

According to the list above, we can get around a half of the data from
`Uniprot`, but unfortunately, we cannot proceed to retrieve the data if
there is no information of ProteinID and ProteinIDsource.

``` r
# Get the info of the effector data that can be retrieved from Uniprot
IDs_from_Uniprot <- effector_list_not_in_phibase_fasta %>% 
  dplyr::filter(ProteinIDsource == "Uniprot") %>%
  dplyr::select(PHIMolConnID, ProteinIDsource, ProteinID, PathogenID)
  
# Get the IDs of effector Unipot
effector_IDs_retrieve_from_uniprot <- IDs_from_Uniprot %>%
  select(PHIMolConnID) %>% 
  unique() 

# Number of ProteinIDs that can be used to retrieve data 
effector_IDs_retrieve_from_uniprot %>% 
  nrow()
```

    ## [1] 38

``` r
# Save the IDs into CSV data
# write_csv(effector_IDs_retrieve_from_uniprot, "effector_IDs_retrieve_from_uniprot.csv",col_names = FALSE)
```

``` r
# Get the info of the effector data that can not be retrieved from Uniprot
IDs_from_somewhere_else <- effector_list_not_in_phibase_fasta %>% 
  dplyr::filter(ProteinIDsource == "") %>%
  dplyr::select(PHIMolConnID, ProteinIDsource, ProteinID,  PathogenID)

IDs_from_somewhere_else %>%  
  select(PHIMolConnID) %>% 
  unique() %>% 
  nrow() %>% 
  knitr::kable()
```

|    x|
|----:|
|  107|

Getting the sequence from Uniprot
---------------------------------

Using the ProteinIDs above to get the seqeunce from Uniprot.

``` r
# Read the file obtained after mapping the IDs to Uniprot
details_from_uniprot <- data.table::fread("../../../../data/getting-data-new/binary-class-data/retrieving_from_uniprot_details.csv")
```

Unfortunately, from 38 ProteinIDs mapped, only 37 were successfully
mapped and some of them are deleted, as shown below.

``` r
# Show which proteinIDs are deleted from Uniprot

details_from_uniprot %>% 
  filter(`Protein names` == "Deleted.") %>% 
  knitr::kable()
```

| yourlist   | Entry      | Entry name        | Status | Protein names | Gene names | Organism |  Length|
|:-----------|:-----------|:------------------|:-------|:--------------|:-----------|:---------|-------:|
| Q4PET8     | Q4PET8     | Q4PET8\_USTMA     |        | Deleted.      |            |          |      NA|
| E2MMQ4     | E2MMQ4     | E2MMQ4\_PSEUB     |        | Deleted.      |            |          |      NA|
| E2MM96     | E2MM96     | E2MM96\_PSEUB     |        | Deleted.      |            |          |      NA|
| Q8XQK6     | Q8XQK6     | Q8XQK6\_RALSO     |        | Deleted.      |            |          |      NA|
| Q8XQK7     | Q8XQK7     | Q8XQK7\_RALSO     |        | Deleted.      |            |          |      NA|
| A0A2P4G5R7 | A0A2P4G5R7 | A0A2P4G5R7\_PSEA0 |        | Deleted.      |            |          |      NA|
| F3DUQ1     | F3DUQ1     | F3DUQ1\_PSEA0     |        | Deleted.      |            |          |      NA|
| A0A0F7A7U6 | A0A0F7A7U6 | A0A0F7A7U6\_PSESY |        | Deleted.      |            |          |      NA|

### Parsed the fasta data from Uniprot

Based on the previous data details, for the data which are not deleted,
we can still get the fasta data from Uniprot, that we can transform to
dataframe using the function below:

``` r
parse_fasta_data_uniprot <- function(file_path) {
  # Read FASTA file
  fasta_data <- seqinr::read.fasta(file_path)
  # Number of entries
  num_data <- fasta_data %>% length()


  # Create empty data frame
  parsed_data <- data.frame(
    protein_id = rep(NA, num_data),
    pathogen = rep(NA, num_data),
    sequence = rep(NA, num_data)
  )

  for (i in 1:num_data) {
    # Read 'Annot' attribute and parse the string between 'OS=' and 'OX='
    pathogen <- fasta_data[[i]] %>%
      attr("Annot") %>%
      sub(".*OS= *(.*?) *OX=.*", "\\1", .)

    protein_id <- fasta_data[[i]] %>%
      attr("name") %>%
      strsplit("\\|") %>%
      unlist() %>%
      .[[2]]

    # Concatenate the vector of the sequence into a single string
    sequence <- fasta_data[[i]] %>%
      as.character() %>%
      toupper() %>%
      paste(collapse = "")

    # Input values into data frame
    parsed_data[i,] <- cbind(protein_id, pathogen, sequence)
  }

  return(parsed_data)
}
```

``` r
# Specify the path of the fasta data
uniprot_path <- "../../../../data/getting-data-new/binary-class-data/uniprot.fasta"

# Parsed data
uniprot_parsed <- parse_fasta_data_uniprot(uniprot_path) 

# Raname the dataframe
uniprot_parsed  <- uniprot_parsed %>% 
  # Rename the variables to be consistent to the sequence data from Phi-base fasta
  rename(ProteinID = "protein_id",
           Pathogenspecies = "pathogen",
           Sequence = "sequence") 

uniprot_parsed  %>% 
  knitr::kable()
```

| ProteinID  | Pathogenspecies                                                                                   | Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|:-----------|:--------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| A0A0D1C5E3 | Ustilago maydis (strain 521 / FGSC 9021)                                                          | MHRPTSLYVTLICLLGTVMSVRAATQRVGDSCSYKQNCQDWAEGVGPDWAKGAITCAVPQDGKPEKCGSGDKKGRDEFSGVCKAVGTLSATEYGCACGVHKADCPETSPFSRIFWPQDWLETAWSQVHH                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Q2QCI9     | Pseudomonas syringae pv. tomato                                                                   | MAGINGAGPSGAYFVGHTDPEPASGGAHGSSSGASSSNSPRLPAPPDAPASQARDRREMLLRARPLSRQTREWVAQGMPPTAEAGVPIRPQESAEAAAPQARAEERHTPEADAAASHVRTEGGRTPQALAGTSPRHTGAVPHANRIVQQLVDAGADLAGINTMIDNAMRRHAIALPSRTVQSILIEHFPHLLAGELISGSELATAFRAALRREVRQQEASAPPRTAARSSVRTPERSTVPPTSTESSSGSNQRTLLGRFAGLMTPNQRRPSSASNASASQRPVDRSPPRVNQVPTGANRVVMRNHGNNEADAALQGLAQQGVDMEDLRAALERHILHRRPIPMDIAYALQGVGIAPSIDTGESLMENPLMNLSVALHRALGPRPARAQAPRPAVPVAPATVSRRPDSARATRLQVIPAREDYENNVAYGVRLLSLNPGAGVRETVAAFVNNRYERQAVVADIRAALNLSKQFNKLRTVSKADAASNKPGFKDLADHPDDATQCLFGEELSLTSSVQQVIGLAGKATDMSESYSREANKDLVFMDMKKLAQFLAGKPEHPMTRETLNAENIAKYAFRIVP                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Q88A08     | Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000)                                    | MTIVSGHIGKHPSLTTVQAGSSASVENQMPDPAQFSDGRWKKLPTQLSSITLARFDQDICTNNHGISQRAMCFGLSLSWINMIHAGKDHVTPYASAERMRFLGSFEGVVHARTVHNFYRTEHKFLMEQASANPGVSSGAMAGTESLLQAAELKGLKLQPVLEDKSNSGLPFLIACKQSGRQVSTDEAALSSLCDAIVENKRGVMVIYSQEIAHALGFSVSSDGKRATLFDPNLGEFHTHSKALADTIENISSADGLPLIGVQVFASKIH                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| G7TEI1     | Xanthomonas oryzae pv. oryzicola (strain BLS256)                                                  | MTASPSKLAQLRELSVVVADTGDYDAIKRLQPVDCTTNPTLVKKALDLPVYADLLERELAWGRAHGGDDRTTTVDEVADRLTIGVGVKLSALVPGRVSTEVDADLAHDTQATIAKARKFIAMYAERGVPKDKILIKIAATWEGIEAARQLQLEGIDCNLTLIFNRAQALACAEANVFLISPFVGRILDYYVAQGQTPASIDEDPGVVFVRTVYNAFKQRGSSTVVMGASFRSTAQIEALAGCDRLTISPDLLEKLDAEHGELPRKLSPGKADNAQITPIDSDSFASGLAADPMATEKLASGIDTFAKDLHALRKTIADKLAG                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| G1FM79     | Xanthomonas oryzae pv. oryzae                                                                     | MDPIRPRAPSPAREVLPGPQPDRVQPTADRGVSAPAGSPLDGLPARRTMSRTRLPSPPAPLPAFSAGSFSDLLRQFDPSLLDTSLFDSMPAVGTPHTEAAPAEGDEVQSALRAADDPPPTVRVAVTAAQVDLRTLGYSQQQEKIKPNVRSTVAQHHEALVGHGFTHAHIVALSRHPAALGTVAVKYQDMIAALPEATHEDIVGVGKQCSGARALEALLTVAGELRGPPLQLDTGQLVKIAKRGGVTAVEAVHASRNALTGAPLNLTPAQVVAIASNSGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPAQVVAIASNSGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNIGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNGGGKQALETVQRLLPMLCQAHGLTPEQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPAQVVAIASNIGGKQALETVRRLLPVLCQAHGLTPAQVVAIANNNGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNIGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQAHGLTLEQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPAQVVAIACNIGGKQALETVRRLLPVLCQAHGLTPAQVVAIANNNGGKQALETVQRLLPVLCQAHGLTPAQVVAIASNGGKQALETVQRLLPVLCQAHGLTPAQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPEQVVAIASNGGGKQALESIVAQLSRPDPALAALTNDHLVALACLGGRPALDAVKKGLPHAPELIRRVNSRIAERTSDRVTDYAQVVRVLEFFQCHSHPAHAFDEAMTQFGMSRNGLLQLFRRVGVTELEACGGTLPPASQRWHRILQASGMKSAKPSCASAQTPDQASLHAFADSPERDLDAPSPMHEGDQTRASSRKRSRSDRAVTGPSAQQAVEVRVPEQRDALHLPLSWSVKRPRTRIGGGLPDPGTPMAADLAASSTLMWEQDVDHFAGAADDFPAFNEEELAWLRELLPQ                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| B2SU53     | Xanthomonas oryzae pv. oryzae (strain PXO99A)                                                     | MDPIRSRTPSPARELLPGPQPDRVQPTADRGGAPPAGGPLDGLPARRTMSRTRLPSPPAPSPAFSAGSFSDLLRQFDPSLLDTSLLDSMPAVGTPHTAAAPAECDEVQSGLRAADDPPPTVRVAVTAARPPRAKPAPRRRAAQPSDASPAAQVDLRTLGYSQQQQEKIKPKVGSTVAQHHEALVGHGFTHAHIVALSRHPAALGTVAVKYQDMIAALPEATHEDIVGVGKQWSGARALEALLTVAGELRGPPLQLDTGQLVKIAKRGGVTAVEAVHASRNALTGAPLNLTPAQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPAQVVAIASHDGGKQALETMQRLLPVLCQAHGLPPDQVVAIASNIGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQTHGLTPAQVVAIASHDGGKQALETVQQLLPVLCQAHGLTPDQVVAIASNIGGKQALATVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTQVQVVAIASNIGGKQALETVQRLLPVLCQAHGLTPAQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTQEQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPAQVVAIASNIGGKQALETVQRLLPVLCQDHGLTLAQVVAIASNIGGKQALETVQRLLPVLCQAHGLTQDQVVAIASNIGGKQALETVQRLLPVLCQDHGLTPDQVVAIASNIGGKQALETVQRLLPVLCQDHGLTLDQVVAIASNGGKQALETVQRLLPVLCQDHGLTPDQVVAIASNSGGKQALETVQRLLPVLCQDHGLTPNQVVAIASNGGKQALESIVAQLSRPDPALAALTNDHLVALACLGGRPAMDAVKKGLPHAPELIRRVNRRIGERTSHRVADYAQVVRVLEFFQCHSHPAYAFDEAMTQFGMSRNGLVQLFRRVGVTELEARGGTLPPASQRWDRILQASGMKRAKPSPTSAQTPDQASLHAFADSLERDLDAPSPMHEGDQTGASSRKRSRSDRAVTGPSAQHSFEVRVPEQRDALHLPLSWRVKRPRTRIGGGLPDPGTPIAADLAASSTVMWEQDAAPFAGAADDFPAFNEEELAWLMELLPQSGSVGGTI                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| A9T728     | Physcomitrella patens subsp. patens                                                               | MAASSTRRCTIPDGLSLRQCEVYQTLIHVAEGREERGQGPRPIVVISDVGKDYDDAAALLVLKEFHRLGHVELRAVVANLMPADKRTRLARAWLDALGLQNVPVGRGTRGKPDEEEELELEYEFSWNDFVMPADVPQRDGQDLLLEAYQHAKAKGEKLYLLCLSSLQDIHKFASAYPDLVAHYTAEVHMQGGNYISSEGKLEPDRSAANNRYNWEAARAWHSFLQKNALPSYTYTKAVAFAAALSSEVFVELEASGHPIGSYLRRVQVEQDLAFYKQACESDPEKRFAPFMDQQWFLVHKTNWQQRPNADLGAGLPIGEEVIPYLSKIVLYDVHGALGVSGEDVIEKLGIFRSRDRSIEVELKSTGKVVRVHHWIADDEQSVAVPDKVHAVLSALLKGSLLDALRYSNL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| G7TBV9     | Xanthomonas oryzae pv. oryzicola (strain BLS256)                                                  | MRVQRSTAQHLVTRDSAAVETPPAAESPQDREGWPLGLNVLQKRPLPKRTHRDESERRPHHGDANARVTEAEAVYDIFVVAKYLADPLHGTADETAELIDDIPPRGEPGHSPESGESAHPPADDPPADHPPADHPPADHPPADHPPAPDGHHLPPPPHAAHPAQPHPAHAPHAPTPEGHHGPHAGAGHPADPHAVGHAAEDGILHALAEEGAPAGAFDLAMSLSISGAMLPLSSLAIYAAYKETREVAEQRTHLRQQERRLRSEQARLQSGLDTTTPAGAADAYGHALSEAIDTNAYLQRRNARDGAIAATSMASAGVIFTKAASELGIQGGLAIASKSANAAGLVGHSAAAAGAASAACIASSFVLAPLASVAATALGGAFLHQSRREKTRVAVDVGRVQRFLQDLEPGELSPGAQRYQHFVSTKLGEYYSFARQFNQCNKGFVVGGITYTASTLTKVGVSAAVLAGATVAAPVGTGLIIGAGLLGAATMGVGSHQFLLAHGKQKRYRRYQTEDMPGVDRALLAVADLLPAPEAESAASRGAPNVADQLQPTAVDVALTPVAEEAAGSDEDDIVAEPDPGLEPKHRQGASAARDAQPAMPEPTASTADRSSARAADKVVPHHGFELRSALYACIDGQEKALEAFLQSSADDLQKLHAVKARSTDQSKDSGQSAPAAPLRRRIDAALRAAANYERAVLTGKPQEALGKAMRSHAKHSPALTETALAQWLQTPGSVKSQIAYMQCCLELQKKYLDAKLSAQADLPQPDVHDDPHSHGETPVATVDQAPERLLAQLHQTRTRDDTQQRAVGLMLQELQMAEQEQASPDARSQARAAARMPGLQQRLIRVLTSNVIVPQEGVAGFARFCMKQARQHTTHVRGTLLATEIQAARIREQAAAAADASTT                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Q3BY51     | Xanthomonas campestris pv. vesicatoria (strain 85-10)                                             | MKAELTRSHSLSSLERIRHPESDMGNDCDLEVGAQPKHSGFCQSEDSRLPTTRPPRRRSTSSGATSPKSTVQGSGHTDNSGLETPLLISLTGLQKPRHMGLVRRESSRLVSADPVVHALLSFAQVDQPFPPQVASTDGVRLELVSRRDPEKALEEFKDAFTVETAQLMPAANSSERTAEQIDADIHIPLLLKAIERGSAAFGPSALIEMADGSQISAKAFLASCAPDVMSNDDVLSAFINQKLKGDEDLQVRLGAQELLYVATKKEFQLGGLAGSIGVSSVIGSVWELGASELLKKAIFGKNFSPSQYALQLAGIDSVPPLIIETMDTMCVLAIIKGMKGEDWSLRDLLPKAMKAGAISSIVSFPNNVLQYTGFKSGLGDLAANMTTTEAAIFGAASGIPPEVKESEETMRAGLFQSIKDGVMACPSEKMGPEEAIEQMARHALDVAPGESTAVKSMGLAAIVGMIPLIAGSKATGLLSEKVLRIFRSTVFNPIEAIALNALALGGRIHVPGLFESDNAKHARVVQTILARASEQMESEAREITAEELHHILAPRSEFLRHVGTAIVKGMNASFEVLPALARKLGYGETPLNKRIPYQDLAVSNTSHQPAPEP                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Q3C000     | Xanthomonas campestris pv. vesicatoria (strain 85-10)                                             | MSDMKVNFSSKIIDSTPSEEEVATQQDSYTKSGLVAPSLDSQALKKAPRKRVIKENIAALHTSSLERVHQKKVLVQNLAQLQRGLAKINGRVELEELIDGFSVKELLIKRNPKIAEEYGEGNPLMIRSLRFSNPQEVTSKLGAEGKTPAKREVDTICNKSTLHDIVMTPASLVKKEVRMNLISEVPRAKDKQKYRGLPSVVYGQSSRRSESDYLTSRNGFGDVHSFKSNNAFNSDYEKICGSLSHAEKLGLIERNLTPFIRHDPDRISTDFVHSIEELAEHQMLLQSRKPASALRHNEYCTKLELWDAKAIAVGESRALAVATLIEFNLEMLSIAQEIDDDGHKSKMVADFIERQLSWLGPQTALDSKSTLERVSAVTIQEREFIANEISRSLRQGVSLCTYDKDEAGSHIREMSLLDFRVEEIIEGISIFISSKLLHVTNAGEA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| G0T341     | Xanthomonas euvesicatoria                                                                         | MKNFMRSLGFGSSRSSRSSSSNWNEQQADNDEQTPASSPSTSPSQTSSAFSGLPERPRKKAIALEESLNSSNNIPYEMRMYAEAALSAAKHGSSEAITKADVENKYYLAHAYNERFPELHLSCHDSAQSFFSEFMTSEKQAWRSIVRLSPSSMHHAAIDVRFKDGKRTMLVIEPALAYGMKDGEIKVMAGYETLGKNVQNCLGENGDMAVIQLGAQKSLFDCVIFSLNMALCAYQKDSVFDNLHDCLRRNVRCFSSGERKSILHKNIEFIEGDKFLPPIFYKHSHSRGVVGEFISNQPEYAHKNVSTGRTNPSEDLSERVENFRVRRGDLSYSMSIEASRLRKIRKTIES                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Q48BD8     | Pseudomonas savastanoi pv. phaseolicola (strain 1448A / Race 6)                                   | MHRPITAGHTTSRLILDQLKQISRTPSESSVQSALSQQASMSSPVLERSKSAPALLTAAQRTMLAQVGACNSHLTSDENMAINELRLHKPRLPKDTWFFTDPNKDPDDVVTYTLGKQLQAEGFVHITDVVATLGDAEVRSQRAEMAKGVFNKLGLHDVHVSRGRDYAMNSLQSKEHAKFLLEGHALRAGPGEIHRDSLQDMSRRLARAPHGVGIVVIAGMSDINALITTCPDMVRERVDDITIMGGVEPLKDADGFVQPDARAYNNATDMDAARSLYRKAQELGIPLRIVTKEAAYKTAVSPSFYEGIAGSGHPVGHYLRDVQKSALKGLWEGIQAGLLPGLDDSWFFRTFMPNAQIEAAQLDKNKESSFEDIWPKVTKLNLYDPLTLLASVPGAAKLLFKPKAIHTEGFGVVEQVGPDDVTHPEKAKLLMSALAKSALVQSTVAPD                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Q48Q39     | Pseudomonas savastanoi pv. phaseolicola (strain 1448A / Race 6)                                   | MVKVTSSGFTANPLPHHADSVSPANSPPQLPEPLPLVDESPSSHKGGMRNRAHASLNSQVLGLQAVPLQRGRHVRIRSHTDGVSVINDWLAKRPSVQSETSLDNNGKLVRYTAMHHEPLAPRNEAFFTSVPGMLMAVLTVHPDIEYGISGEITADAVAARLAEPPIGLLTGVWQSSHDRAYIERGGVVHTANMEERWTSLTLPGADPREPLRMIGLQADGDVYLHNGRQLWRATDTSAESVPTERLPEGAAVRIGAGREVQWLHEGAVHSNGISRPVELLRPEASGLGVEQSPARPVDLLPLPRGNAALILDDKGRIYQADLKGIGPVEAHRLKLPGDFAQGKGWAVTTMGLSRDNTVHLMLQDQNGRRMSLQRAQNEVLFRPAYLLDRPLLLLYTEGLHVPSETAVQSHVQLDGHAQVGHIDGVLHYQPAPDQPWERLRQPGGEPLTGVTALYSSALGFIDRKPVFALLGDARQVVELKLEGRTSWLPTDAELPRHPAGGPLAVIPDTVELRARLIAQFDEPVQALAVHNNRQSVALMESGLLIAADADGKVRRLPLLQRPIALAIGLNDQLLVLHHPRNQRPQLKRLSAKDDWEPVPIILPGIDHPSNLRATRTGQIQLQLGDNWHTLLPAMTSHDNRPLPARVKPEPEADELPSENFLAGTNARVNQQQASRISTPHHDASAVTTLMGTAANNPLTMTSSVQAVVDTTRAQVGALARDVLGGVTKNTMRAMAHKLGVVLPPTTQERRLASFHYEAKQAYTSGKTLFEQLPTLAQVRVASAVGPSDGEKFGLSHQQTQRLLNLREEKLEALLRDLRKIGFHEGVIMGDMGDGESAGDLTSTTSTPTFRLAELWRRQHSQVNKALSSVGLSRSEDILPNLNQSIKALAGGAALHADRMGEREAELLSVLCEVSEKIMRAGVRLPADDGSADSASSHAPHGLRTSGLMAGLVDYDALLTSTDTQALEMAERLQQDARLPGLCKLGLSSWVQLAAFDDVVTTFREQISLPGSARRTQLLKNLGLPPDAAPDEMAARMSDLFLDLFNRSTFFSTQSRGLEIRGSLGGADWKHLNAFSVGITGEALQVLGVERIGDGKDGDAGLVAFFVRHAKASVSATSGVGIDFKPGLGTGGRVLDSRPGRSMNSTWGGSANLGISGTYQHGQGAAVIIAPSTISDFVRLLFDVNHPDSTQILRTGVNGGSIGLDLFETNLNGSVGANVSVSPFSLSQKYGPQKPSADAAASGADNRRSTASGALSVGATAQAGAHWGQMELHLDHAWAEIIGLEFQGRTDFNLEFNSNLNLGGALSSALGDRPQKLINASTGNGNLQLAGIRVASSDVQLPTDSVVDDKRRGPFLSTASYKRTFDTQVAKPITADEWSQMRQRLASAFPDNIGELGALAYSTSPSERIATINQMIDRIQSSKARNVEAGGALDGNALRRQRLDAAREMSNAGNSVWQASSEIDRASVVEMLHQLRQQEQSAVQHHARAIPGARVEFNLFGRESLETVVFHAIGHIGLGSKLNDMAELRRKTPGLDQVMRCFQSLPNINQVRYVFEMRPQARFAINDALLAREQQASARALGLPGSSASELNWRGVLDKIKNTPDLYRLAAIAVHNTDENPVTSRIGLPLLNVSATGATSHQLFEAEIQFRYGLYDGLQGIEVLEAGSRALQSPLRVLQQSGIQALGQRTQAGDVPYGPPSPRKEATMRTAGEAAALTTNDVWRQLDGKIQRMNSAHEREANAIGSFQHAYGTASAHVDRLLLRIPELPAPEIDDRNADGESVRGAFASLKRNHQALDEDVKAMRQASEMVYSIPGERTRQTETSALAHVLSVEKSRRSLGHAMEILAGKGVEAGTATGLELNRVSSQVNALVVRRDELLTQLERGAQEGGSNSEEMAMELQQTTSLLQRLRADLLGERQAMDATAKRLDQANRVALEGERSFSDAVRDMTWGELDKMPQ |
| Q888Y1     | Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000)                                    | MVKVTSSGFTANPLSHHADSVSPANSPPQLPEPVHLVDLSESSRKGGMRNRPHASLNSQVLELQAVPSQRGKHVRVRSHADGESVINAWLAKRPSVQSETSLDNDGKLVRYTPVNHEPLAPRNEAFFTSVPGMLMAVLTVHPEMEHGISGDITADAVAARLAEPPIGLLTGIWQSSHDRAYLERGGVVHTANMEERWAPLTLPGINPREPLRMAGLQADGGVYLHNGSQLWRLTETAAESVTTENLPEGAAVRIGAGGEVHGLHEGALHSNGISRPIELWRPKAGAPGREQSPARPVDLLPLPGGTAALILDDKGRIYHADLKGTGAVEAHRLKLPADFAQGKGWAVTAMGLSRDDTVHLMLQDQNGRRMSLQRAPGEALFRPAYLLDRPLLLLYTEGLHVPSEAAVQSHVQLDGHAQLGHIDGVLHYKAAPDQSWERLKQSGGEPLTGLTALYSSPLGFIDRKPVFALVGDARQVVELKLEGRTSWLPSDAELPRHPAGGPLAVIPDTVALRTSPIAQFDEPVQALAVHGNRRVVALTDSGRLMAADADTPARRLPTLQRPIAIAVGLNDQLLVLHHPHSQRPQLKRLSAKDDWEPVPIILPGIVHPSSLRATRTGQIQVQLGENWHTLLPSMTSHDNQRLPARVKPEPEGDEAPSANFLAGSNALANQQQASRISTPHHDASVVTTLAGTTANNPLTMASSLQAVVDTTRAQVGALARDVVGAAANSTMRAMAHTLGVVLPPTPQEKRLASFHNEAKQAYTSGKILFEHLPSLAQVRVASAVGPSDGERFGLSHQQTQRLLTLREGKLEALLRDLRKIGFHEGVIMGDMGDSDSAHGLVSTTSTPTFRLAELWRRQHSRVDKALSSAGLSRSEDIFPDLNLSINALAGGAALNADRMSEREAELLSVLCEVSEKMMRAGVRLPADDGSVDSAHSQAPYGLRTAGLIAGLVDYDALLSSTDAQALEMAERLQQDARLAALCKLGLSSWGQLAAFDDVVTTFREQISLPGSARRTQLLKNLGLPPDAAPDEMAARMSDLLLDLFNRSTFFSTQSRGLELRGSLGSADWKHLNAFSVGVTGEALQVLGVERIGDGKDGDAGLVAFFVRHAKASVSATSGIGIDFKPGPGTGGRVIDSRPGRSMNSTWGGSTNLGISGAYQHGQGAAVIIAPSTISDFVRLLFDVNHPDTTQILRTGVNGGSIGLDLFETNVNASVGANVSVSPFSLSQKYGPQKPTADAAVSGPDNRRSTASGSLSVGGTAQAGAHWGQMELHLDHAWADIIGLEFQGRTDFNLEFNSGLNLGGALSSALGDNPQKLINASTGNGNLQLAGIRVASSDVQLPTDAVVDDKRRGPFLSTASYKRTFDTEVAKPVTAGEWSQMRQRLAKAFPDNIAELGALDYPTRPGERIATIKQVIDRIQGAKARSVEAVGAMDGKALHRQRFDAAREMSNAGNSVWRASSEIERASIVEMLHQLRQQEQSAVQNHARAIPGARVEFNLFGRESLETVVFHAIGHLGLGSKLNDLAELRRKVPGLDQVMLSFQSLPKVNQVRYVFEMRPQARFAINDALLAREQQASARALGLQGPSGSELNWRGVLDKIKTTPDLYRLAAIAVHNTDENPVTSRIGLPLLNVSATGATSHQLFEAEIQFRYGLYDGLQGVELLEAGNRALQSPLRALQQSGIQALGQRTQAGEVAYGPPSPRKESPLRTAVDAAALTTSDIARQLEVKVQRMNTAHEREANAISSFQQAYGIASAHLDRLLLRIPELPLPEIDDRDVDGGRVRGTFASLQRHHQALDDAISAMHQASEKVYTIPGKQATQEQDPALAQLLSVEKRRRSLGHALETLAGRGVEAGTATGLELNRVSSQVNDLVARRDALLRQRESGVQEGGLDSEELEMELQLTTSVLQRLRADLLGERQAMEATAKRLDQASRAALEGERSFSDAVRDRAWGELDNV   |
| Q888Y7     | Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000)                                    | MHRPITAGHTTSRLILDQSKQISRTPSESSAQSALSQQASMSSPVLERSKSAPALLTAAQRTMLAQVGACNAHLTSDENMAINELRSHKPLLPKDTWFFTDPNKDPDDVVTYTLGKQLQAEGFVHITDVVATLGDAEVRSQRAEMAKGVFNKLELHDVHVSRGRDYAMNSLQSKEHAKFLLEGHALRAGPGEIHRDSLQDMSRRLARAPHGVGIVVIAGMSDINALITTCPDMVRERVDDITIMGGVEPLKDADGFVQPDARAYNNATDMDAARSLYRKAQELGIPLRIVTKEAAYKTAVSPSFYEGIAGSGHPVGHYLRDVQKSALKGLWEGIQAGLLPGLDDSWFFRTFMPNAQIEAAQLDKNKESSFEDIWPKVTKLNLYDPLTLLASVPGAAKLLFKPKAIHTEGFGVVEQVGPDDVTHPEKAKLLMSALAKSALVQSTVAPD                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Q1WEM6     | Erwinia amylovora                                                                                 | MKVSHLTSPAPVVIEHQPRQSEKVSRDGDVIKPWSQLPAAAPSFGGCFGKSKKSRGYDSGSSSGSRSNAGFRLNHVPYVSQQNERMGCWYACTRMLGHSISSGPRLGLPELYDSSGPQGLQQREDVLRLMRNENLAEVSLPESRQFSANELGNLLSRHGPIMFGWQTPAGSWHMSVLTGIDKPNDAIIFHDPQRGPDLTMPLDSFNQRLAWRVPHAMLYSEN                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| P78737     | Pyrenophora tritici-repentis                                                                      | MRSILVLLFSAAAVLAAPTPEADPGYEIVKLFEAANSSELDARGLSLDWTLKPRGLLQERQGSCMSITINPSRPSVNNIGQVDIDSVILGRPGAIGSWELNNFITIGLNRVNADTVRVNIRNTGRTNRLIITQWDNTVTRGDVYELFGDYALIQGRGSFCLNIRSDTGRENWRMQLEN                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Q9C173     | Pyrenophora tritici-repentis                                                                      | MAPIFKTTMLLAVAILPAALVSANCVANILNINEAVIATGCVPAGGELRIFVGSSHSYLIKATSSCGLSLTNQVFINGESVQSGGRC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| Q48M21     | Pseudomonas savastanoi pv. phaseolicola (strain 1448A / Race 6)                                   | MNQVINFLNMIALSAMRRSELVGAFFVIAIVFMMITPLPTGLVDVLIAVNICISCLLIMLAMHLPRPLAFSTFPAVLLLTTMFRLALSISTTRLILLNQDAGHIVEAFGQFVVGGNLAVGLVIFLILTVVNFLVITKGSERVAEVGARFTLDAMPGKQMSIDSDLRANLITVQEARKRRAELNKESQLFGAMDGAMKFVNGDAIASLIIVAINMIGGISIGVLQHNMSAGDALQLYTVLTIGDGLIAQIPALLISVTCGMIITRVPNTEASAEANIGREIAEQITSQPKAWIIAAVAMLGFAALPGMPTGVFITIAIICGAGGMLQLQRARPKPDEQRATAVAPEMNGKEDLRTFSPSRQFVLQFHPGQDSAQIEALISEIRRRRNRLVVQYGLTLPSFIIEHADHLPPDEFRFTVYDVPMLKATFTQSHVAVDARQLGGENLPAAIPGNTDRQEDQWVWLPAEQCGELNPVSAMTLIVERMERALQSCAPQFIGLQETKAILSWLESEQPELAQEIQRVLTLTRFSAVLQRLASERVPLRAIRLIAETLIEHCQHERDTNVLTDYVRIALKSQIYHQYCGAEGLQVWLLTPESEGLLRDGLRQTQTETFFALSNEISQMLVQQMNIAFPVRAPEQAVLLVAQDLRSPLRTLLREEFYHVPVLSFAEISNAAKVKVMGRFDLEDDLEPIDNEHAA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| Q8RP13     | Pseudomonas syringae pv. maculicola                                                               | MGNICIGGPRMSQQVYSPERADTPQRNALREQTEADIEQNSESMRLQQKIDDLKPYVSNATGPMKAYGQAAMDRASGKKTYVSFAELDAEHLDAMVDVENQRNPDLNLRYFKHHKEFIQALESDGPSSFRAIYPLTKPRTGQAAKHHVMADVRLTPCEPPSIVITEPGVIVGKRNKQLHRHNQTLEDLSESGVQLSQVAIIETQAQKTPDDCVMYCLNYAIKAHKNADKFDDIHHGLQRGTFVNRINGRRVQDADHCGHF                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Q4ZX47     | Pseudomonas syringae pv. syringae (strain B728a)                                                  | MNISGPNRRQGTQAENTESASSSSVTNPPLQRGEGRRLRRQDALPTDIRYNANQTATSPQNARAAGRYESGASSSGANDTPQAEGSMPSSSALLQFRLAGGRNHSELENFHTMMLNSPKASRGDAIPEKPEAIPKRLLEKMEPINLAQLALRDKDLHEYAVMVCNQVKKGEGPNSNITQGDIKLLPLFAKAENTRNPGLNLHTFKSHKDCYQAIKEQNRDIQKNKQSLSMRVVYPPFKKMPDHHIALDIQLRYGHRPSIVGFESAPGNIIDAAEREILSALGNVKIKMVGNFLQYSKTDCTMFALNNALKAFKHHEEYTARLHNGEKQVPIPATFLKHAQSKSLVENHPEKDTTVTKDQGGLHMETLLHRNRAYRAQRSAGQHVTSIEGFRMQEIKRAGDFLAANRVRAKP                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Q4ZPS9     | Pseudomonas syringae pv. syringae (strain B728a)                                                  | MGLCISKHSGSSYSYSDSDRWEVPVNPSNVRSASSHQTVSASDRASDKVDERPATFSHFQLARCGEDYTLSMVSLAAYQAERRHRGNLIKDRSQSALPWVQVYHSETGLDYSFQIDRTTTVKVAGFNYNVPNDGETRHLYSAGTSQVNMPVITDNMSACIAVACAAENVDAGTGERRPGAKVRVFHLLPFRREDLMPEQVLASVRDYLRDIKEQGLTMRVALHGGNREGDFSVSTAEALKSLFADEGIPIEFDETCANRTSETLLGAVILNDNSTQFIKHLVAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| A0A0N8SZV2 | Pseudomonas syringae pv. syringae                                                                 | MPINRPAFNLKLNTAIAQPTLKKDAGAELRRLNQSEVRANTQTRFAVNHRAPTYDVAQSALGENHGGWTAANHFKMTGSEVFIHMDRLEPNCKGEFAGDKIHLSVAPEDVPDAFNAIGKILQASDSPVDSWKVTDMKCLQAEMPAAKQRVALGAQFTIYAKPDREDNTYSPEYMGKMRGMISSIEQELSAAGVRQSSHRPDSDVSPGHWSYASYRNEHKSNRSGTSNQHRNLEAEPFFQLVSFSDGASGSSRSSADHQALLPPPWAR                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Q4ZX82     | Pseudomonas syringae pv. syringae (strain B728a)                                                  | MIGTRVGGSGSTEIVQANQPQPSAAVAQAHPHAVSPSSNPPLTASQSAAQAPESSAAGAARLPVAPRHLPTLEKFRAEQPTVQGTSTPTISANAALLIGSLLQSEKLPFEVMAARLSPERYALQQFHGSDLQQMLGRFAEPGHLPGKAETEQLIKGFARSLADQLEHFQLMHDATAEAFGPGGLRDRNTLAVSQAALGEYAGRASKSIEAGLNHSLAVLDERIAALDSQLEGATEDSRPVLLMDRQALETARAMLSDLHVDFCKSPEAKRLSAVAAHTQMDALIDKLNVDRSSVGGWKGIGPIVAAAVPQFMVSMLHLGYIRTATSDAMKDAVPEKSADASMKRALAVGLTAGVAHEGVTNLLKPMVQAGFQKAGLNERLNMVPLKGIDTDSVIPDPFELKNDNGALVRKTPEEAAEDKAFVASERAVLNQKKVQVSSTHPLGEMIPYGAFGGGQAVRQMLNDFNLLNGQTLSARAVTSGIAGAISATTQTIAQLNSTYVDPRGRKIPVFTPDRANADLGKDLAKGLDLREPAVRTAFYSKAVSGVQSAALNGALPSVAVQPQGASGTLSAGNIMRNMALAATGSVSYLSTLYANQSVTAEAKALKEAGMGGATPMVARTETALSNIRHPDRASLPHTFQPDTLGGVPRAVENAYHMARGALQLPTQVVVDTVRVVEDGVASGVSSLRDAHKPAETSSPTADDAAAVELTAMEEGRRR                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Q4ZX48     | Pseudomonas syringae pv. syringae (strain B728a)                                                  | MRIHSSGHGISGPVSSAETVEKAVQSSAQAQNEASHSGPSEHPESRSCQARPNYPYSSVKTRLPPVASAGQSLSETPSSLPGYLLLRRLDRRPLDQDAIKGLIPADEAVGEARRALPFGRGNIDVDAQRSNLESGARTLAARRLRKDAETAGHEPMPENEDMNWHVLVAMSGQVFGAGNCGEHARIASFAYGASAQEKGRAGDENIHLAAQSGEDHVWAETDDSSAGSSPIVMDPWSNGPAVFAEDSRFAKDRRAVERTDSFTLSTAAKAGKITRETAEKALTQATSRLQQRLADQQAQVSPVEGGRYRQENSVLDDAFARRVSDMLNNADPRRALQVEIEASGVAMSLGAQGVKTVVRQAPKVVRQARGVASAKGMSPRAT                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Q4ZV89     | Pseudomonas syringae pv. syringae (strain B728a)                                                  | MITPSRYPGIYIAPLSSEPTAAHTFKEQAEKALDHISAGPSGRELLQEISKLASKKDRKVTLKEIEINNQCYTDAVLSRRQLEKYEPEDFNENRRIASRLSRKGTFTKGEGSNAIIGWSPDKASIRLNQNGSPLRLGMDNDDKITTLAHELVHVRHVLSGSSLADDGDRYNPRTGSGKEELRAVGLDKYSYSLTKEPSENSIRAEHGLPLRMKYRPHQ                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Q8P8Y0     | Xanthomonas campestris pv. campestris (strain ATCC 33913 / DSM 3586 / NCPPB 528 / LMG 568 / P 25) | MSDMKVNFSSKIIDSTPSEEEVATQQDSYTKSGLVAPSLDSQALKKAPRKRVIKENIAALHTSSLERVHQKKVLVQNLAQLQRGLAKINGRVELEELIDGFSVKELLIKRNPKIAEEYGEGNPLMIRSLRFSNPQEVTSKLGAEGKTPAKREVDTICNKSTLHDIVMTPASLVKKEVRMNLISEVPRAKDKQKYRGLPSVVYGQSSRRSESDYLTSRNGFGDVHSFKSNNAFNSDYEKICGSLSHAEKLGLIERNLTPFIRHDPDRISTDFVHSIEELAEHQMLLQSRKPASALRHNEYCTKLELWDAKAIAVGESRALAVATLIEFNLEMLSIAQEIDDDGHKSKMVADFIERQLSWLGPQTALDSKSTLERVSAVTIQEREFIANEISRSLRQGVSLCTYDKDEAGSHIREMSLLDFRVEEIIEGISIFISSKLLHVTNAGEA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| A0A2P1C6A6 | Puccinia striiformis f. sp. tritici                                                               | MIFNYRSIALLLLAAASEVASLKAGECGRPGDKAVCFSEAQYKGKDDLINAKYQEDGKICLNFDTPDFARTRCCKKAVLDKLKPDGNKVIKIPRASLVHGDINAPDQGDGNDCI                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |

In order to get the non-effector data later on, we need `Pathogen ID`,
then now we need to get the *PathogenID* for each PorteinIDs above
(mapping to phi-base data).

### Mapping the Uniprot-parsed data to Phi\_base CSV

Goal: to get column PathogenIDs

``` r
uniprot_parsed_with_pathogen_id <- uniprot_parsed %>% 
  left_join(., phi_base, by = "ProteinID") 

uniprot_parsed_with_pathogen_id <- uniprot_parsed_with_pathogen_id %>% 
  select(ProteinID, 
         PathogenID,
         Pathogenspecies.x, 
         Sequence) %>% 
  group_by(ProteinID, PathogenID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  rename(Pathogenspecies = "Pathogenspecies.x") 
```

We need unique PathogenID for every ProteinID, therefore we need to make
sure as follows

``` r
# Print the ProteinID with more than different PathogenID
uniprot_parsed_with_pathogen_id %>% 
  select(ProteinID, PathogenID, Pathogenspecies) %>% 
  group_by(ProteinID, PathogenID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(ProteinID) %>% 
  filter(n() > 1)
```

    ## # A tibble: 2 x 3
    ## # Groups:   ProteinID [1]
    ##   ProteinID PathogenID Pathogenspecies                                 
    ##   <chr>          <int> <chr>                                           
    ## 1 G7TBV9           317 Xanthomonas oryzae pv. oryzicola (strain BLS256)
    ## 2 G7TBV9           347 Xanthomonas oryzae pv. oryzicola (strain BLS256)

“ProteinID == G7TBV9” has different two different PathogenID, its likely
to be incosistency of Phibase database, now we can just check using
`lib(taxize)`.

``` r
taxize::get_uid(sciname = "Xanthomonas oryzae") %>% 
  knitr::kable()
```

    ## ══  1 queries  ═══════════════

    ## 
    ## Retrieving data for taxon 'Xanthomonas oryzae'

    ## ✔  Found:  Xanthomonas+oryzae
    ## ══  Results  ═════════════════
    ## 
    ## ● Total: 1 
    ## ● Found: 1 
    ## ● Not Found: 0

| ids | class | match | multiple\_matches | pattern\_match | uri                                                                                                           |
|:----|:------|:------|:------------------|:---------------|:--------------------------------------------------------------------------------------------------------------|
| 347 | uid   | found | FALSE             | FALSE          | <a href="https://www.ncbi.nlm.nih.gov/taxonomy/347" class="uri">https://www.ncbi.nlm.nih.gov/taxonomy/347</a> |

According to the checking above, then the correct IDs are 347. Then now,
we remove the data with incorrect PathogenID.

``` r
uniprot_parsed_with_pathogen_id %>% 
  nrow()
```

    ## [1] 29

``` r
# Remove the ProteinID with wrong PathogenID
uniprot_parsed_with_pathogen_id <- uniprot_parsed_with_pathogen_id %>% 
  filter(ProteinID != "G7TBV9" | PathogenID != "317") 
```

Combine all of the effector data collected from Uniprot dan Phi-base
--------------------------------------------------------------------

We can then combine both data frame using rbind()

``` r
effector_phibase_fasta_all %>% 
  class()
```

    ## [1] "rowwise_df" "tbl_df"     "tbl"        "data.frame"

``` r
uniprot_parsed_with_pathogen_id %>% 
  class()
```

    ## [1] "tbl_df"     "tbl"        "data.frame"

``` r
effector_phibase_fasta_all %>% 
  nrow()
```

    ## [1] 379

``` r
uniprot_parsed_with_pathogen_id %>% 
  nrow()
```

    ## [1] 28

``` r
effector_data <- effector_phibase_fasta_all %>% 
  rbind(., uniprot_parsed_with_pathogen_id) 

effector_data %>% 
  head(10) %>% 
  knitr::kable()
```

| ProteinID  | PathogenID | Pathogenspecies         | Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|:-----------|:-----------|:------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| A0A023UJQ9 | 5499       | Passalora\_fulva        | MKSPIVITILATALGALGSYDALPINCRDTTNYCFNGNGRHEVCSYCNQAKEEPLKLGRRGGQRDCGVAGSQCNDVDHQQCDARCCSKIGSPTFYGVRCPYPY                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| A0A059UDR8 | 34373      | Blumeria\_graminis      | MQSVLLLTVLTQSFIATASPLVERSTPLSFAEKRPQKVSYDWTTPYVKDFTIHGSCNATQTNVLRRGLEDAVTLAQHAKEHILVHGKESPIYQKYYGALPTGAVIGWFDTIATANRAGVTFRCDDPDKKCATENRWAGHWRGKDAPSETVICDVSFFERLPLEDLCSRGYKIATGKVYSYWGADLIHRMFHVDIVGQNAITHASHGYQDALNLAAGQNYTQTATNTDSLIYFAVEAYAFDISVPGEGCAGQASSIPDTKVTPTDPNALPDSIVIPSPKEPEKKDEPEKKAEPEKKAEPEALGKNCHTHDDGEVHCV                                                                                                                                                                                                                                                                                                                                                                                                          |
| A0A0A0S3X0 | 5022       | Leptosphaeria\_maculans | MRLANFLFYLAPMIVSSLAFDFVPLSGELDFSQEMVFINLTQQQFSELHLQHQQWHQKNILKRYTLTELDEICQQYNANFRFNSGFCSGKDRRWDCYDLNFPTTQSERRVQRRRVCRGEHQTCETIDVINAFGAHARFPQCVHRFELPINDPIPYKDSYQGQYTVEKALDDSWEDILANTGGSHVDFSYQSGTQHYQGYGLTFACIHCIGGSILRMIHANDPARATVTIKFH                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| A0A0A2ILW0 | 27334      | Penicillium\_expansum   | MFFPSLILAAGSLSTLIQAIPHGAKHHHSLHRRAAATYAVMGGDGEASDGWPTISQWSEYETLWGLNQILIAASCDNSDDETSDINTSIKSIASETGVDARFILAIIMQESKGCVRVQSTNNGVENTGLMQSHDGEGSCNKDGSKTTPCPSSMITQMIQDGTAGTTQGDGLKQCYEAQTGGTAAKYYKAARTYNSGSIASSGNLGQGGATHCYASDIANRVRGWAGDVSECVEATIGTITSGVESALGGDDGSSSTSTSTTAAQSTETAEPVQTSSSAAEQPVTTEPIQTSSAPAQAAETSSAASSATSTETTSVAPAPTWTPSSNVQVAAQTTTPTPSWTTKSAPAATTAPAASSSASGTAPLYPYASSSCQKYYTVKAGDFCDKVTEAVGISFLDLRSLNPGLDEKCSDLWLGYQYCIKA                                                                                                                                                                                                                                                                                                |
| A0A0A2JB06 | 27334      | Penicillium\_expansum   | MDKMLFSFLRLCFVLLTVWTTWVAADFALYENYDEDALIAGLALSSTCLAALNTTVTCNETAVGLLGHGADIHFWTTTDVNNLCTTDCVSSLSSWKDNVATVCAEETTIQGNVVVKARALPLTFTYNSDLVCMQDSSSNWCFLDSQTWQGSDYIRWDPTMCFTNGDDNSTVALQCADPDFDIGDISDDMAALTNLYDKELFCNECFLNLYRQRLLDPWLPVTNFTDYLIDQFDLVQANCSTTLPYSTSASTLYVGTATATTTTSATTTTDSTTTATCLGQMVQPLQNWLTCNDLCDTYNISTGDARVITGDYACYFNQATCFPLPCEIDTVWDTPSCDELATRYSNSTYSVTTAQFLSWNTNVQGSCSGIAAGQRARREKSAPGGTFPKPNATITAPGATGEPTYYTAATAAYPTQTGTISECGDYYLVVAGDDCATVDLRFGLNYTQLQEYNTYLDATCSNLWLNYDICVAPVTAQTVSTDVVILIMAARLQMESVARILPETRPALVLSLELVVPYLVIVEAPVTTVLGQTAIAVLPYQS                                                                                                                                                                        |
| A0A0A2JVC2 | 27334      | Penicillium\_expansum   | MMAPKSLQTGLLILLLAKLKLAWGMTLVPRQVTCDYEAAASSGDTCTSFAAEWGLTEETFASLNPSAACPSLVAGQNYCMVGTISAASTTSSSSSTTSSSTTSSSTTSSSTTSSSTTTSSFTTTTASETTSTAANGVTTPMPVQSGMIDSCSKFDLIKSGDTCASIASTYNIDLSSFYSWNPAVGSSCAYLDVGYYVCVDAVPSPVQSGITESCSKFDLVESGDTCASIVSTYNIPLSSFYAWNPAVGSSCAYLDLGYYVCVGTDSSTVATSTTSAGNGITTPTPDEPGMVSDCTKFWLVGSDDTCTSIASSEDITVADIEKWNPKVGSTCTDVWLGDYICVGV                                                                                                                                                                                                                                                                                                                                                                              |
| A0A0A2JY92 | 27334      | Penicillium\_expansum   | MGLFRIDFLALWLLLLRLATAAAPVPTKDGFCFKYIIQGYDTCALIAKAHGITEADIESFNKNTWAWLGCGRLYQGDFICLSTGKPPMPMALPQATCGPQVPGTTRPANWANLAGMNPCAESKCCAFWGQCGVTDIYCQNCRAPPAGPTATVKGATEAPKSQTSSNTKAGGANQSVNTAAGSKPTSNPNSKSNTTTVKATKITTITKSSTTTKAKAAPEPSVSQSPWTAPWEITLYSKMGCEGDYYHLEGYNKEFLDNKGCLNLHGGLNSKFTETGVTCKWWTDDGFTWSSCDSSKLMKPQSWILKNGYCSAFMFD                                                                                                                                                                                                                                                                                                                                                                                                          |
| A0A0A7DLN4 | 318829     | Magnaporthe\_oryzae     | MQFSTIQLFALMAVGAMANTVAPADVAAVDAPRAPLLVRRNCEGKNTQAECERFKWLTGCTWLRGSKYCSST                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| A0A0A7DM22 | 318829     | Magnaporthe\_oryzae     | MKFIYSTALISVLAVAIPTGRGLTSGLMRRESGEVSTNYIFTDDPAEGGKTTQKKKRESGEVSTNYIFTDDPAEGGKATQ                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| A0A0D1E5A8 | 5270       | Ustilago\_maydis        | MRFLLAALLLILLSATCQAEFLREHGTHITFGHSPLQKRIVSEAAKTLLDTETPERLKRWAIFQDGNRFREAIRSGYIKLLRRVDKIREPYAFSKGVKKQWIYEPDLRKHYFVPRVLRVPAILQDSPIQHDEAFAAFIKGTNQDFGLQDVVVTHVPELSKFKEAETGKESSLSATKSESQEPLTEQNPPGESMRRRRRTRGYHRRQSFNQFRSATNDAASTEKSPAAEIESEGSLAHQQAGPTQNPEQPAAERVSNAERLSNTLVEMEERFNAMKKFKETQMHKSAMRYVGHPQYRAPVEQESDPAVTAGAPPGLKSIRIRPAVDSQDVTPETTPRQPVDLAGTDADASNPITPHNNWDGTFLYTRFKTPFSDTPRAPSPLPDPSPLVSTGGGMTPRTGDNTGATSPNPFRLPVDYDELLHSARLQRGQNSLHTSHSKPSLPLSSDADSKSDSGTTFGSSDGTASRESAHKLQETSRLNRVDYTVYDHASLDGRRQFASTLFDRLHPYSDEEPVSDALSPSVKRKVVAKEAASFVRSARRPFRSPTNRQFARQTPAKEATHSYSSQPQPVDVEMQDASMSVIQTKASQPDETLLTGATGSKRFKRPPALDKAFVRSQSFLFPTDTRVDPASPFLEGFKPSVKADADSSAELSPMLPQLSPNQLQRLHERLKDMPSTGSSDAVIQAFRSHRIPKSLHSDRFDVRSTFSTS |

``` r
effector_data %>% 
  nrow()
```

    ## [1] 407

``` r
protein_IDs_unique <- effector_data %>% 
  dplyr::select(ProteinID) %>% 
  unique()
  
protein_IDs_unique %>% 
  nrow()
```

    ## [1] 402

There are 402 ProteinID, however if we pay attention of the number of
unique ProteinIDs from separate data, the total is 407 ProteinID,
therefore there must a duplicate from both data.

``` r
effector_data %>% 
  group_by(ProteinID) %>% 
  arrange(ProteinID) %>% 
  filter(n() > 1)
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

    ## # A tibble: 10 x 4
    ## # Groups:   ProteinID [5]
    ##    ProteinID PathogenID Pathogenspecies              Sequence                   
    ##    <chr>     <fct>      <fct>                        <fct>                      
    ##  1 B2SU53    347        Xanthomonas_oryzae           MDPIRSRTPSPARELLPGPQPDRVQP…
    ##  2 B2SU53    347        Xanthomonas oryzae pv. oryz… MDPIRSRTPSPARELLPGPQPDRVQP…
    ##  3 G1FM79    347        Xanthomonas_oryzae           MDPIRPRAPSPAREVLPGPQPDRVQP…
    ##  4 G1FM79    347        Xanthomonas oryzae pv. oryz… MDPIRPRAPSPAREVLPGPQPDRVQP…
    ##  5 Q4ZX47    317        Pseudomonas_syringae         MNISGPNRRQGTQAENTESASSSSVT…
    ##  6 Q4ZX47    317        Pseudomonas syringae pv. sy… MNISGPNRRQGTQAENTESASSSSVT…
    ##  7 Q4ZX82    317        Pseudomonas_syringae         MIGTRVGGSGSTEIVQANQPQPSAAV…
    ##  8 Q4ZX82    317        Pseudomonas syringae pv. sy… MIGTRVGGSGSTEIVQANQPQPSAAV…
    ##  9 Q888Y7    317        Pseudomonas_syringae         MHRPITAGHTTSRLILDQSKQISRTP…
    ## 10 Q888Y7    317        Pseudomonas syringae pv. to… MHRPITAGHTTSRLILDQSKQISRTP…

Now if there is duplicate of ProteinID, we just need to select one of
them:

``` r
effector_data <- effector_data %>% 
  group_by(ProteinID) %>% 
  slice(1)
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

``` r
# Check the number of unique ProteinIDs
effector_data %>% 
  group_by(ProteinID) %>% 
  unique() %>% 
  nrow()
```

    ## [1] 402
