Getting the data and data analysis
==================================

This report will present how we get and clean the data.

### Getting the Protein IDs

The data are the protein amino acids sequence data. The sequences can be
retrieved from `uni-prot.org` with the protein IDs obtained from
`phi-base.org` by putting keyword *effector (plant avirulence
determinant)*. However, due to the unavailability of then download
option on `phi-base.org` website, then complete data can be obtained
from `phi-base` github repository on
`https://github.com/PHI-base/data/tree/master/releases`. Then we can
import and clean the data on `R` using `tidyverse` library. And from the
data we obtained which is `phi-base_v4-6_2018-12-05.csv`, we can filter
the data only for *effector (plant avirulence determinant)*. The data
are saved into `phi-base_current.csv`.

``` r
library(tidyverse)

phi_base <- data.table::fread("data/phi-base_current.csv", header = TRUE)

# filter all of the data with 'plant avirulence determinant' information
phi_small <- phi_base %>%
  dplyr::filter_all(any_vars(str_detect(., 'plant avirulence determinant')))

# select only the protein ID data
phi_proteinID <- phi_small %>%
  dplyr::select(`Protein ID`)


# find the unique values and remove all of the rows that the IDs are not available
proteinID_unique <- unique(phi_proteinID) %>%
  dplyr::filter_all(any_vars(!str_detect(., 'no data found')))
```

### Retriving the protein sequence

After we made sure that there is no redundancy in our data, we can use
the protein IDs data to retrieve the sequence amino acids sequence data
on `uni-prot.org`, and we will get data in `.fatsa` file format. The
protein IDs that we obtained from
’phi-base.org`is around 496 protein IDs, however only 482 that are succesfully mapped on`uni-prot.org\`.

![Uniprot Screnshoot](/data/images/uniprot-screenshoot.png)

### Reading and cleaning the data

Using R, we can read and clean the data. In order to read the sequence
data in `.fasta` format file on `R`, we can use package `seqinr`
together with `tidyverse`.

``` r
library(seqinr)

# Read FASTA file
fasta_data <- seqinr::read.fasta("../../data/uniprot-data-mapped.fasta")
```

After we read the `.fasta` data, we can see how the data look like, as
follows.

``` r
fasta_data[[1]]
```

    ##  [1] "m" "k" "l" "s" "l" "l" "s" "v" "e" "l" "a" "l" "l" "i" "a" "t" "t"
    ## [18] "l" "p" "l" "c" "w" "a" "a" "a" "l" "p" "v" "g" "l" "g" "v" "g" "l"
    ## [35] "d" "y" "c" "n" "s" "s" "c" "t" "r" "a" "f" "d" "c" "l" "g" "q" "c"
    ## [52] "g" "r" "c" "d" "f" "h" "k" "l" "q" "c" "v" "h"
    ## attr(,"name")
    ## [1] "sp|P22287|AVR9_PASFU"
    ## attr(,"Annot")
    ## [1] ">sp|P22287|AVR9_PASFU Race-specific elicitor A9 OS=Passalora fulva OX=5499 GN=AVR9 PE=1 SV=1"
    ## attr(,"class")
    ## [1] "SeqFastadna"

The data is saved in attribute format, therefore, we can do the
following to get the data that we need (the name of pathogen and their
sequence) in dataframe.

``` r
# Number of entries
num_data <- fasta_data %>% length()


# Create empty data frame
parsed_data <- data.frame(
  pathogen = rep(NA, num_data),
  sequence = rep(NA, num_data)
)

for (i in 1:num_data) {
  # Read 'Annot' attribute and parse the string between 'OS=' and 'OX='
  pathogen <- fasta_data[[i]] %>%
    attr("Annot") %>%
    sub(".*OS= *(.*?) *OX=.*", "\\1", .)

  # Concatenate the vector of the sequence into a single string
  sequence <- fasta_data[[i]] %>%
    as.character() %>%
    toupper() %>%
    paste(collapse = "")

  # Input values into data frame
  parsed_data[i,] <- cbind(pathogen, sequence)
}

# Save data frame into CSV file
write.csv(parsed_data, "data/uniprot-data-mapped.csv", row.names = FALSE)
```

We have data frame with two columns, the first column is the organism
pathogen name and the second one is the sequence data.

``` r
uniprot_data <- data.table::fread("../../data/uniprot-data-mapped.csv")
```

### Analysing the data

As we analyze further the data, the name of pathogen is not unique,
since some of the organism with the same species name are the strains
(those with *pv*).

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0     ✔ purrr   0.2.5
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.8
    ## ✔ tidyr   0.8.2     ✔ stringr 1.3.1
    ## ✔ readr   1.2.1     ✔ forcats 0.3.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::count()  masks seqinr::count()
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
# view the head of the data with restricted number of strings
uniprot_data[40:50,] %>%
  mutate(sequence = substr(sequence, 1, 30)) %>% 
  head(20)  %>%
  knitr::kable()
```

| pathogen                                                        | sequence                       |
|:----------------------------------------------------------------|:-------------------------------|
| Pseudomonas syringae                                            | MGNVCFRPSRSHVSQEFSQSEFSAASPVRT |
| Pseudomonas syringae pv. syringae (strain B728a)                | MGCVSSKASVISSDSFRASYTNSPEASSVH |
| Pseudomonas savastanoi pv. phaseolicola (strain 1448A / Race 6) | MGCITSKPLVSSPQWHNSATNSENLETGQR |
| Pseudomonas savastanoi pv. glycinea                             | MQDLSFSTIENHLGPAKDRFFGDGFKHVEY |
| Pseudomonas savastanoi pv. glycinea                             | MQDLSFSTIENHLGPAKDHFFGDGFKHVEY |
| Pseudomonas syringae                                            | MQSPSIHRNTGSIIQPTVTPDARAATDLQE |
| Pseudomonas syringae                                            | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000)  | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae pv. syringae (strain B728a)                | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae                                            | MGCVSSTSRSTGYYSGYENHEEPRVASSPT |
| Pseudomonas syringae                                            | MTRISTSSVNSSFSYSAPAEEAQNRVSSAP |

Let us take a look further the data we have right now. Here we have many
organism name which start with `Pseudomonas syringae`. Considering the
perfomance of deep learning that will be better if we have more data
training, then we need to be careful in defining the problem.

``` r
library(tidyverse)

# view the head of the data with restricted number of strings
uniprot_data  %>% 
   dplyr::filter(str_detect(pathogen, "Pseudomonas syringae"))  %>%
   head(20)  %>%
   mutate(sequence = substr(sequence, 1, 30)) %>% 
   knitr::kable()
```

| pathogen                                                       | sequence                       |
|:---------------------------------------------------------------|:-------------------------------|
| Pseudomonas syringae                                           | MGNVCFRPSRSHVSQEFSQSEFSAASPVRT |
| Pseudomonas syringae pv. syringae (strain B728a)               | MGCVSSKASVISSDSFRASYTNSPEASSVH |
| Pseudomonas syringae                                           | MQSPSIHRNTGSIIQPTVTPDARAATDLQE |
| Pseudomonas syringae                                           | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae pv. syringae (strain B728a)               | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae                                           | MGCVSSTSRSTGYYSGYENHEEPRVASSPT |
| Pseudomonas syringae                                           | MTRISTSSVNSSFSYSAPAEEAQNRVSSAP |
| Pseudomonas syringae pv. tomato                                | MKIAPVAINHSPLSREVPSHAAPTQAKQTN |
| Pseudomonas syringae pv. syringae                              | MNPIHARFSSVEALRHSNVDIQAIKSEGQL |
| Pseudomonas syringae pv. maculicola                            | MINLSNLSAVATSVARAALGESNTSKINHA |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MISSRIGGAGGVELSRVNQQHDTVPAQTAH |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MYIQQSGAQSGVAAKTQHDKPSSLSGLAPG |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MHRPITAGHTTSRLILDQSKQISRTPSESS |
| Pseudomonas syringae pv. maculicola                            | MKIHNAGLTPPLPGISNGNVGKAAQSSITQ |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MHINRRVQQPPVTATDSFRTASDASLASSS |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MHINQSAQQPPGVAMESFRTASDASLASSS |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MAGINRAGPSGAYFVGHTDPEPVSGQAHGS |
| Pseudomonas syringae                                           | MHANPLSSFNRAQHGNLTNVEASQVKSAGT |
| Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000) | MNPLQPIQHSITNSQMSGGQQLEAEGSQAH |

One alternative we can do is defining the problem as predicting the
organism of effector proteins (inclusing the strains). By defining as
preceding, then now we can neglect the pathovar for each species by
taking the first two words of each pathogen species, by doing the
following.

``` r
num <- nrow(uniprot_data)
num

# Create empty data frame
new_parsed_data <- data.frame(
  new_pathogen = rep(NA, num),
  new_sequence = rep(NA, num)
)

for (i in 1:num){
  new_pathogen <- word(parsed_data[['pathogen']][i],1,2, sep=" ")
  
  new_sequence <- uniprot_data[['sequence']][i]
  
  new_parsed_data[i,] <- cbind(new_pathogen, new_sequence)
  
}
```

Now, we can load the data, and take an analysis further

``` r
new_parsed_data <- data.table::fread("../../data/new_parsed_data.csv", header = TRUE) %>% 
  dplyr::select(-V1)
```

``` r
# view the head of the data with restricted number of strings
new_parsed_data[40:50,] %>% 
  head(10) %>% 
  mutate(new_sequence = substr(new_sequence, 1, 30)) %>% 
  knitr::kable()
```

| new\_pathogen          | new\_sequence                  |
|:-----------------------|:-------------------------------|
| Pseudomonas syringae   | MGNVCFRPSRSHVSQEFSQSEFSAASPVRT |
| Pseudomonas syringae   | MGCVSSKASVISSDSFRASYTNSPEASSVH |
| Pseudomonas savastanoi | MGCITSKPLVSSPQWHNSATNSENLETGQR |
| Pseudomonas savastanoi | MQDLSFSTIENHLGPAKDRFFGDGFKHVEY |
| Pseudomonas savastanoi | MQDLSFSTIENHLGPAKDHFFGDGFKHVEY |
| Pseudomonas syringae   | MQSPSIHRNTGSIIQPTVTPDARAATDLQE |
| Pseudomonas syringae   | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae   | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae   | MGNICVGGSRMAHQVNSPDRVSNNSGDEDN |
| Pseudomonas syringae   | MGCVSSTSRSTGYYSGYENHEEPRVASSPT |

By doing previous steps, we have more training data for some pathogen
species. However, some species only have single sequence which is not
enough to be learned by deep learning models.

``` r
library(tidyverse)
count_new_pathogen <- new_parsed_data %>%
  group_by(new_pathogen) %>%
  summarise(count = n())

knitr::kable(count_new_pathogen)
```

| new\_pathogen                  |  count|
|:-------------------------------|------:|
| Acinetobacter baumannii        |      1|
| Aeromonas hydrophila           |      2|
| Aeromonas salmonicida          |      1|
| Beauveria bassiana             |      2|
| Blumeria graminis              |      2|
| Botryotinia fuckeliana         |      1|
| Brucella abortus               |      1|
| Burkholderia glumae            |      1|
| Burkholderia mallei            |      3|
| Burkholderia pseudomallei      |     19|
| Campylobacter jejuni           |      1|
| Cercospora apii                |      1|
| Cercospora beticola            |      1|
| Citrobacter rodentium          |      1|
| Clavibacter michiganensis      |      2|
| Colletotrichum orbiculare      |      1|
| Coxiella burnetii              |      3|
| Cystobacter fuscus             |      1|
| Dothistroma septosporum        |      1|
| Edwardsiella ictaluri          |      8|
| Erwinia amylovora              |      9|
| Escherichia coli               |      6|
| Francisella tularensis         |      1|
| Fusarium oxysporum             |      9|
| Globodera rostochiensis        |      1|
| Helicobacter pylori            |      1|
| Heterodera glycines            |      1|
| Hyaloperonospora arabidopsidis |     36|
| Hyaloperonospora parasitica    |      6|
| Legionella pneumophila         |      7|
| Leptosphaeria maculans         |      6|
| Listeria monocytogenes         |      2|
| Macrosiphum euphorbiae         |      2|
| Magnaporthe grisea             |      6|
| Magnaporthe oryzae             |     19|
| Melampsora lini                |      5|
| Mycobacterium tuberculosis     |      1|
| Pantoea stewartii              |      2|
| Passalora fulva                |      7|
| Penicillium expansum           |      4|
| Phaeosphaeria nodorum          |      2|
| Phaseolus vulgaris             |      1|
| Phytophthora cactorum          |      1|
| Phytophthora capsici           |     10|
| Phytophthora infestans         |     15|
| Phytophthora parasitica        |      2|
| Phytophthora sojae             |     23|
| Pseudocercospora fuligena      |      1|
| Pseudomonas aeruginosa         |      2|
| Pseudomonas cichorii           |      1|
| Pseudomonas savastanoi         |     14|
| Pseudomonas syringae           |     42|
| Puccinia striiformis           |      1|
| Pythium aphanidermatum         |      1|
| Ralstonia solanacearum         |     46|
| Rhynchosporium commune         |      2|
| Rhynchosporium secalis         |      1|
| Salmonella enterica            |      2|
| Salmonella typhi               |      2|
| Salmonella typhimurium         |     60|
| Shigella flexneri              |      2|
| Staphylococcus aureus          |      1|
| Toxoplasma gondii              |      1|
| Ustilago maydis                |      2|
| Verticillium dahliae           |      6|
| Vibrio parahaemolyticus        |      1|
| Xanthomonas axonopodis         |      7|
| Xanthomonas campestris         |     15|
| Xanthomonas euvesicatoria      |      1|
| Xanthomonas manihotis          |      1|
| Xanthomonas oryzae             |     22|
| Xylella fastidiosa             |      3|
| Yersinia enterocolitica        |      1|
| Yersinia pestis                |      1|
| Yersinia pseudotuberculosis    |      3|
| Zymoseptoria tritici           |      2|

### Further discussion

By looking at the data, and did a little analysis, we can do probably
several options:

1.  Continue with the data obtained above, but we can limit/threshold
    only for the species with enough data (for instance, more than 40
    sequences). Then the organism with the number sequences under the
    threshold will not be used as data training.
2.  The second option is to re-define the problem, for instance instead
    of identifying the organism as multi-class classification, we can
    change to binary classification (effector or non-effector) by
    getting non-effector proteins data.
