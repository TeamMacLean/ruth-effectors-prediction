Getting secreted data
=====================

Introduction
------------

Ideally, the heatmap plot made in jupyter lab book can show which
positions of the sequences activated the most by convolutional layers,
however after grouping the data based on the pathogen organisms (it
still has not shows similar pattern for each organism). Hypothesis: This
may due to the effector data which is not secreted by the organism.

Aim
---

In this report, several steps will be conducted in order to get the
protein data secreted from each pathogen organisms (to improve the
perfomance of the model in recognising the pattern).

Plan
----

In order to do that, several plans will be done:

1.  Using `rentrez` library, get the genomic sequence ID by inputting
    the GeneID known from data we have (Phi-base data).

2.  If the the genomic sequence ID exist, then:

-   After getting the genomic sequence ID, then get all of protein
    sequences encoded from there.

However, in many cases, the protein does not have genomic sequence ID:

-   If it is the case, then I may need to pick the protein randomly from
    the same pathogen organism and species.

1.  After getting the sequence, we need to use SignalP to identify the
    signal peptide, this way, we know whether a protein is secreted or
    not. The threshold of probability to be accepted to have signal
    peptide = 90%.

Execution
---------

### Load of the data

``` r
# Loading the data

# Load the effector data obtained previously
effector_data <- data.table::fread("../../../data/getting-data-new/binary-class-data/effector_data.csv")

# Load phi base data
all_phi_base <- data.table::fread("../../../data/phi-base-current-data/phi-base_without_column_desc.csv")
```

Get all of the ProteinID from the effector data.

``` r
list_protein_id <- effector_data %>%
  pull(ProteinID)
```

``` r
# Getting the information about the Gene names and Gene IDs, and the protein source
phi_base_prot <- all_phi_base %>%
  dplyr::filter(ProteinID %in% list_protein_id) %>%
  group_by(ProteinID, ProteinIDsource) %>%
  slice(1) %>%
  dplyr::select(ProteinIDsource, ProteinID, GeneIDsource, GeneID, Gene, PathogenID, Pathogenspecies, PathogenstrainID, Pathogenstrain)

# Getting the dictinct value of each data
phi_base_prot %>%
  dplyr::distinct(Gene) %>%
  nrow()
```

    ## [1] 402

From the result, it seems that the they have unique name. NOw, using the
unique name we have above, we want to get the genome name.

The example of the genome I can get from gene is `KGO57615` ==&gt; input
in Gene, and then I can get the info of the Genome, and using that
Genome ID, retrive it in Genome NCBI. Get the whole protein sequences.
In this case, using the gene `LysM1` will not result in a specific one.

### Making sure we connect to internet

``` r
# Function to check whether we connect to internet
connected.to.internet <- function() {
  if (curl::curl_fetch_memory("www.google.com")$status_code == 200) {
    return(TRUE)
  } else {
    stop(
      "It seems that you are not connected to the internet. Could you please check?",
      call. = FALSE
    )
  }
}

connected.to.internet()
```

    ## [1] TRUE

Below is a function to get the information about the genomic sequences
by accessing the html file.

``` r
# get_proteins_in_genome <- function(gene){
#   # URL of the queries:
#   #   https://www.ncbi.nlm.nih.gov/genome/?term=<QUERY>
#   # URL to get the Genome ID if the gene ID is known
#   #   https://www.ncbi.nlm.nih.gov/gene/?term=<GENE_ID>
#   # URL to get all protein data in all genome
#   #   https://www.ncbi.nlm.nih.gov/nuccore/<GNOME_ID>
#
#   # gene_name <- paste0(gene)
#   # Get the gene information
#   # Base URL and queries
#   base_url <- "https://www.ncbi.nlm.nih.gov"
#   search_gene_query <- "/gene?term="
#
#   # Download HTML file
#   search_gene_url <- paste0(base_url, search_gene_query, gene)
#   search_gene_html_file <- textreadr::read_html(search_gene_url)
#
#   # Get the genomic sequence ID
#   genomic_sequence_id <- search_gene_html_file %>%
#   stringr::str_detect("Genomic Sequence") %>%
#   which() %>%
#   `+`(1) %>%
#   search_gene_html_file[.]
#
#   # Get the all of sequence from the whole genomic sequence
#   search_whole_genome_query <- "/nuccore/"
#   search_whole_genome_url <- paste0(base_url, search_whole_genome_query, genomic_sequence_id)
#   search_whole_genome_html_file <- textreadr::read_html(search_whole_genome_url)
#
#   print(search_whole_genome_url)
#   # return(search_whole_genome_html_file)
# }
```

### Checking if a protein has Gene and Gene ID or not

``` r
# Function if protein ID has protein ID
check_if_genomic_seq_exist <- function(gene_term) {
  if (nchar(gene_term) > 0) {
    # Get the Gene ID from NCBI
    ncbi_gene_id <- entrez_search(db = "gene", term = gene_term)$ids

    # Return boolean whether the sequence exists or not
    seq_exists <- (length(ncbi_gene_id) != 0)

    return(seq_exists)
  } else {
    return(NA)
  }
}
```

``` r
df_check_genomic_seq <- phi_base_prot %>%
  ungroup() %>%
  rowwise() %>%
  dplyr::mutate("ncbi_genomic_seq_exist" = check_if_genomic_seq_exist(GeneID))
```

``` r
df_check_genomic_seq %>% 
  group_by(ncbi_genomic_seq_exist) %>% 
  summarise(count = n()) %>% 
  knitr::kable()
```

    ## Warning: Grouping rowwise data frame strips rowwise nature

| ncbi\_genomic\_seq\_exist |  count|
|:--------------------------|------:|
| FALSE                     |    201|
| TRUE                      |    185|
| NA                        |     16|

From the table above, we can see that more than a half of the data does
not have the genomic sequence. Therefore, we will first get the protein
sequences for the data that does have genomic sequence.

We can now get the ID of the gene terms with the genomic sequence exists
in NCBI.

``` r
# Get the list of geneID with the genomic sequence exist in NCBI
geneID_with_genomic_seq <- df_check_genomic_seq %>% 
  dplyr::filter(ncbi_genomic_seq_exist == TRUE) %>% 
  pull(GeneID)

geneID_with_genomic_seq
```

    ##   [1] "KGO57615"     "KGO51991"     "KGO58766"     "KGO59781"     "XP_011389576"
    ##   [6] "KIS70866"     "ABO47652"     "ABO48503"     "CAM98538"     "CAN00074"    
    ##  [11] "ABQ81647"     "ABS50087"     "EDP89767"     "XP_001774397" "ABZ10804"    
    ##  [16] "ACD57163"     "ACN59479"     "EEY65678"     "EEY67147"     "EEY67160"    
    ##  [21] "EEY70623"     "EEY55256"     "EEY55258"     "EEY55635"     "EEY56155"    
    ##  [26] "EEY60877"     "EEY61979"     "EEY66592"     "EEY67823"     "EEY54138"    
    ##  [31] "CBA22834"     "CBA22840"     "CBA22873"     "CBA22883"     "CBA22887"    
    ##  [36] "CBA23250"     "CBA23268"     "CBA23273"     "ADE28519"     "CBJ35355"    
    ##  [41] "AEK81007"     "AEK80501"     "AEK80452"     "AEK80609"     "CBX90026"    
    ##  [46] "um05319"      "um05318"      "um05314"      "um05312"      "um05310"     
    ##  [51] "um05308"      "um05303"      "um05302"      "um05299"      "um05294"     
    ##  [56] "EGY17220"     "EGY23343"     "EGY23742"     "EGY15993"     "EGY15834"    
    ##  [61] "um05309"      "um05301"      "AAO54894"     "AAO54892"     "EHA55842"    
    ##  [66] "EHA55132"     "EHA50277"     "EGZ16757"     "EGZ11358"     "EGZ06166"    
    ##  [71] "EAQ71199"     "AEQ94695"     "AEQ98497"     "AEQ98514"     "AEQ98519"    
    ##  [76] "AEQ94827"     "AEQ97580"     "AEQ95160"     "AEQ96548"     "AEQ96593"    
    ##  [81] "AEQ94400"     "AEQ95731"     "AEQ95741"     "AEQ95751"     "AEQ94587"    
    ##  [86] "CBZ39993"     "CAA65843"     "AAA25726"     "CAA34257"     "CAJ26169"    
    ##  [91] "CAJ26046"     "CAJ25517"     "CAJ24869"     "CAJ24824"     "CAJ23833"    
    ##  [96] "CAJ22828"     "CAJ22827"     "CAJ22212"     "CAJ19916"     "CAJ19898"    
    ## [101] "AAZ37967"     "AAZ37987"     "AAZ37975"     "AAZ34495"     "AAZ37024"    
    ## [106] "AAY50743"     "AAY49067"     "AFP74845"     "AAY39946"     "AAY39686"    
    ## [111] "AAY38843"     "AAY36933"     "AAY36274"     "AAY36240"     "AAY36237"    
    ## [116] "AAY35802.1"   "AAO58160"     "DAA00391"     "AAO28780"     "AAO29541"    
    ## [121] "AAO28848"     "AAO58206"     "AAO58163"     "AAO58156"     "AAO58043"    
    ## [126] "AAO58039"     "AAO58038"     "AAO57782"     "NP_793972"    "AAO57459"    
    ## [131] "AAO55088"     "NP_791204"    "AAO54897"     "AAO54440"     "AAO54439"    
    ## [136] "AAO54435"     "AAO54417"     "AAO54131"     "AAO54046"     "AAO53599"    
    ## [141] "AAO59043"     "AAO59098"     "AAM20937"     "AAM42989"     "AAM41388"    
    ## [146] "AAM40151"     "AAM38068"     "AAM37935"     "AAM35178"     "AAM39311"    
    ## [151] "AAM39261"     "AAM39243"     "AAM39226"     "CAD18611"     "CAD18535"    
    ## [156] "CAD18432"     "CAD18390"     "CAD18281"     "CAD18065"     "CAD18005"    
    ## [161] "CAD17998"     "CAD17973"     "CAD17723"     "CAD17455"     "CAD17447"    
    ## [166] "CAD17369"     "CAD17367"     "CAD17366"     "CAD17364"     "CAD17250"    
    ## [171] "CAD16866"     "CAD17000"     "CAD16962"     "CAD16604"     "AL646052"    
    ## [176] "CAD15541"     "CAD15425"     "CAD15059"     "CAD15058"     "CAD14597"    
    ## [181] "CAD14570"     "CAD14528"     "CAD13849"     "AAM48170"     "AAD47203"

``` r
# Get the list of geneID with the genomic sequence does not exist in NCBI
geneID_without_genomic_seq <- df_check_genomic_seq %>% 
  dplyr::filter(ncbi_genomic_seq_exist == FALSE) %>% 
  pull(GeneID)

geneID_without_genomic_seq
```

    ##   [1] "AHY02126"        "KJ571201.1"      "AIW07468"        "AIP86908"       
    ##   [5] "AIP86909"        "AJR19786"        "PSTG_02549"      "ALD51957"       
    ##   [9] "KPY75424"        "KPY88691"        "BAU98376"        "AVJ35632"       
    ##  [13] "AUD40034"        "AUD40033"        "AUD40036"        "AUD40032"       
    ##  [17] "AUD40035"        "AUD40039"        "AUD40038"        "AUD40037"       
    ##  [21] "ABK13738"        "CAJ90695"        "ABV66276"        "KC312950"       
    ##  [25] "ABZ10809"        "AAS46025"        "ACF19427"        "CAQ53119"       
    ##  [29] "ACN87967"        "BAH47286"        "AB520830"        "AB520831"       
    ##  [33] "ACO58459"        "AB498874"        "AB916598"        "AB498876"       
    ##  [37] "ACN69116"        "ACY39289.1"      "BAI67929"        "BAI67930"       
    ##  [41] "ADE28521"        "ADM80412"        "CBI63251"        "AEF57433"       
    ##  [45] "AEF57434"        "AEF57436"        "AEF57438"        "AEF57440"       
    ##  [49] "AEF57442"        "AEF57443"        "AEF57445"        "AEF57446"       
    ##  [53] "AEF57447"        "AEF57448"        "AEF57449"        "AEF57450"       
    ##  [57] "Mycgr3g105487"   "Mycgr3G111221"   "AEK86698"        "AEK86750"       
    ##  [61] "um10556"         "AEK86768"        "um10554"         "um10553"        
    ##  [65] "AEK86831"        "AEK86877"        "AEK86905"        "AEK86941"       
    ##  [69] "AEK86957"        "AEK86668"        "AEL16453"        "CCD44833"       
    ##  [73] "CCC55811"        "CCC55814"        "CCC55815"        "CCC55817"       
    ##  [77] "CCC55822"        "CCC55825"        "CCC55835"        "CCC55839"       
    ##  [81] "CCC55859"        "CCC55864"        "MGG10097"        "KC312957"       
    ##  [85] "RBL98065"        "AEU17941"        "BAL70272"        "AFD62207"       
    ##  [89] "JN616379"        "AFI43935"        "AFI43936"        "AFV08857"       
    ##  [93] "AEJ88232"        "AEJ88241"        "AGE15681"        "JX134488"       
    ##  [97] "JX134493"        "CCU75836"        "CCU74273"        "CCU76597"       
    ## [101] "CCU83233"        "CCU83284"        "CCU82362"        "CCU83226"       
    ## [105] "CCU82697"        "CCU82707"        "CCU82095"        "CCU82600"       
    ## [109] "EME41286"        "AAB53624"        "AAB53625"        "AAU10320"       
    ## [113] "AAB94815"        "AAA25725"        "AAA88428"        "CAA42824"       
    ## [117] "AAB61464"        "CAA69643"        "CAA78401"        "AAA91019"       
    ## [121] "AAA80239"        "AAB49807"        "AAA86496"        "AAA25728"       
    ## [125] "ABH03482"        "ABH03483"        "ABH07404"        "ABC70473"       
    ## [129] "CAJ29326"        "ACY39283"        "ABB96267"        "ABB96264"       
    ## [133] "ABB96261"        "ABA00713"        "AY842883"        "AAX51198"       
    ## [137] "AAF71492"        "AAF71496"        "AAA83419"        "AAC43431"       
    ## [141] "CAA11940"        "AAA25727"        "CAA59280"        "AAB51082"       
    ## [145] "CAA48009"        "AAB00675"        "CAI72345"        "AAX12108"       
    ## [149] "RCSB10216; 2LAI" "AY785302"        "AAW63763"        "AAT28197"       
    ## [153] "AAQ24627"        "AAA71943"        "AAS66949"        "AAS66948"       
    ## [157] "AAR05401"        "AAQ97593"        "AAQ23181"        "CAG28797"       
    ## [161] "CAG27601"        "CAE55870"        "ACY39282"        "CAE55866"       
    ## [165] "AAF71499"        "AF461560_2"      "AAL84260"        "AAO5441"        
    ## [169] "AAN52341"        "AAN31502"        "AAN31500"        "CAD16675"       
    ## [173] "AAL84247"        "WP_054080048"    "AAL71883"        "CAD18733"       
    ## [177] "CAD18363"        "CAD18173"        "CAD17311"        "CAD17179"       
    ## [181] "CAD16066"        "CAD15839"        "CAD15808"        "CAD15503"       
    ## [185] "CAD15502"        "AAK63068"        "AAK19753"        "AF483831_1"     
    ## [189] "AAK00131"        "AF282857_16"     "AF275317_1"      "AAF71498"       
    ## [193] "AAF67149"        "CAB58262"        "AAD53944"        "AGM16438"       
    ## [197] "AGM61352"        "EPX56639"        "AHB63016"        "AHF65896"       
    ## [201] "EXK24251"

``` r
gene_with_genomic_seq_NA <- df_check_genomic_seq %>% 
  dplyr::filter(is.na(ncbi_genomic_seq_exist)) %>% 
  pull(Gene)

gene_with_genomic_seq_NA
```

    ##  [1] "Avr4"       "AvrPiz-t"   "AvrBsT"     "HaRxL2"     "HaRxL24"   
    ##  [6] "HaRxL18"    "HaRxL36"    "HaRxL68"    "HaRxL70"    "HaRxL73"   
    ## [11] "HaRxL106"   "HaRxLL3a"   "HaRxLL108"  "HaRxLL470b" "hrmA"      
    ## [16] "Brg11"

### Retriving genomic sequence ID for the proteins that have Gene ID from phi\_base data

``` r
# Define a function to get the protein ids if the Gene ID is known
get_protein_NCBI_id <- function(gene_id){
  
  # Retrive the genomic seq NCBI ID id the Gene ID is known 
  get_genomic_id <- entrez_link(dbfrom='protein', id = gene_id , db='all')$links$protein_nuccore
  
  # Get all of protein IDs using 
  protein_ids <- entrez_link(dbfrom = "nuccore", id = get_genomic_id, db = "all")$links$nuccore_protein
  
  return(protein_ids)
}
```

``` r
get_protein_NCBI_id("KGO57615")
```

    ##  [1] "700469201" "700469200" "700469199" "700469198" "700469197" "700469196"
    ##  [7] "700469195" "700469194" "700469193" "700469192" "700469191" "700469190"
    ## [13] "700469189" "700469188" "700469187" "700469186" "700469185" "700469184"
    ## [19] "700469183" "700469182" "700469181" "700469180" "700469179" "700469178"
    ## [25] "700469177" "700469176" "700469175" "700469174" "700469173" "700469172"
    ## [31] "700469171" "700469170" "700469169" "700469168" "700469167" "700469166"
    ## [37] "700469165" "700469164" "700469163" "700469162" "700469161"

### Get the protein sequences from genomic sequences

``` r
# Function to predict the secreted protein, with the probability that the signal peptide exist = 90%
get_signalp_pred <- function(fasta_filename, verbose = FALSE) {
  signalp_path <- "/Users/kristian/Documents/Workspace/Software/signalp/bin"
  verbose_string <- tolower(deparse(substitute(verbose)))
  fasta_path <- fasta_filename

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
get_all_prot_seqs <- function(gene_id) {

  # Get the protein IDs
  ncbi_all_protein_id <- get_protein_NCBI_id(gene_id)

  # Get all of the protein sequences
  prot_fasta <- entrez_fetch(db = "protein", id = ncbi_all_protein_id, rettype = "fasta")

  # Save in temp file
  temp <- tempfile()
  write(prot_fasta, temp)
  parsed_fasta <- ape::read.FASTA(temp, type="AA")

  # Get the SignalP prediction from the fasta file we get
  pred <- get_signalp_pred(temp)

  # Get the information about fasta file if it is secreted
  pred_signalP_exist <- pred %>%
    # dplyr::filter(Prediction == "SP(Sec/SPI)") %>%
    dplyr::mutate(Gene = gene_id)

  return(pred_signalP_exist)
}
```

``` r
get_all_prot_seqs("EXK24251")
```

    ## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
    ## This warning is displayed once per session.

    ## # A tibble: 2 x 6
    ##   ID        Prediction  `SP(Sec/SPI)` OTHER   `CS Position`              Gene   
    ##   <chr>     <chr>       <chr>         <chr>   <chr>                      <chr>  
    ## 1 EXK24252… OTHER       0.000764      0.9992… ""                         EXK242…
    ## 2 EXK24251… SP(Sec/SPI) 0.984236      0.0157… CS pos: 22-23. GAA-LP. Pr… EXK242…

``` r
list_ids_without_NA <- c(geneID_with_genomic_seq, geneID_without_genomic_seq)
```

``` r
# Initialize 
out = NULL
# Create for loop to get the protein contained in the same genomic sequence
for (i in 1:3){
  
    get_prot <- get_all_prot_seqs(list_ids_without_NA)
    
    out = rbind(out, get_prot)
}
```

Get information the negative sequence for protein data that do not have neither genomic sequence in NCBI nor the geneID on the list
-----------------------------------------------------------------------------------------------------------------------------------

``` r
eff_info_without_NCBI_GeneID <- phi_base_prot %>% 
  dplyr::filter(Gene %in% gene_with_genomic_seq_NA) 

eff_info_without_NCBI_GeneID
```

    ## # A tibble: 16 x 9
    ## # Groups:   ProteinID, ProteinIDsource [16]
    ##    ProteinIDsource ProteinID GeneIDsource GeneID Gene  PathogenID
    ##    <chr>           <chr>     <chr>        <chr>  <chr>      <int>
    ##  1 Uniprot         A0A1B0RF… ""           ""     Avr4      685502
    ##  2 Uniprot         C6ZEZ6    ""           ""     AvrP…     318829
    ##  3 Uniprot         G0T341    ""           ""     AvrB…        339
    ##  4 Uniprot         G3C9L0    ""           ""     HaRx…     272952
    ##  5 Uniprot         G3C9L1    ""           ""     HaRx…     272952
    ##  6 Uniprot         G3C9N6    ""           ""     HaRx…     272952
    ##  7 Uniprot         G3C9P0    ""           ""     HaRx…     272952
    ##  8 Uniprot         G3C9Q4    ""           ""     HaRx…     272952
    ##  9 Uniprot         G3C9Q5    ""           ""     HaRx…     272952
    ## 10 Uniprot         G3C9Q7    ""           ""     HaRx…     272952
    ## 11 Uniprot         G3C9R4    ""           ""     HaRx…     272952
    ## 12 Uniprot         G3C9S3    ""           ""     HaRx…     272952
    ## 13 Uniprot         G3C9S7    ""           ""     HaRx…     272952
    ## 14 Uniprot         G3C9U1    ""           ""     HaRx…     272952
    ## 15 Uniprot         Q08370    ""           ""     hrmA         317
    ## 16 Uniprot         Q8XYE3    ""           ""     Brg11        305
    ## # … with 3 more variables: Pathogenspecies <chr>, PathogenstrainID <int>,
    ## #   Pathogenstrain <chr>

``` r
eff_info_without_NCBI_GeneID <- eff_info_without_NCBI_GeneID %>% 
  ungroup() %>% 
  select("PathogenID", "PathogenstrainID")

eff_info_without_NCBI_GeneID <- eff_info_without_NCBI_GeneID %>% 
  dplyr::mutate(pathogen_strain_id = ifelse(!is.na(PathogenstrainID), PathogenstrainID, PathogenID))

eff_info_without_NCBI_GeneID_ids <- eff_info_without_NCBI_GeneID %>% 
  pull(pathogen_strain_id)
```

``` r
eff_info_without_NCBI_GeneID_ids
```

    ##  [1] 685502 318829 316273 559515 559515 559515 559515 559515 559515 559515
    ## [11] 559515 559515 559515 559515    321 267608
