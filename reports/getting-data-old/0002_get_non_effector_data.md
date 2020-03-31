Getting the Non-effector data
=============================

Getting the data from NCBI
--------------------------

In order to get data from NCBI, we can do the following simple query by
specifying the name of organism we want to get.

``` sql
"(Name of Organism)"[Organism]
NOT patent US 
NOT virulence[All Fields] 
NOT effector[All Fields]
NOT elicitor[All Fields]
NOT partial[All Fields] 
NOT multispecies[All Fields] 
AND ("2000"[SLEN] : "2500"[SLEN]) 
NOT "Unknown"[Organism]
NOT hypothetical[All Fields]
NOT uncharacterized[All Fields]
NOT unnamed[All Fields]
NOT putative[All Fields]
```

The tables below depicts the list of the number of sequence data for
each organism and the categorization of them.

### Bacteria data

| Bacteria                  | Number |     | Bacteria                    | Number |
|---------------------------|:------:|-----|-----------------------------|--------|
| Acinetobacter baumannii   |    1   |     | Mycobacterium tuberculosis  | 1      |
| Aeromonas hydrophila      |    2   |     | Pantoea stewartii           | 2      |
| Aeromonas salmonicida     |    1   |     | Pseudomonas aeruginosa      | 2      |
| Brucella abortus          |    1   |     | Pseudomonas cichorii        | 1      |
| Burkholderia glumae       |    1   |     | Pseudomonas savastanoi      | 14     |
| Burkholderia pseudomallei |   22   |     | Pseudomonas syringae        | 42     |
| Campylobacter jejuni      |    1   |     | Ralstonia solanacearum      | 46     |
| Citrobacter rodentium     |    1   |     | Salmonella enterica         | 64     |
| Clavibacter michiganensis |    2   |     | Shigella flexneri           | 2      |
| Coxiella burnetii         |    3   |     | Staphylococcus aureus       | 1      |
| Cystobacter fuscus        |    1   |     | Vibrio parahaemolyticus     | 1      |
| Edwardsiella ictaluri     |    8   |     | Xanthomonas axonopodis      | 8      |
| Erwinia amylovora         |    9   |     | Xanthomonas campestris      | 16     |
| Escherichia coli          |    6   |     | Xanthomonas oryzae          | 22     |
| Francisella tularensis    |    1   |     | Xylella fastidiosa          | 3      |
| Helicobacter pylori       |    1   |     | Yersinia enterocolitica     | 1      |
| Legionella pneumophila    |    7   |     | Yersinia pestis             | 1      |
| Listeria monocytogenes    |    2   |     | Yersinia pseudotuberculosis | 3      |

### Fungi data

| Fungus                    | Number |     | Fungus                    | Number |
|---------------------------|:------:|-----|---------------------------|--------|
| Beauveria bassiana        |    2   |     | Melampsora lini           | 5      |
| Blumeria graminis         |    2   |     | Parastagonospora nodorum  | 2      |
| Botrytis cinerea          |    2   |     | Passalora fulva           | 7      |
| Cercospora apii           |    1   |     | Penicillium expansum      | 4      |
| Cercospora beticola       |    1   |     | Pseudocercospora fuligena | 1      |
| Colletotrichum orbiculare |    1   |     | Puccinia striiformis      | 1      |
| Dothistroma septosporum   |    1   |     | Rhynchosporium commune    | 3      |
| Fusarium oxysporum        |    9   |     | Verticillium dahliae      | 6      |
| Leptosphaeria maculans    |    6   |     | Ustilago maydis           | 2      |
| Magnaporthe oryzae        |   25   |     | Zymoseptoria tritici      | 2      |

### Oomycetes and others data

| Oomycetes                      | Number |     | Others                  | Number |
|--------------------------------|:------:|-----|-------------------------|--------|
| Hyaloperonospora arabidopsidis |    1   |     | Globodera rostochiensis | 1      |
| Phytophthora cactorum          |    2   |     | Heterodera glycines     | 1      |
| Phytophthora capsici           |   10   |     | Macrosiphum euphorbiae  | 2      |
| Phytophthora infestans         |   15   |     | Toxoplasma gondii       | 1      |
| Phytophthora parasitica        |   17   |     |                         |        |
| Phytophthora sojae             |   23   |     |                         |        |
| Pythium aphanidermatum         |    5   |     |                         |        |
| Plasmopara halstedii           |   10   |     |                         |        |
| Phytophthora megakarya         |   11   |     |                         |        |

Note here that the category others is used for the organism which are
neither bacteria, fungi, nor oomycetes.

Function to get the unique data and to create dataframe from fasta data
-----------------------------------------------------------------------

``` r
parse_fasta_data_ncbi <- function(file_path) {
  # Read FASTA file
  fasta_data <- seqinr::read.fasta(file_path)
  # Number of entries
  num_data <- fasta_data %>% length()


  # Create empty data frame
  parsed_data <- data.frame(
    protein_id = rep(NA, num_data),
    protein_fun = rep(NA, num_data),
    pathogen = rep(NA, num_data),
    sequence = rep(NA, num_data)
  )

  for (i in 1:num_data) {
    # Read 'Annot' attribute and parse the string between 'OS=' and 'OX='
    pathogen <- fasta_data[[i]] %>%
      attr("Annot") %>%
      sub(".*\\[ *(.*?) *\\].*", "\\1", .)

    protein_id <- fasta_data[[i]] %>%
      attr("name")

    protein_fun <- fasta_data[[i]] %>%
      attr("Annot") %>%
      stringr::str_remove(protein_id) %>%
      sub(".*> *(.*?) *\\[.*", "\\1", .)

    # Concatenate the vector of the sequence into a single string
    sequence <- fasta_data[[i]] %>%
      as.character() %>%
      toupper() %>%
      paste(collapse = "")

    # Input values into data frame
    parsed_data[i,] <- cbind(protein_id, protein_fun, pathogen, sequence)
  }

  return(parsed_data)
}
```

``` r
get_different_sequence_samples <- function(data, samples, seq_chars = 10) {
  data <- data %>%
    mutate(short_sequence = stringr::str_sub(sequence, 1, seq_chars)) %>%
    group_by(short_sequence) %>%
    slice(1) %>%
    ungroup() %>%
    select(-short_sequence) %>% 
    slice(1:samples)
  
  return(data)
}
```

``` r
# Run function ---------------------------------------------

# Paths
ncbi_others_path <- "../../../data/getting-data-old/ncbi-others.fasta"
ncbi_fungus_path <- "../../../data/getting-data-old/ncbi-fungus.fasta"
ncbi_oomycetes_path <- "../../../data/getting-data-old/ncbi-oomycetes.fasta"
ncbi_bacteria_path <- "../../../data/getting-data-old/ncbi-bacteria.fasta"
ncbi_bacteria_salmonella_path <- "../../../data/getting-data-old/ncbi-bacteria-salmonella.fasta"
ncbi_bacteria_xanthomonas_path <- "../../../data/getting-data-old/ncbi-bacteria-xanthomonas.fasta"

# Parsed data
ncbi_others_parsed <- parse_fasta_data_ncbi(ncbi_others_path)
ncbi_fungus_parsed <- parse_fasta_data_ncbi(ncbi_fungus_path)
ncbi_oomycetes_parsed <- parse_fasta_data_ncbi(ncbi_oomycetes_path)
ncbi_bacteria_parsed <- parse_fasta_data_ncbi(ncbi_bacteria_path)
ncbi_bacteria_xanthomonas_parsed <- parse_fasta_data_ncbi(ncbi_bacteria_xanthomonas_path)
ncbi_bacteria_salmonella_parsed <- parse_fasta_data_ncbi(ncbi_bacteria_salmonella_path)

# get the data from Salmonella and Xanthomonas
unique_bacteria_xanthomonas_parsed <- get_different_sequence_samples(ncbi_bacteria_xanthomonas_parsed, 22)
unique_bacteria_salmonella_parsed <- get_different_sequence_samples(ncbi_bacteria_salmonella_parsed, 64)

# merge all of data by rows
ncbi_noneffector_parsed <- ncbi_bacteria_parsed %>% 
  rbind(., unique_bacteria_xanthomonas_parsed) %>% 
  rbind(., unique_bacteria_salmonella_parsed) %>% 
  rbind(., ncbi_fungus_parsed) %>% 
  rbind(., ncbi_oomycetes_parsed) %>% 
  rbind(., ncbi_others_parsed)

# Save data frames into CSV files
write.csv(ncbi_noneffector_parsed, "../../../data/getting-data-old/ncbi_noneffector_parsed.csv", row.names = FALSE)
```

``` r
# add label on the data frame 
noneffector <- ncbi_noneffector_parsed %>% 
  select(sequence) %>% 
  mutate(label = as.factor(0))

write.csv(noneffector, "../../../data/getting-data-old/noneffector.csv", row.names = FALSE)
```

Now, we can get the data of non-effector with the information of the IDs
and also the organism

``` r
library(tidyverse)
noneffector <- data.table::fread("../../../data/getting-data-old/ncbi_noneffector_parsed.csv") 

noneffector_with_IDs_organism <- noneffector %>% 
  mutate(pathogen_short = word(pathogen, 1, 2, sep = " ")) %>% 
  select(sequence, pathogen_short, protein_id)

# write.csv(noneffector_with_IDs_organism, "../../../data/getting-data-old/noneffector_with_IDs_organism.csv", row.names = FALSE)
```

View the effector and noneffector data
--------------------------------------

Afer all of the process done previously, finally we have the noneffector
data, and we can view both effector and noneffector data as follows.

``` r
effector <- data.table::fread("../../../data/getting-data-old/effector.csv", header = TRUE)
noneffector <- data.table::fread("../../../data/getting-data-old/noneffector.csv", header = TRUE)
```

``` r
# view the first ten data of effector
effector %>% 
  head(10) %>%
  mutate(sequence = substr(sequence, 1, 30)) %>%
  knitr::kable()
```

| sequence                       |  label|
|:-------------------------------|------:|
| MKLSLLSVELALLIATTLPLCWAAALPVGL |      1|
| MHYTTLLLSTLLVGTALAQPTNPPAKTPKK |      1|
| MVQFKTIFLSTALAALFSTGSSSPATKNNV |      1|
| MKFLVLPLSLAFLQIGLVFSTPDRCRYTLC |      1|
| MKFNKTIPLYILAFFSTAVIAGGRKWTNKV |      1|
| MKCNNIILPFALVFFSTTVTAGGGWTNKQF |      1|
| MLFNAAAAAVFAPLLVMGNVLPRNAGNSPG |      1|
| MNFRALFAATVAALVGSTSATTCTTSQQTV |      1|
| MLFYSLFFFHTVAISAFTNIGTFSHPVYDY |      1|
| MRDEMWNTATEPIAIIGSGCKFPGGSTTPS |      1|

``` r
# view the first ten data of non-effector
noneffector %>% 
  head(10) %>%
  mutate(sequence = substr(sequence, 1, 30)) %>%
  knitr::kable()
```

| sequence                       |  label|
|:-------------------------------|------:|
| MPVLTGTAEANSTISIFDGTTLLGTTTADA |      0|
| MLPVPPGSIDDRTVLYTYNVLNQKTSMTRV |      0|
| MSALRDAMSHAALGWVKPELDETLRQARNE |      0|
| MNSKLYKLIFCRRLGCLIAVGEFTRSYGRA |      0|
| MISGAPSQDSLLPDNRHAADYQQLRERLIQ |      0|
| MQEAKVDGGYVKFVNLFKLSDLEQALLPYK |      0|
| MWWSKQKSRITEGQRTLAAPMIMSLEPRML |      0|
| MSSQEFQSIIAGLQHRQITPAEARQRLQRL |      0|
| MHHPSSAQPESASLEAVGDVSAARPASSVP |      0|
| MEIAMAAIWAQVLGIQRVGRQDNFFELGGH |      0|

Splitting the data into training, development, and testing data sets
--------------------------------------------------------------------

In order to splitting data sets, there are several steps need to be
done:

1.  Combine the effector and noneffector data
2.  Using scikit learn package to split the data

``` r
library(tidyverse)
all_data <- effector %>% 
  rbind(., noneffector)

all_data %>% 
  head(10) %>% 
  knitr::kable()
```

| sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |  label|
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------:|
| MKLSLLSVELALLIATTLPLCWAAALPVGLGVGLDYCNSSCTRAFDCLGQCGRCDFHKLQCVH                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |      1|
| MHYTTLLLSTLLVGTALAQPTNPPAKTPKKAPKTQPYNPCKPQEVIDTKCMGPKDCLYPNPDSCTTYIQCVPLDEVGNAKPVVKPCPKGLQWNDNVGKKWCDYPNLSTCPVKTPQPKPKKGGVGGKKASVGHPGY                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |      1|
| MVQFKTIFLSTALAALFSTGSSSPATKNNVNQPLDNISRRSEWKSVQISPVKEHSAKTADNTENNHNLEKRVFTSPHMKRTFTLALENTFYAMAWLIDFSFSDDGEPHFSYKLQSFNHEDNPPKILADVVVPLITTSYNSANQYRGKANVELLCHLAKEYVHVYFSVEVFASGASFVIGKIIDYPTVYVNNQFRKVVKFDIAGAI                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |      1|
| MKFLVLPLSLAFLQIGLVFSTPDRCRYTLCCDGALKAVSACLHESESCLVPGDCCRGKSRLTLCSYGEGGNGFQCPTGYRQC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |      1|
| MKFNKTIPLYILAFFSTAVIAGGRKWTNKVIYNDKGPREGSISIRKGAEGDFNCGPGYPGGPDRMVRVHEDNGNIRGMPPGYRLGPDDKEDKGDNQYYSRNGYHVGDGPAEYQNHGGGQWGDGYYGPPGQITNQHGKRQGDQGCHIM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |      1|
| MKCNNIILPFALVFFSTTVTAGGGWTNKQFYNDKGEREGSISIRKGSEGDFNYGPSYPGGPDRMVRVHENNGNIRGMPPGYSLGPDHQEDKSDRQYYNRHGYHVGDGPAEYGNHGGGQWGDGYYGPPGEFTHEHREQREEGCNIM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |      1|
| MLFNAAAAAVFAPLLVMGNVLPRNAGNSPGSNRCDASTFNNGQDFDIPQAPVNDCRQMVENINRDSQFSVSHSWARPFGGYGDCAFNVRVIAGWRNGLVGGADAVDLLTDSVKNFGEANKVSSKGTYNQIVSAEGEVTCDSVDRGGQVRVQWIVASSSYNPSNDD                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |      1|
| MNFRALFAATVAALVGSTSATTCTTSQQTVAYVALVSILSDTSFNQCSTDSGYSMLTATSLPTTEQYKLMCASTACKTMINKIVSLNAPDCELTVPTSGLVLNVYSYANGFSSTCASL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |      1|
| MLFYSLFFFHTVAISAFTNIGTFSHPVYDYNPIPNHIHGDLKRRAYIERYSQCSDSQASEIRAALKSCAELASWGYHAVKNDNRLFRLIFKTDSTDIQNWVQKNFNEIYKECNRDADEISLTCHDKNVYTCVREGVHNLAYALINEKEIVICPPFFNNPVNSREITAGNQDTVILHEMVHIILKEWKDYGYEWDGIHKLDSTESIKNPDSYAIFAQCARYKYC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |      1|
| MRDEMWNTATEPIAIIGSGCKFPGGSTTPSKLWELLKDPKDIVSEIRPDRFDVDKYFHPDHKHHGTSNVRHSYFLEENFKHFDAKFFGIRPQEAMAMDPQQRFLLETVYESLEAAGITISDLKGSQAGVFVGNMGVDYSELLSQDIDAFPTYFAPGTARSILSNRISYFFDLHGPSVTVDTACSSSLVAVHQAVQSLRLGETPVAIVCGANLLLGPAQYIAESKLQMLSPNGRSRMWDASADGYARGEGFASIVLKPLSVALANGDHIECIIRETGCNQDGRTKGITMPSPLAQCKLIQETYKRAGLDLSKSSDRPQYFEAHGTGTPAGDPVEAEAISTAFFGPESGFRRTSHDPKLYVGSVKTVIGHTEGTAGLAGLIKASLAMKAKSIPPNLHLERVNPAVQPFYGNLEIPTRLMDWPEPAPGQPLRASVNSFGFGGANAHVILESYTPAAEVAMVTPTAAAGPVFSPFVFSASSDKALASMLSAYSDYLSLNPTVDLRSVAYTLSQHRSVFDKRAAISAPDLDTLKTKLKARSEEASPSGKTAAVQSLERRPRYLGVFTGQGAQWARMGVDVINASPAARAIFEDLEQSLKTLPEEDRPSWSMLEELLAPPETSRVYQANISQTVCTAVQVMMVQLLRAAGIEFSCVVGHSSGEMAAAYTAGYLSARDAVRAAYFRGVHSQLAKGSNGQPGGMIAVGTNFEDAEELCELDDFKGRLCVAASNSAELVTLSGDLDAVQEVKKILDAEEKFNKQLQVDKGYHSHHMLPCSEPYVASLQKCGIQAQVPGDATACRWISSVYVDDMTNLDCRVQDRYWIENLAKPVMFSQALSHALGGDDKFDSVIEVGPHPALKGPASQTIQACLGERLPYFGCLSRGTDSNEAFAEFLGGVWSTFGSSAVDLAAYERFATGGCDQRLVKGLPSYTWDHDVEHYFQSRLSKVVLHRSTPPNELLGTRLPDDTAGEVRWRNSLHPGELPWLLQHSAQGQTVFPGTGYIATTLEAVKQLFDSSGVQTVEIRDMVIGNALVIEANTGVETLFSLTFINTQTDRITAHFSFCSQQGGSTKLVENASGDLVVLLGEPSEDALPRSFHPGTQMKDIDEERFYEAIDKLGYGYEGPFRALSQLQRRMGAATGLVAIPEKTKHFDQMVLHPAALDAMVQTVLLAYCYPGDTRLQGISLPTGIDCIRFNYGMLSEAARPGCQLPFLSCTAFEGDDVLGGVGGDVGGDVDVFSEDKRFALIQLQGLHTKPLSPPSAATDLQIFSEMEWKTASPEGADMEVRGEKRAYVADLYSSMERVAYFYMRHVDREIGKDRSGLAPHQVRFLEWVDHMCGRVEAGTLPHISRKWDHDTRQDILKIIAKYPDSIDLELMHAVGENLCSVFRGEMNALEPMVKKNMLNRFYSDALGMSPYTEDLARMVGHITHRYPHMNILEVGAGTGGATKVMLRRLQDAFASYTYTDISSGFFADARQVFKAHESKMLFKTLDIEKDIVDQGYEENSFDLVIANLVVHATADLDATMGRLRRLVKPGGHLVLLEITTNDPLRFGFIFGPLPGWWLGGEDGRVHSPCVDVEWWDRVMKRNGFSGADIVTPHHTLGPLSVIMTQAVDHRVQLLRQPTSADFGDFTIDPERLTIVGGVKPLAEGLEQLLKPRYQSVAWIPTLEEVSSHSLPVMGSVLSLVELDEPLFKDMTAQTLEGFKFVFQQSRSVYWITCGASGANPYSNMAAGVARTVALEMRHLRLGFLDFEDAKDATVQRLADRFLEFEILGTLEQQGKLDHLTWYQEPELRFDGKNLLVPRMKLSKDRNGRYNSRRRQLTKNVNPREVPVSLVPTTSGKDFVLKESLSSSSTKHGAQDTVSLRVHYASQRSLRLESSDYLFLVLGTNLSSGEAMFALADSNRSIVHVDRQWTTSYLGNLDHGRHALADLYTQIMASTVVAGLSAGDSLVVLDAETPLSQALSARCAAKGVRLTLLSTTTATSHSEADGTNKTNVRIHPLESRRSIESKLPSNATCFLDLSTNNGSEAAAVINSYIPAQCRVETRDTLTATACQVTRSTSTGGLGPAVGDVLPACWANVEAAGRDLSFFSAAVVTPTELTAAAGNGKTSAPRVGDDALLLITDWTAEAEVGVLVQPADSMVRFRQDKTYWLVGLTGGLALSLCRWMVNRGARYVVMTSRNPKIDKEWLQGVESCGATVKIFSNDVTDRAAVNSAYRTISATLPPIAGVVQGAMVLRDTMFAETTMETIESILGPKVRGSIYLDEIFYSTPLDFFVFLSSVTATSGNPGQSIYAGANMFMNSLAAQRRKRGVAGSSVEIGCIMGNGSVTTILSYEHQKYLFSVGNTWLAEQDFLTMFGEAVLASPPDAPDSVTSVTGLRLQFNDDKPDITWFSNPIFQHLVLQSGNAMQTSLSVARQGTPVKSLLQEAKSSEEVLDILKDAFTAKLVSSLQADPDSNLLEVDLETLGMDSLVAVDLRSWFLAELSVDVPVLKILNGSTARLLLEFVQGLIPASMTPKLDGSDGADAAAQEAPPVAPPVTKPKPDVSVKVPPPHQPVASLKPSGPASPTSPSSATASPGRSRSVASPVTADTPVSPTTSASMASLNDSRKLIRTVPVSFGQSRFWFLGSYNPDPLAFNITSLMRISGPLRTNDFGKAVDKVLNHHEALRTSFVSENDAPVQKIWSSPAFALEQRKIADDESEVVKAYTEVQNTRYNLEAGQTMRIMLLTKSPTKHVLVLGYHHINMDGVSFEVLFSDIEKAYNRTPLDRSVMQFPDFTIREAGEYKSGAWRSELQYWQSKFTSLPEPTPLLSVSKRRTRPVNLSYTTHSVSRRINAEQSQAIHTVGRKFKATPFHFYLSVFKTLIARFSGADDFCIGIADANRKEDKVMGAVGLYLNLLPLRVRSALGQTFGETLADMKKVSQEAFANSKVPFDVLLNELSVPRSSSQTPLFQTFVNYRRGVSEERSFCGCTGAGELISGGQIGYDISLDIVENPGGDALVTLSVQKDLYNVDMANLLLDSYFRLVDSFAKNPATSLNRPAIYDPVAVDKALTLGCGPTLEDSSWPETLIHRIENMSVKYATKFALRNGQNGGLTYSQMIARINDIAAKLIDAKVGTGIVGVMQASTMDFICSILAVWKAGAIYTPLDPRLNSTDRLKAVVDECQPACILVDATTKPLFDSLATNAVQIDVSMVQSSKTLEASPKVAIHAKAPSAAAVFYTSGSTGVPKGITLSHASLTYNIMAATRQFGFKEGVDIMLQQSSFSFDMALAQMLTSLSNGGTLVVVPSHLRGDALGLSQLIVAENVSIVQASPTEYKSLIGVNAQHLKTSKWRVALSGGENMTQSLLEVFRSLGKPDLVLFNGYGPTEATINANTRIVPYHEPNSNPDLPLLTWPNYSISIVDLELNPVPVGVFGEVCIGGAGVGLGYFKNDELTAKAFVADKTAPAEFVAKGWKTKFRTGDLGRLSPDGGLIIEGRIDGDTQVKLRGMRIDLKNIESAILQAGAGKIIDAAVSVRRGGADESEPQYLVGHVVLDADQTPEDSQQDFLAQLIPRLRLPRHMKPSLLVPIRALPQTASHKLDRRALQQLPISDAGQIAKQSQQGAELGSDQARMWKLWKQVIPRDVVSQYSITPQSDFFHVGGTSLLLVNLQSLIAREHGRAPPLHAMFESSTVAAMTDLVLSDDASGSTALIDWEQETSIPTLPPHIIPGGAGNKVSVPPRVVLLTGATGFLGRQLMAFLLRQPSVKRIHCLAVRGGAPPSSAAPFSDPRVSIHAGDLNAPHLGLGEAVAELLFAQADVIIHNGADVSFLKTYATLRATNVGSTRELARLAAPRRIPFHFVSSASITQLTGLDEFGEASMAAWAPPADPRGMSGGYAAAKWASEVLLEKAARAWGLPVVIHRPSSITGEGTNSLDLMGNMFKYIEQLEAVPESDSWKGNFDFVSVENVAADIVQAVVAANVVAAGGVKFIYEAGDIVYPLSMVKDMSEGGAKLPVKTMPLAKWVEKAAEKGLDSMLAEYLIKAASTGTSLAFPRLLKDGN |      1|

``` r
# splitting the sequence and label in different dataframe
all_sequence <- all_data %>%
  select(sequence)

all_label <- all_data %>% 
  select(label)
```

``` r
library(reticulate)
reticulate::use_condaenv("tensorflow", conda = "/anaconda3/bin/conda")
```

``` python
import pandas as pd
# import numpy as np
```

``` python
import os
cwd = os.getcwd()
print(cwd)
```

    ## /Users/kristian/Documents/Workspace/ruth-effectors-prediction/scripts/r-scripts/getting-data-old

``` r
py_all_sequence <- r_to_py(all_sequence)
```

``` python
r.all_sequence
```

    ##                                               sequence
    ## 0    MKLSLLSVELALLIATTLPLCWAAALPVGLGVGLDYCNSSCTRAFD...
    ## 1    MHYTTLLLSTLLVGTALAQPTNPPAKTPKKAPKTQPYNPCKPQEVI...
    ## 2    MVQFKTIFLSTALAALFSTGSSSPATKNNVNQPLDNISRRSEWKSV...
    ## 3    MKFLVLPLSLAFLQIGLVFSTPDRCRYTLCCDGALKAVSACLHESE...
    ## 4    MKFNKTIPLYILAFFSTAVIAGGRKWTNKVIYNDKGPREGSISIRK...
    ## 5    MKCNNIILPFALVFFSTTVTAGGGWTNKQFYNDKGEREGSISIRKG...
    ## 6    MLFNAAAAAVFAPLLVMGNVLPRNAGNSPGSNRCDASTFNNGQDFD...
    ## 7    MNFRALFAATVAALVGSTSATTCTTSQQTVAYVALVSILSDTSFNQ...
    ## 8    MLFYSLFFFHTVAISAFTNIGTFSHPVYDYNPIPNHIHGDLKRRAY...
    ## 9    MRDEMWNTATEPIAIIGSGCKFPGGSTTPSKLWELLKDPKDIVSEI...
    ## 10   MRLVHAVLLPGIIVFVSNGNLLHAHALHEDETGVTAGRQLRAAASE...
    ## 11   MAPYSMVLLGALSILGFGAYAQEAAVREPQIFFNLTYTEYLDKVAA...
    ## 12   MKLFILTFIWLLTASEVIAAAKKLPGCDKDPCKVKEKSGKYKLKIG...
    ## 13   MRLAIMLSATAVAINFATSSAIDQTKVLVYGTPAHYIHDSAGRRLL...
    ## 14   MQFPTPHLLLTTLLATITTADFSRDCPPGSGVGNQAEWSARGVDGT...
    ## 15   MRLSFVLSLVVAIGYVVTCNATEYSDETNIAMVESPDLVRRSLRNG...
    ## 16   MRVCYFVLVPSVALAVIATESSETSGTIVHVFPLRDVADHRNDALI...
    ## 17   MKINLSWKTYAVLSIVSLGGIHAMEHVPAELTRVSEGYTRFYRSPT...
    ## 18   MKINLSWKTYAVLSIVSLGGIHAMEHVPAELTGVSEGYTRFYRSST...
    ## 19   MMMISTKSKLRNIFLLLAIANVNFCVKGHPMNSAKLAEEVKDGVNQ...
    ## 20   MLIMMGIPLRKIVLLSLIFGISIIPMMCEFLEDARDIQGFSRKSGS...
    ## 21   MLFKQCTALKFLIFILGFSIIAAQYVVDPGFGEIECMCGQIARLTQ...
    ## 22   MQDLSFSTIENHLGPAKDRFFGDGFKHVEYSARHVNLTESEANASI...
    ## 23   MKPVQSVSATQHYNPLTNLGQEASCSHEISGITQTEGMQATAATTQ...
    ## 24   MKVFPALTSALVALGTAGVEAEHVQRSLVMGGGTVPVGAKTYTVGL...
    ## 25   MKVTATIAAASMAIAAASADADTTSRQLILGGSIIPSGQKTYSVGI...
    ## 26   MVTLYCVVVGVAGSAFPVDIDENKSVGHLKDAIKEKNASTITCDAK...
    ## 27   MVKLVCAIVGVAGSAFPVDTDASQLVGDLKKAIKAENAMTFTGDAK...
    ## 28   MNLRPALLATLASFAYVSASVINHDQVVPFTQPTPTTALQQAAVKY...
    ## 29   MAALCDISLVIDSQNITKRTLALEKNWVPNESRSCCTVCSRKFHKL...
    ## ..                                                 ...
    ## 934  MLPAKVKLDTNMTKLAERRTAWTFIRALLWKNWLIKRRHPMATACE...
    ## 935  MPGPDDDAVEQTAAVRTNSSGYGLKAAPSFSTSGLGSVASVNSYRS...
    ## 936  MSSTQFSYSSARLRRVKYLQFGVFSPEEIRAMSVTKQTKVNDRIIP...
    ## 937  MGQHFLFHRDLVAATAVHHAVHARLSSDLQSDLVLLGPQWLRVYRV...
    ## 938  MTTICGAFRQLTSACVSAADSPPALKAGAGSSLSSSRELALFSGPS...
    ## 939  MDQPLTVNAPASPPAPAKRSNLRFVPTLLRKNWLLKRKHPVALFFE...
    ## 940  MGLTGAGIIASVVGILGGISLSCGGWSSLSLGARSLFVTTQFLSAF...
    ## 941  MGLTGAGIIASVVGILGGISLSCGGWSSLSLGARSLFVTTQFLSAF...
    ## 942  MGLTGAGIIASVVGILGGVSLSCGGWSSLSLGARSLFVTTQFLSAF...
    ## 943  MAYRQPPQGGARIPEDDAYYPNETPQAGLLMSDQARYNDTIGAATQ...
    ## 944  MGLTAVNIIAAVAEIIGGLSLAFGGWSALSIGARSLFSTTQFLSGF...
    ## 945  MAHQADEYEMYATPSAMEAEISRRGTENNVRMTEGRPVPQQGYVDP...
    ## 946  MFGRPSDKRSLLDEQQQEADYDLHATPAAVEARPSGAYKEFEITRG...
    ## 947  MKVFRYGFIAAFVFATQYTAAINTGTAPGLASGATGGGDTDPVYPT...
    ## 948  MGTYKTQVGGYKATAASKFDVASGNFGVGDLGGNSVSQSSGSGAVV...
    ## 949  MVNVYRALEKVQLPKGLHHERIQWLYKGPKALPHNGYLLYWMQTSV...
    ## 950  MNIELKFLMWRKRRGSLGGSNTPRESPNLMNSNSLTSTPSSMSITS...
    ## 951  MSPYAKISYASKPVSATAYKAGQGATFLKDFDAAAAPVESDGFCER...
    ## 952  MRANSDAPADATATPFNPLVALQTLLSADVMSVKQRIKAAQQLEEY...
    ## 953  MIRITREELDDRVQRARRENASQSHEESRKERHLRKQRLRSAEARY...
    ## 954  MASTQFSYSSARLRRVKYLQFGVFSPEEIRAMSVTKQTKVNDRIIP...
    ## 955  MINDTPAWKALAEHAAEIKSTHLRELLNDDARNAAMRTEQQGIYLD...
    ## 956  MEYAELGEDATSATTAQEQNESAQEVEPEALDYVTSQLVRPLSNWE...
    ## 957  MSAQDEDKNGDGVYVTKSRSLFSMWLHGKAAPSKAHPAVVFRSADV...
    ## 958  MASDAMPVTADAVAIALPTLPPIPPTDQHPIVVPSTSPPSSKVVAS...
    ## 959  MSLPPPRRTVAPAQLAAGAQAHPHAHLLHAQVRDPRGDGGSAGAAG...
    ## 960  MHYAKLFFVILPIFLLHPEQSNGCFSRGNSGNPSDGNPEAQGLKSE...
    ## 961  MEMPSCFFLLFFLMLFVSPSRHQLVTVSNSSSSPIGTTVAFGTPSP...
    ## 962  MTDRCWIPHPTEVWGIAVPAGAPGTCTYKLVEVTDDGTDPRKLEVT...
    ## 963  MNALALLGAPSASRDSARWSSSESSRSWASPEKTEFSVTRACLPLT...
    ## 
    ## [964 rows x 1 columns]
