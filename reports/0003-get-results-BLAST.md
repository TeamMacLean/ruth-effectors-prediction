Getting the identical protein using BLAST
=========================================

Introduction
------------

After getting all of data sets (training, validation / development, and
testing), we need to make sure that there is no identical protein
between them, otherwise it would not help the deep learning model to
learn the data. In order to identify the identical protein sequence
between the data sets, we can use BLAST.

Using BLAST
-----------

We can use terminal to run BLAST as follows

-   Defining database

Now let us BLAST both validation or development datasets against the
entire training datasets. First, we need to tell BLAST that the protein
sequences in training datasets are the database. That’s done by calling
‘makeblastdb’ (creating database using local file):

``` bash
$ makeblastdb -in blast_train.fasta -dbtype prot
```

-   Check the query to the source data base

Next, we call BLAST to do the search.

``` bash
$ blastp -query blast_val.fasta -db blast_train.fasta -out blast_train_x_val.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

``` bash
$ blastp -query blast_test.fasta -db blast_train.fasta -out blast_train_x_test.tsv -outfmt "6 qseqid qlen sseqid slen length nident mismatch positive"
```

Here, `-outfmt` indicates how we specify the output format. The option
`6` will give us the output in tabular, with the supported format
specifiers:

-   qseqid means Query Seq-id
-   qlen means Query sequence length
-   sseqid means Subject Seq-id
-   slen means Subject sequence length
-   qstart means Start of alignment in query
-   length means Alignment length
-   nident means Number of identical matches
-   mismatch means Number of mismatches
-   positive means Number of positive-scoring matches

Results Analysis using R
------------------------

### Defining the function

``` r
# function to read the results and get the list of row index
blast_results <- function(result_path){

  # Read the results and turn it into dataframe
  # qseqid means Query Seq-id
  # qlen means Query sequence length
  # sseqid means Subject Seq-id
  # slen means Subject sequence length
  # length means Alignment length
  # nident means Number of identical matches
  # mismatch means Number of mismatches
  # positive means Number of positive-scoring matches
  df_results <- data.table::fread(result_path) %>%
    setNames(c("qseqid", "qlen", "sseqid", "slen", "length", "nident", "mismatch", "positive")) %>%
    rowwise() %>%
    mutate(
      percent_indentical = (nident/max(qlen, slen))*100, # The percentage of identical sequence over the longer sequence
      percent_positive = (positive/max(qlen, slen))*100 # The percentage of positive sequence over the longer sequence
      )
  
  # Get the data frame where the percent identical > 90
  df_identical_protein <- df_results %>% 
    filter(percent_indentical > 90)

  # Get the row indices of the source data for all of the identical percentage > 90%
  source_index_list_to_remove <- df_results %>%
    filter(percent_indentical > 90) %>%
    select(sseqid) %>%
    unique() %>%
    unlist()

  # Get the row indices of the query data for all of the identical percentage > 90%
  query_index_list_to_remove <- df_results %>%
    filter(percent_indentical > 90) %>%
    select(qseqid) %>%
    unique() %>%
    unlist()

  # Make list of all of the values
  list_results <- list(
    df = df_results,
    df_identical_protein = df_identical_protein,
    source_index = source_index_list_to_remove,
    query_index = query_index_list_to_remove
  )

  # Return the lists
  return(list_results)
}


# Function to remove from the list

remove_from_data_sets <- function(input_seq_path, label_seq_path, drop){

  # Read the data from the .csv file
  df_input <- data.table::fread(input_seq_path, header = FALSE) %>%
    setNames("sequence")
  df_label <- data.table::fread(label_seq_path, header = FALSE) %>%
    setNames("label")

  # Combine the input data and the label, then drop the rows based on the list
  df_new <- df_input %>%
    cbind(df_label) %>%
    filter(!row_number() %in% drop)
  
  # Combine the input data and the label, then take the data that should be removed
  removed_rows <- df_input %>%
    cbind(df_label) %>%
    filter(row_number() %in% drop)

  # Get the information about the removed columns
  removed_rows_freq <- df_label %>%
    filter(row_number() %in% drop) %>%
    select(label) %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("label", "freq")) %>%
    mutate(label = ifelse(label == 1, "effector", "noneffector"))

  # Create list for the results
  results_list <- list(
    df = df_new,
    removed_freq = removed_rows_freq, 
    removed_rows = removed_rows
  )

  return(results_list)
}


# Function to get information from all data
get_info <- function(blast_train_x_val_path,
                     blast_train_x_test_path,
                     input_train_path,
                     label_train_path,
                     input_val_path,
                     label_val_path,
                     input_test_path,
                     label_test_path){

  # Getting the results from comparing validation dataset and training dataset
  blast_results_train_x_val <- blast_results(blast_train_x_val_path)
  df_train_x_val <- blast_results_train_x_val[["df"]]
  source_index_train_x_val <- blast_results_train_x_val[["source_index"]]
  query_index_train_x_val <- blast_results_train_x_val[["query_index"]]

  # Getting the results from comparing testing dataset and training dataset
  blast_results_train_x_test <- blast_results(blast_train_x_test_path)
  df_train_x_test <- blast_results_train_x_test[["df"]]
  source_index_train_x_test <- blast_results_train_x_test[["source_index"]]
  query_index_train_x_test <- blast_results_train_x_test[["query_index"]]

  # Remove all of the rows of the training data
  # Since there might be intersection between index on the training when comparing them with validation and testing sets,
  # we need to find the intesections

  intersec_training_index <- c(source_index_train_x_val, source_index_train_x_test) %>%
    unique() %>%
    unlist()

  results_after_removed_training <- remove_from_data_sets(input_train_path, label_train_path, intersec_training_index)
  removed_rows_training_freq <- results_after_removed_training[["removed_freq"]] %>%
    mutate(type = "training")
  removed_rows_training <- results_after_removed_training[["removed_rows"]]

  results_after_removed_validation <- remove_from_data_sets(input_val_path, label_val_path, query_index_train_x_val)
  removed_rows_val_freq <- results_after_removed_validation[["removed_freq"]] %>%
    mutate(type = "validation")
  removed_rows_val <- results_after_removed_validation[["removed_rows"]]
  
  
  results_after_removed_testing <- remove_from_data_sets(input_test_path, label_test_path, query_index_train_x_test)
  removed_rows_test_freq <- results_after_removed_testing[["removed_freq"]] %>%
    mutate(type = "testing")
  removed_rows_test <- results_after_removed_testing[["removed_rows"]]

  all_removed_rows_freq <- removed_rows_training_freq %>%
    rbind(., removed_rows_val_freq) %>%
    rbind(., removed_rows_test_freq)
  
  all_removed_rows <- removed_rows_training %>%
    rbind(., removed_rows_val) %>%
    rbind(., removed_rows_test)

  all_ <- reshape2::dcast(all_removed_rows_freq, label ~ type, value.var = "freq") %>%
    as.data.frame()

  results_list <- list(
    all_val = blast_results_train_x_val,
    all_test = blast_results_train_x_test,
    removed_training = results_after_removed_training, 
    removed_validation = results_after_removed_validation, 
    removed_testing = results_after_removed_testing,
    all_removed_rows_freq = all_removed_rows_freq, 
    all_removed_rows = all_removed_rows,
    all = all_
  )

  return(results_list)
}
```

### Use the function defined

``` r
# define all the path of the results and the datasets in .csv
blast_train_x_val_path = "../../data/BLAST-data/blast_train_x_val.tsv"
blast_train_x_test_path = "../../data/BLAST-data/blast_train_x_test.tsv"
input_train_path = "../../data/BLAST-data/blast_train.csv"
label_train_path = "../../data/BLAST-data/blast_label_train.csv"
input_val_path = "../../data/BLAST-data/blast_val.csv"
label_val_path = "../../data/BLAST-data/blast_label_val.csv"
input_test_path = "../../data/BLAST-data/blast_test.csv"
label_test_path = "../../data/BLAST-data/blast_label_test.csv"


all_results <-  get_info(blast_train_x_val_path,
                         blast_train_x_test_path,
                         input_train_path,
                         label_train_path,
                         input_val_path,
                         label_val_path,
                         input_test_path,
                         label_test_path)
```

#### Results for comparing the validation data set and training dataset

Here is the result of the data that need to be removed since they are
considered as identical (90%):

``` r
df_iden_train_x_val <- all_results[["all_val"]][["df_identical_protein"]]

df_iden_train_x_val %>% 
  select(-c(mismatch, positive, percent_positive)) %>%
  rename("Query Seq-ID" = qseqid,
         "Query Length" = qlen,
         "Subject Seq-ID" = sseqid,
         "Subject Length" = slen,
         "Alignment Length" = length,
         "No Identical Matches" = nident,
         "% Identical Matches" =  percent_indentical
          ) %>%
  head(10) %>% 
  knitr::kable()
```

|  Query Seq-ID|  Query Length|  Subject Seq-ID|  Subject Length|  Alignment Length|  No Identical Matches|  % Identical Matches|
|-------------:|-------------:|---------------:|---------------:|-----------------:|---------------------:|--------------------:|
|             2|          2022|             109|            2015|              2024|                  2005|             99.15925|
|             2|          2022|             412|            2001|              2022|                  1998|             98.81306|
|             2|          2022|             184|            2023|              2047|                  1965|             97.13297|
|             3|          2057|             528|            2057|              2057|                  2048|             99.56247|
|             3|          2057|             106|            2035|              2035|                  2034|             98.88187|
|             6|          2290|             548|            2301|              2301|                  2097|             91.13429|
|            16|          1548|             222|            1548|              1546|                  1438|             92.89406|
|            26|          2048|             273|            2047|              2048|                  2040|             99.60938|
|            26|          2048|              84|            2048|              2048|                  2041|             99.65820|
|            27|          2339|              51|            2339|              2339|                  2339|            100.00000|

Below result is one of the example for pairwise for each sequence:

``` bash
Query= 2

Length=2022
                                                                      Score        E
Sequences producing significant alignments:                          (Bits)     Value

109                                                                   3934       0.0   
412                                                                   3915       0.0   
184                                                                   3823       0.0   
381                                                                   579        2e-173
221                                                                   576        4e-173
192                                                                   574        3e-172
455                                                                   570        5e-171
307                                                                   570        1e-170
437                                                                   569        2e-170
427                                                                   569        2e-170
353                                                                   569        2e-170
528                                                                   566        2e-169
106                                                                   565        5e-169
55                                                                    550        1e-163
520                                                                   540        4e-160
338                                                                   534        6e-158
573                                                                   524        1e-154
508                                                                   512        1e-150
423                                                                   512        1e-150
367                                                                   510        4e-150
524                                                                   478        2e-139
56                                                                    478        2e-139
420                                                                   262        2e-71 
285                                                                   262        3e-71 
38                                                                    261        4e-71 
100                                                                   249        2e-67 
253                                                                   239        4e-64 
58                                                                    234        9e-63 
115                                                                   232        3e-62 
116                                                                   231        6e-62 
235                                                                   230        1e-61 
47                                                                    230        1e-61 
175                                                                   230        1e-61 
568                                                                   225        5e-60 
120                                                                   223        2e-59 
364                                                                   221        1e-58 
44                                                                    202        3e-53 
284                                                                   190        3e-49 
241                                                                   188        9e-49 
548                                                                   187        2e-48 
494                                                                   186        5e-48 
533                                                                   178        8e-46 
398                                                                   178        9e-46 
512                                                                   157        4e-39 
372                                                                   151        1e-37 
532                                                                   146        5e-36 
478                                                                   136        4e-33 
213                                                                   117        3e-27 
559                                                                   110        4e-25 
311                                                                   110        4e-25 
417                                                                   70.5       6e-13 
340                                                                   68.6       3e-12 
193                                                                   55.8       2e-08 
451                                                                   48.9       2e-06 
43                                                                    39.7       0.002 
541                                                                   33.9       0.073 
444                                                                   32.7       0.15  
176                                                                   31.2       0.47  
521                                                                   29.6       1.5   
187                                                                   29.3       1.8   
431                                                                   27.3       6.9   
85                                                                    26.9       9.3   


>109
Length=2015

 Score = 3934 bits (10202),  Expect = 0.0, Method: Compositional matrix adjust.
 Identities = 2005/2024 (99%), Positives = 2005/2024 (99%), Gaps = 11/2024 (1%)

Query  1     MSFEDVLETCRERQIELWCDGEQLRYRAPAGALDAPLAARINAQRAAFVRFLRDGAWRAQ  60
             MSFEDVLETCRERQIELWCDGEQLRYRAPAGALDAPLAARINAQRAAFVRFLRDGAWRAQ
Sbjct  1     MSFEDVLETCRERQIELWCDGEQLRYRAPAGALDAPLAARINAQRAAFVRFLRDGAWRAQ  60

Query  61    PSRTHERFALTPVQAAYVLGRHAAFEFGGSACHLYVEYRESADLDVRRFEAAWNACVARH  120
             PSRTHERFALTPVQAAYVLGRHAAFEFGGSACHLYVEYRESADLDVRRFEAAWNACVARH
Sbjct  61    PSRTHERFALTPVQAAYVLGRHAAFEFGGSACHLYVEYRESADLDVRRFEAAWNACVARH  120

Query  121   PMLRAIVEDNAWQRILPDVPWQTLAVHDLRDANAAAFDAHVRRVRERLDHAVHALDRWPV  180
             PMLRAIVEDNAWQRILPDVPWQTLAVHDLRDANAAAFDAHVRRVRERLDHAVHALDRWPV
Sbjct  121   PMLRAIVEDNAWQRILPDVPWQTLAVHDLRDANAAAFDAHVRRVRERLDHAVHALDRWPV  180

Query  181   LQPEVSIGPDGAVLHLSVDFTLIDYASLQLLLAEWRRRHDDPQWKPAPLDVTFRDYVVNE  240
             LQPEVSIGPDGAVLHLSVDFTLIDYASLQLLLAEWRRRHDDPQWKPAPLDVTFRDYVVNE
Sbjct  181   LQPEVSIGPDGAVLHLSVDFTLIDYASLQLLLAEWRRRHDDPQWKPAPLDVTFRDYVVNE  240

Query  241   ARAGRRAQHARDAAWWLARIDGLPGRPDLPVLPRPPEGRADARPRFTHRHARLDRARWDA  300
             ARAGRRAQHARDAAWWLARIDGLPGRPDLPVLPRPPEGRADARPRFTHRHARLDRARWD 
Sbjct  241   ARAGRRAQHARDAAWWLARIDGLPGRPDLPVLPRPPEGRADARPRFTHRHARLDRARWDV  300

Query  301   LVAFASRFGLSPAGVALAAFAEVVGRWSQSPAFCLNLTVLNRPPVHPQIDAVLGDFTALS  360
             LVAFASRFGLSPAGVALAAFAEVVGRWSQSPAFCLNLTVLNRPPVHPQIDAVLGDFTALS
Sbjct  301   LVAFASRFGLSPAGVALAAFAEVVGRWSQSPAFCLNLTVLNRPPVHPQIDAVLGDFTALS  360

Query  361   LLAVDVAQGRDFAERARAIGAQMFDDLDHRAFTGVEVMRELARRRGKDAALMPVVFTSGI  420
             LLAVDVAQG DFA RARAIGAQMFDDLDHRAFTGVEVMRELARRRGKDAALMPVVFTSGI
Sbjct  361   LLAVDVAQGSDFAARARAIGAQMFDDLDHRAFTGVEVMRELARRRGKDAALMPVVFTSGI  420

Query  421   GSVGRLLGEHGEHGARAAPPCYMISQTPQVWLDCQVTDQFGGLEIGWDVRDDLFPAGMPA  480
             GSVGRLLGEHGEHGARAAPPCYMISQTPQVWLDCQVTDQFGGLEIGWDVRDDLFPAGMPA
Sbjct  421   GSVGRLLGEHGEHGARAAPPCYMISQTPQVWLDCQVTDQFGGLEIGWDVRDDLFPAGMPA  480

Query  481   AMFDAYVALLARLASDAAWWARPGDIVLPSGPPVACAHPDADRHLAAGFAAQARRTPDAT  540
             AMFDAYVALLARLASDAAWWARPGDIVLPSGPPVACAHPDADRHLAAGFAAQARRTPDAT
Sbjct  481   AMFDAYVALLARLASDAAWWARPGDIVLPSGPPVACAHPDADRHLAAGFAAQARRTPDAT  540

Query  541   VVIDAAGAHTYRDVAQRAAAVRAALERAGVAPGDKVAVRMPKGANQLVAVLGIVQAGAAY  600
             VVIDAAGAHTYRDVAQRAAAVRAALERAGVAPGDKVAVRMPKGANQLVAVLGIVQAGAAY
Sbjct  541   VVIDAAGAHTYRDVAQRAAAVRAALERAGVAPGDKVAVRMPKGANQLVAVLGIVQAGAAY  600

Query  601   VPIDYRQPALRRRAILRNAQVGAIVTERALDGEPDGCARIDVDALAPDPRWPPRDAHPLE  660
             VPIDYRQPALRRRAILRNAQVGAIVTERALDGEPDGCARIDVDALAPDPRWPPRDAHPLE
Sbjct  601   VPIDYRQPALRRRAILRNAQVGAIVTERALDGEPDGCARIDVDALAPDPRWPPRDAHPLE  660

Query  661   GDALAYVIYTSGSTGEPKGVMVSHAAVCNTLADINARHAVGAGDAVLGLAELSFDLSVYD  720
             GDALAYVIYTSGSTGEPKGVMVSHAAVCNTLADINARHAVGAGDAVLGLAELSFDLSVYD
Sbjct  661   GDALAYVIYTSGSTGEPKGVMVSHAAVCNTLADINARHAVGAGDAVLGLAELSFDLSVYD  720

Query  721   LFGATARGARVVLPDPARGNDPSHWAELIARHGVTLWNSVPAQGQMLIDYLESEPARAMP  780
             LFGATARGARVVLPDPARGNDPSHWAELIARHGVTLWNSVPAQGQMLIDYLESEPARAMP
Sbjct  721   LFGATARGARVVLPDPARGNDPSHWAELIARHGVTLWNSVPAQGQMLIDYLESEPARAMP  780

Query  781   GPRCVMWSGDWIPVSLPTRWWRRWPDSRLFSLGGATEASIWSIEHPIRPEDTRLVSIPYG  840
             GPRCVMWSGDWIPVSLPTRWWRRWPDSRLFSLGGATEASIWSIEHPIRPEDTRLVSIPYG
Sbjct  781   GPRCVMWSGDWIPVSLPTRWWRRWPDSRLFSLGGATEASIWSIEHPIRPEDTRLVSIPYG  840

Query  841   RALTGQTVDVLDALGRPCPPGVRGEIHIGGVGLATGYANDPARTAERFVRHADGRRLYRT  900
             RALTGQTVDVLDALGRPCPPGVRGEIHIGGVGLATGYANDPARTAERFVRHADGRRLYRT
Sbjct  841   RALTGQTVDVLDALGRPCPPGVRGEIHIGGVGLATGYANDPARTAERFVRHADGRRLYRT  900

Query  901   GDLGRRRADGSLEFLGRQDDQVKIRGYRIELAEIDAALSAHPRVASAATIVLGDAAQRRL  960
             GDLGRRRADGSLEFLGRQDDQVKIRGYRIELAEIDAALSAHPRVASAATIVLGDAAQRRL
Sbjct  901   GDLGRRRADGSLEFLGRQDDQVKIRGYRIELAEIDAALSAHPRVASAATIVLGDAAQRRL  960

Query  961   ASFVTLHGAAPDPRRRDAQLLAIAQRVRDALAAERWPARAEIRRSVGQLDAACVASLASW  1020
             ASFVTLHGAAPDPRRRDAQL AIAQRVRDALAAERWPARAEIRRSVGQLDAACVASLASW
Sbjct  961   ASFVTLHGAAPDPRRRDAQLHAIAQRVRDALAAERWPARAEIRRSVGQLDAACVASLASW  1020

Query  1021  IVRAGGLTRDASVDFATLAERLRVPAARHGLLRHWLALLEARGVLCTDTAGPVAGVDVTE  1080
             IVRAGGLTRDASVDFATLAERLRVPAARHGLLRHWLALLEARGVLCTDTAGPVAGVDVTE
Sbjct  1021  IVRAGGLTRDASVDFATLAERLRVPAARHGLLRHWLALLEARGVLCTDTAGPVAGVDVTE  1080

Query  1081  SENSAASRGAADADAFEAFEAFEAFEAFEAFEAFEAFEAFEAFEAFEAPEAPEAPEAPEA  1140
             SENSAASRGAADAD      AFEAFEAFEAFEAFEAFE FEAFEA EAPEAPEAPEAPEA
Sbjct  1081  SENSAASRGAADAD------AFEAFEAFEAFEAFEAFEVFEAFEAPEAPEAPEAPEAPEA  1134

Query  1141  PEAPEAPEAPEAPEAPEAPEATEATEATEATEAPEATEATEATEAIETIETIETIETIET  1200
             PEAPEAPEAPEAPEAPEAPEATEATEATEATEAPEA EATEATEA   IETIETIETIET
Sbjct  1135  PEAPEAPEAPEAPEAPEAPEATEATEATEATEAPEAPEATEATEA---IETIETIETIET  1191

Query  1201  IETIETIETIETPNAADVTKATKSAGTTDTTTSPAGWRVARAHAARATDAAGAVPLSDAD  1260
             IETIETIETIETPNAADVTKATKSAGTTDTTTSPAGWRVARAHAARATDAAGAVPLSDAD
Sbjct  1192  IETIETIETIETPNAADVTKATKSAGTTDTTTSPAGWRVARAHAARATDAAGAVPLSDAD  1251

Query  1261  ACWDAFARDASPALWPAVLIDYFRASAACLGEQVDDRVSPAGLMFPQGAPHVAEAMYSDG  1320
             ACWDAFARDASPALWPAVLIDYFRASAACLGEQVDDRVSPAGLMFPQGAPHVAEAMYSDG
Sbjct  1252  ACWDAFARDASPALWPAVLIDYFRASAACLGEQVDDRVSPAGLMFPQGAPHVAEAMYSDG  1311

Query  1321  MHARALHRAMAEAVSGIVAREPRRAWRILEIGAGTAAATRAIVDALAPHAESGVRIDYLF  1380
             MHARALHRAMAEAVSGIVAREPRRAWRILEIGAGTAAATRAIVDALAPHAESGVRIDYLF
Sbjct  1312  MHARALHRAMAEAVSGIVAREPRRAWRILEIGAGTAAATRAIVDALAPHAESGVRIDYLF  1371

Query  1381  SDVSSYFLAAARERFAAYPWVRFIRFDMNAPLDAQGVAPHSVDLIASSGALNNARDTVAL  1440
             SDVSSYFLAAARERFAAYPWVRFIRFDMNAPLDAQGVAPHSVDLIASSGALNNARDTVAL
Sbjct  1372  SDVSSYFLAAARERFAAYPWVRFIRFDMNAPLDAQGVAPHSVDLIASSGALNNARDTVAL  1431

Query  1441  VAGLRALSSADAWWVIQELTAEHPEISISQGLMMEPPRDARAASHRLFVHRAQWLEWLRA  1500
             VAGLRALSSADAWWVIQELTAEHPEISISQGLMMEPPRDARAASHRLFVHRAQWLEWLRA
Sbjct  1432  VAGLRALSSADAWWVIQELTAEHPEISISQGLMMEPPRDARAASHRLFVHRAQWLEWLRA  1491

Query  1501  AGDRALGCVEPDSPLDALGYDILLARVKADAARVDADELLAFAAERVPRYMVPSQLCALE  1560
             AGDRALGCVEPDSPLDALGYDILLARVKADAARVDADELLAFAAERVPRYMVPSQLCALE
Sbjct  1492  AGDRALGCVEPDSPLDALGYDILLARVKADAARVDADELLAFAAERVPRYMVPSQLCALE  1551

Query  1561  RLPVTANGKIDRRALAQIAGVHEPRAAHGRAVRPCDAEPLLARLIGLWEAVLDTQGVAAD  1620
             RLPVTANGKIDRRALAQIAGVHEPRAAHGRAVR CDAEPLLARLIGLWEAVLDTQGVAAD
Sbjct  1552  RLPVTANGKIDRRALAQIAGVHEPRAAHGRAVRACDAEPLLARLIGLWEAVLDTQGVAAD  1611

Query  1621  QDFFAAGGDSLLIAQLVSRVRAEEPLARAHPFDRLLRWALAQPTPAGLAQCLRDAKAQAQ  1680
             QDFFAAGGDSLLIAQLVSRVRAEEPLARAHPFDRLLRWALAQPTPAGLAQCLRDAKAQAQ
Sbjct  1612  QDFFAAGGDSLLIAQLVSRVRAEEPLARAHPFDRLLRWALAQPTPAGLAQCLRDAKAQAQ  1671

Query  1681  SSTAAPSR--TATATATATATATATAATAAARAAPAASIVHAPAPRGRIDVRPLRLAPGV  1738
             SSTAAPSR  TATATATATATATATAATAAARAAPAASIVHAPAPRGRIDVRPLRLAPGV
Sbjct  1672  SSTAAPSRTATATATATATATATATAATAAARAAPAASIVHAPAPRGRIDVRPLRLAPGV  1731

Query  1739  GVPRVVVHEGLGTVHAYRPIVPALARLGPVLGFAVRDAQDYLDLPAQHLSATLGARYAYA  1798
             GVPRVVVHEGLGTVHAYRPIVPALARLGPVLGFAVRDAQDYLDLPAQHLSATLGARYAYA
Sbjct  1732  GVPRVVVHEGLGTVHAYRPIVPALARLGPVLGFAVRDAQDYLDLPAQHLSATLGARYAYA  1791

Query  1799  LSREGIGEVDVLGYCSGGLIALEMAKTLVQLGVAVRSLDIVSSYRIPYLIEDERLVLFNF  1858
             LSREGIGEVDVLGYCSGGLIALEMAKTLVQLGVAVRSLDIVSSYRIPYLIEDERLVLFNF
Sbjct  1792  LSREGIGEVDVLGYCSGGLIALEMAKTLVQLGVAVRSLDIVSSYRIPYLIEDERLVLFNF  1851

Query  1859  AATLGLPLDALGFPAAHLLADALADALKADPARVAPGSIQAQLEAFGDRCAPLDILRRRV  1918
             AATLGLPLDALGFPAAHLLADALADALKADPARVAPGSIQAQLEAFGDRCAPLDILRRRV
Sbjct  1852  AATLGLPLDALGFPAAHLLADALADALKADPARVAPGSIQAQLEAFGDRCAPLDILRRRV  1911

Query  1919  LRVSAGLPADGDEPHPLVDERERLYRLFMHSVHASHWAAHAPYAGPLRLFVPERCNPLIP  1978
             LRVSAGLPADGDEPHPLVDERERLYRLFMHSVHASHWAAHAPYAGPLRLFVPERCNPLIP
Sbjct  1912  LRVSAGLPADGDEPHPLVDERERLYRLFMHSVHASHWAAHAPYAGPLRLFVPERCNPLIP  1971

Query  1979  QQRAALFDYWRDQALGGIALVDVPGGHFDCLTAAFVDTHLKEAR  2022
             QQRAALFDYWRDQALGGIALVDVPGGHFDCLTAAFVDTHLKEAR
Sbjct  1972  QQRAALFDYWRDQALGGIALVDVPGGHFDCLTAAFVDTHLKEAR  2015
```

And the following table shows the identical protein between training
data sets and testing datasets.

``` r
df_iden_train_x_test <- all_results[["all_test"]][["df_identical_protein"]]

df_iden_train_x_test %>% 
  select(-c(mismatch, positive, percent_positive)) %>% 
  rename("Query Seq-ID" = qseqid,
         "Query Length" = qlen,
         "Subject Seq-ID" = sseqid,
         "Subject Length" = slen,
         "Alignment Length" = length,
         "No Identical Matches" = nident,
         "% Identical Matches" =  percent_indentical
          ) %>% 
  head(10) %>% 
  knitr::kable()
```

|  Query Seq-ID|  Query Length|  Subject Seq-ID|  Subject Length|  Alignment Length|  No Identical Matches|  % Identical Matches|
|-------------:|-------------:|---------------:|---------------:|-----------------:|---------------------:|--------------------:|
|             1|          1926|             226|            1926|              1926|                  1926|            100.00000|
|             4|          1126|              81|            1163|              1163|                  1108|             95.27085|
|             4|          1126|             346|            1164|              1164|                  1094|             93.98625|
|             4|          1126|               3|            1096|              1130|                  1063|             94.40497|
|            14|          2059|             350|            2059|              2059|                  2059|            100.00000|
|            15|          2260|             340|            2260|              2260|                  2256|             99.82301|
|            15|          2260|             417|            2260|              2260|                  2231|             98.71681|
|            16|           125|               8|             126|               124|                   115|             91.26984|
|            26|          2273|              60|            2273|              2282|                  2255|             99.20810|
|            26|          2273|              96|            2323|              2332|                  2253|             96.98666|

Lastly, we can get the result of the summary of total number of sequence
that has to be removed from each data sets:

``` r
summary <- all_results[["all"]]

summary %>% 
  knitr::kable()
```

| label       |  testing|  training|  validation|
|:------------|--------:|---------:|-----------:|
| effector    |       11|        21|           6|
| noneffector |       28|        76|          36|

``` r
all_removed_rows <- all_results[["all_removed_rows"]]

removed_all_effector <- all_removed_rows %>% 
  filter(., label == 1)

removed_all_noneffector <- all_removed_rows %>% 
  filter(., label == 0)

removed_all_noneffector %>% 
  nrow()
```

    ## [1] 140

View how the data looks like

``` r
bacteria <- c("Acinetobacter baumannii", 
"Aeromonas hydrophila", 
"Aeromonas salmonicida", 
"Brucella abortus", 
"Burkholderia glumae", 
"Burkholderia pseudomallei", 
"Campylobacter jejuni", 
"Citrobacter rodentium", 
"Clavibacter michiganensis", 
"Coxiella burnetii", 
"Cystobacter fuscus", 
"Edwardsiella ictaluri", 
"Erwinia amylovora", 
"Escherichia coli", 
"Francisella tularensis", 
"Helicobacter pylori", 
"Legionella pneumophila",  
"Listeria monocytogenes", 
"Mycobacterium tuberculosis", 
"Pantoea stewartii", 
"Pseudomonas aeruginosa", 
"Pseudomonas cichorii", 
"Pseudomonas savastanoi", 
"Pseudomonas syringae", 
"Ralstonia solanacearum", 
"Salmonella enterica", 
"Shigella flexneri", 
"Staphylococcus aureus", 
"Vibrio parahaemolyticus", 
"Xanthomonas axonopodis", 
"Xanthomonas campestris", 
"Xanthomonas oryzae", 
"Xylella fastidiosa", 
"Yersinia enterocolitica", 
"Yersinia pestis", 
"Yersinia pseudotuberculosis")

fungi <- c("Beauveria bassiana",
"Blumeria graminis",
"Botrytis cinerea",
"Cercospora apii",
"Cercospora beticola",
"Colletotrichum orbiculare",
"Dothistroma septosporum",
"Fusarium oxysporum",
"Leptosphaeria maculans",
"Magnaporthe oryzae",
"Melampsora lini",
"Parastagonospora nodorum",
"Passalora fulva",
"Penicillium expansum",
"Pseudocercospora fuligena",
"Puccinia striiformis",
"Rhynchosporium commune",
"Verticillium dahliae",
"Ustilago maydis",
"Zymoseptoria tritici")

oomycetes <- c("Hyaloperonospora arabidopsidis", 
"Phytophthora cactorum", 
"Phytophthora capsici", 
"Phytophthora infestans", 
"Phytophthora parasitica", 
"Phytophthora sojae", 
"Pythium aphanidermatum", 
"Plasmopara halstedii", 
"Phytophthora megakarya")
```

``` r
effector_data <- data.table::fread("../../data/effector_with_IDs_organism.csv")

effector_removed_with_species <- removed_all_effector %>%
  left_join(effector_data, by = "sequence") %>% 
  select(sequence, label, pathogen_short) %>% 
  mutate(category = ifelse(pathogen_short %in% bacteria, "bacteria", 
                    ifelse(pathogen_short %in% fungi, "fungus", 
                    ifelse(pathogen_short %in% oomycetes, "oomycetes", "others")))) %>% 
  group_by(category) %>% 
  summarise(count = n())

noneffector_data <- data.table::fread("../../data/noneffector_with_IDs_organism.csv")

noneffector_removed_with_species <- removed_all_noneffector %>%
  left_join(noneffector_data, by = "sequence") %>% 
  unique() %>% 
  select(sequence, label, pathogen_short) %>%
  mutate(category = ifelse(pathogen_short %in% bacteria, "bacteria",
                    ifelse(pathogen_short %in% fungi, "fungus",
                    ifelse(pathogen_short %in% oomycetes, "oomycetes", "others")))) %>% 
  group_by(category) %>% 
  summarise(count = n())

effector_removed_with_species %>% 
  knitr::kable()
```

| category  |  count|
|:----------|------:|
| bacteria  |     27|
| fungus    |     11|
| oomycetes |      6|

``` r
noneffector_removed_with_species %>% 
  knitr::kable()
```

| category  |  count|
|:----------|------:|
| bacteria  |    122|
| fungus    |      2|
| oomycetes |     16|

Identify the Classification of the sequences that are removed
-------------------------------------------------------------
