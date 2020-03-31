
# December 2019 Lab Report <img src="figures/tsl-logo.png" align="right" width="120" />

### [Ruth Kristianingsih](https://github.com/ruthkr)

## 02 December 2019 (Monday)

Commits: 0

#### Worked on

  - Heatmap and reports.

#### Weekly catch up meeting

##### Effector and non-effector prediction

  - Based on the results of CNN-LSTM and CNN-GRU, the model can learn
    quite good for each separated, and they were not optimised models,
    then we can find the optimised one using hyperparameter scans for
    each

  - About the strange behaviour on CNN-LSTM layer on fungi data, check
    again the code, it might be simple mistake on the code (for instance
    can be indexing), also can check the sequence length of fungi data

  - using the signalP, try to get the secreted protein sequences of
    organisms and include them on the datasets (also train them on the
    current models, just to see if the convolutional network actually
    recognises the signal peptide of the sequences)

## 03 December 2019 (Tuesday)

Commits: 0

#### Worked on

  - Created some scripts of effector prediction and run in `tsl-gpu` and
    `tsl-long`:
    
      - CNN-LSTM models for each bacteria, fungi, and oomycete data.
      - CNN-GRU models for each bacteria, fungi, and oomycete data.

  - Installed the SignalP 5.0, and read the instruction.

## 04 December 2019 (Wednesday)

Commits: 0

#### Worked on

  - All jobs killed at 9 am due to the service discruption.
  - Run all jobs again, and the jobs got killed again at around 12 pm
    (HPC is down).
  - Read about the secretory pathway protein - the process and the
    relationship between effector protein and secreted protein.
  - Checking the bug on my heatmap code.
      - Check the length of fungi sequences (apparently, the length of
        fungi dataset is very low, see the
        0010\_getting\_info\_data\_each\_class.Rmd), so its definitely
        the problem with code or method.

## 05, 06 December 2019 (Thursday, Friday)

Commits: 0

#### Worked on

  - Focus to see what is wrong with my code, I have found the bugs but
    seems it is not the only bugs I have.
  - Created a function in R to run SignalP.

-----

## 09 December 2019 (Monday)

Commits: 0

#### Weekly catch up meeting

##### Effector and non-effector prediction

  - Finding bugs:
    
      - Try to put everything in one function to maintain consistency
        (DRY-Dont Repeat Yourself), use docstring so that if one try to
        get `help()`, it will print the documentation written in
        docstring.
      - Try to use Unit-test.
      - Make a python module - to organise python function in one file.

  - Regarding to the protein data:
    
      - Use SignalP 2 not signalP 5.
      - Try to get secreted data forthe organism.

## 10, 11, 12 December 2019 (Tuesday, Wednesday, Thursday)

Commits: 0

#### Worked on

  - Downloaded SignalP 2.0, but I could not install it on my laptop.
    There is no SignalP 2.0 in cluster.
  - Learning to understand properly how to get protein from the same
    genome in NCBI data (tried to undestand the process, and if it is
    possible to retrieve the data esily)
  - Wrote a script to retrive the protein data using
    `library("rentrez")`, I tried created create a script to access the
    NCBI site, and read the html file using `readLines()`, but I found
    that using `rentrez` library is easier.

## 13 December 2019 (Friday)

Commits: 2

#### Worked on

  - investigated the bug on heatmap scripts (found that only conv1d\_4
    behaves differently, Ho: this is due to the fact that last
    convolutional layer is concatenation layer)
  - continued the process on getting the secreted protein data for
    different cases where the protein data does not have the Gene term
    and also if the protein data does not have genomic sequence in NCBI
    database. The script report is available on
    [0001\_get\_secreted\_data.md](https://github.com/TeamMacLean/ruth-effectors-prediction/blob/master/reports/getting-data-secreted/0001_get_secreted_data.md)

## 16 December 2019 (Monday)

Commits: 9

#### Worked on

  - completing report of heatmaps investigation for the strange
    behaviour on the last concatenation layer (`conv1d_4`)
    
      - [Heatmaps for CNN-GRU
        model](https://github.com/TeamMacLean/ruth-effectors-prediction/blob/master/scripts/jupyter-note/heatmaps/heatmap_cnn_gru.ipynb)
      - [Heatmaps for CNN-LSTM
        model](https://github.com/TeamMacLean/ruth-effectors-prediction/blob/master/scripts/jupyter-note/heatmaps/heatmap_cnn_lstm.ipynb)
      - [Heatmaps for fungi data on different layers of CNN-LSTM
        models](https://github.com/TeamMacLean/ruth-effectors-prediction/blob/master/scripts/jupyter-note/heatmaps/fungi-test/heatmap_cnn_lstm_fungi_test.ipynb)

  - updated new function to retrive protein datasets

## 17 December 2019 (Tuesday)

Commits: 1

#### Worked on

  - added the function to get protein data from NCBI using geneID and
    get the genome sequence ID ===\> this is wrong step, needs to get
    the data from EnsEMBL instead and also it has signal Peptide
    prediction already, so that perhaps does not need to use SignalP to
    predict the signal peptide
  - tried to install SignalP 2.0 with Martin’s help, but it can’t be
    installed (the version is too old)
  - installing the biomaRt package and reading the documentation of
    biomaRt
  - understanding how the ensEMBL database works for each organism

#### Weekly catch up meeting

##### Effector and non-effector prediction

  - For the transfer learning, Dan agreed that I can go in getting data
    from previous research about effector prediction (effectorP for
    fungi), and others.  
  - Think about what it means with why some convolutional layers are
    differently activated.
  - For retrieving data or in general, keep it in mind that if what we
    are doing is not specific in this topic, then there must be
    something has been developed before (need to develop a sense like
    this) – ensEmbl database by EMBL.

## 18 December 2019 (Wednesday)

Commits: 1

  - apparently we can’t use `biomaRt` package in R in ensEmbl genomes,
    we need to write query to retrieve the protein data
  - `biomaRt` can be only used for `ensemble.org`

## 19, 20 December 2019 (Thursday, Friday)

Commits: 0

  - retrieve data using from ensEmbl using bioMart
  - the bioMart tool can be used only for retrieving the sequences one
    by one, thats why I need to do query in order to retrieve the data

-----
