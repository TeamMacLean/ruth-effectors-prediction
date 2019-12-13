
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

## 05, 06 December 2019 (Thu, Fri)

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

## 10, 11, 12 December 2019 (Tue, Wed, Thu)

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

Commits: 0

#### Worked on

-----
