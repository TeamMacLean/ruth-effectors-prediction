Daily lab book June
====================


03 - 04 June 2019 (Monday - Tuesday)
------------------------------------

- made report of the models for both effector--noneffector and subcellular protein data
- wrote scripts hyperparameter scans for effector and noneffector protein data


05 June 2019 (Wednesday)
-------------------------

- had meeting at 10 am
- made the hyper-parameter scripts for the Conv1D model 

**Weekly catch up minutes**

- Project 1 (Subcellular Localisation)

    - Try to think alternative of other encoding method for the amino acids protein sequence (K-mer, CGR, or n-grams from K-mer)
    - For the CGR, due to the memory limitation, then probably the best mac depth is 5
    
- Project 2 (Effector - Noneffector protein)

    - Since the results are too good, too good to be true, I need to do further analysis on my data sets
    - Present to Dan for next meeting about the training data set I have
    - I can do the bootstrap random to select control set for the non-effector data sets. 

06 June 2019 (Thursday)
------------------------

- finished the hyperparameter scan scripts for Conv1D + LSTM
- find out other options of encoding amino acids
- tried to download the data as non-effector control sets

07 June 2019 (Friday)
------------------------

- testing new models for subcellular localisation if it is using more filter
- read papers about possible encoding methods
- getting data for control sets

10 June 2019 (Monday)
----------------------

- made the report
- continued the encoding parts of the control datasets
- read about the encoding method alternatives

**Weekly catch up minutes**

- Project 1

    - trying the parallel running: reduce the hyperparameters into several groups
    - using intuition: pick one of the encoding method that is do-able
    - if I want to try something with model: read about benchmarking for code run
    
- Project 2

    - get again the data from NCBI not from from Phi-base
    - the next step is doing multiclass classification
    - recognizing signal peptide carrying using `SignalP`
    
11, 12 June 2019 (Tuesday and Wednesday)
-----------------------

- retrieve data 
- Chaos Game Representation (CGR)
- investigate the encoding method that may be do-able

13 June 2019 (Thursday)
-----------------------

- sick day

14 June 2019 (Friday)
----------------------

- retrieving data

17 -- 21 June 2019 (A week)
----------------------

- continue the report of retrieving data 
- getting data from the beginning. (All was wrong since the first step of filtering)
- making the parallism code for the subcellular localisation

**Weekly catch up minutes**

- Project 1: prioritize the parallelism (naively) -- may need CSV file to break down the hyperparameters, CGR, and getting one encoding method 

- Project 2: present to Dan better about the data obtained

24 -- 28 June 2019 (A week)
-------------------------------

- getting again effector data (step by step)
- made report of the best model in the test data (not inlcuded in the previous scripts)


