Daily lab book May
==================

01 May 2019 (Wednesday)
----------------------

- figured out how to get the information / data from the `.log` file of CNN model using sep convolutional 1D


02 - 03 May 2019 (Thursday and Friday)
----------------------

- continue the extract the information from the `.log` file
- getting the non-effector data from the NCBI (per organism, put inside the collection)

06 May 2019 (Monday -- Bank holiday)
------------------------------------

07 May 2019 (Tuesday)
---------------------

- add the verbose to get more info in the log file of the scripts that will be run in the GPU cluster
- install BLAST and learn how it works

**Weekly catch up minutes**

- Make sure focus on what the objectives are (not to overdoing something)
- While waiting for the first model run to finish, I can create another model for sub-cellular localisation (LSTM + CNN model)
- Get the training, validation, and testing data. Compare the sequence in those set of data using BLAST, and report to Dan on Monday which data will get rid of. 

08, 09May 2019 (Wednesday and Thursday)
-----------------------

- play a lil bit with BLAST (try using some data in a tutorial)
- started to work with my data 

10 May 2019 (Friday)
--------------------

- understanding the results and made the report using R

13 May 2019 (Monday)
--------------------

- learned RNN from `Deep Learning with Python`
- tried to retrieve data from NCBI

**Weekly catch up meeting's minutes**

- Draft model will be expected next monday
- Continue finding sequences using API / Batch Entrez. Can't be done using SRA

14 May 2019 (Tuesday)
---------------------

- Effector Classification Project

    - specified the categories of data that have been removed from the datasets (whether it is bacteria, fungi, or oomycetes)
    - get new data from NCBI from Bacteria
    

- Subcellular Localisation Project

    - learned about word Embedding, sequence preprocessing in R
    - SimpleRNN, and the concepts of LSTM
    
15, 16, 17 May 2019 (Wednesday, Thursday. and Friday)
------------------------------------------------------

- read and learn about LSTM, GRU architectures
- create the draft of architecture models
- added the additional effector and noneffector data to dataframe without identical protein sequence
- check the models that have been running in GPUs

20 May 2019 (Monday)
--------------------

- updated the report
- checked how much progress that has been achieved of the model running in cluster

**Weekly catch up meeting's minutes**

- think about grow structure, think carefully how each layer actually impact on the data and results (end to end process) >> balancing the layers
- try to work with epoch max 40 (in the range 20 -- 40)
- it is better to have simple model that run quickly and can be analysed and give results quickly

**Next steps**

- create simpler model with considering the grow structure (for CNN and also CNN + LSTM)
- start to encode the data that has been ready 


