Daily lab book April
====================

10 April 2019 (Wednesday)
-------------------------

-   met Dan for the research plan
-   tried to complete the script of hyperparameters scan
-   copied some protein IDs

11 April 2019 (Thursday)
------------------------

-   copied manually the protein IDs until get the idea to find the `phi-base` github repo
-   started to clean the data to get the protein IDs
-   after got the protein IDs, then used it to get the sequence data


12 April 2019 (Friday)
----------------------

-   analysed and cleaned the sequence data from `uni-prot.org`

15 April 2019 (Monday)
----------------------

- finished reading Dan's paper
- tried to complete hyperparameter scan Model 1

**Minutes of weekly catch up meeting**

- finish the hyperparameter scan! (time management and prioritizing task skill should be improved!)
- next task: get the data of non-effector from NCBI 

16 April 2019 (Tuesday)
------------------------

- finished hyper scan Model 1 and 2 scripts (run in GPU)

17 April 2019 (Wednesday)
-------------------------

- tried to finish the Model 3 scripts
- getting the non-effector data from NCBI and phi-base

18 April 2019 (Thursday -- Half Day)
------------------------------------

- continue with getting the non-effector data

**Note**

- the first plan is to get the non-effector from `phi-base.org`, but some of the organism does not have sample non-effector data in `phi-base.org`. Solution: getting the sample from NCBI (using the query)

- some of the organism names in `uniprot` have different name with some in `phi-base.org`, then I did mapping to get the original organism name. For example the protein ID `G2XWG3` in `phi-base`: *Botrytis cinerea*, whereas in  `uniprot`: *Botryotinia fuckeliana*. And when we search the non-effector in `phi-base` using *Botryotinia fuckeliana*, we will not find it, because it is not there. Then we need to map it back to get the original name. (The `phi-base` data is the main data for finding the non-effector data)

23 April 2019 (Tuesday)
------------------------

- obtained the intersect pathogen name of effector data from phi-base and uni-prot, identified the non-intersect one, and replace with the pathogen name from the phi-base.
- added the running time in SBATCH scripts in tsl-gpu (one job got canceled bc of the time limit)

24 April 2019 (Wednesday)
--------------------------

- focus on getting the non-effector data from phi-base big data



