
# February 2020 Lab Report <img src="figures/tsl-logo.png" align="right" width="120" />

### [Ruth Kristianingsih](https://github.com/ruthkr)

## 03 February 2020 (Monday)

Commits: 4

#### Weekly meeting minutes:

  - need to run all of the four models we got for all different datasets
    and see the results
  - after getting all of the results, then find the correlation between
    the prediction result for each models using
    Pearson-correlation-matrix ==\> this will summarise the behaviour of
    the models (whether they are overlapping of completely separated)
  - get another ensemble model (which using voting classifier,
    consenseus)
  - run the multi-class classification data on the CNN-GRU model
  - putting sigmoid in after calculating weighted average

## 04–07 February 2020

Commits: 0

#### Worked on

  - created scripts for all models in different datasets : `LSTM_emb`
    and `GRU-emb` and get it running in cluster
  - made report once models’ training done
  - created a script for CNN-GRU for multi-class classification
    subcellular-localisation
  - created a function to make the pearson correlation matrix (ggplot)
  - made function of weighted average and voting classifier ensemble

-----

## 10 February 2020 (Monday)

#### Activities:

  - Going for a German visa application in London.

## 11 February 2020 (Tuesday)

Commits: 4

#### Weekly meeting minutes:

  - The next step of the model training is to manually tuning each model
    and plot using graph to see if its overfitting or not (comparing the
    training and validation data)
  - the sigmoid function of the data is with range\[0,1\] not \[-4, 4\]
    ==\> shifted sigmoid
  - for the pearson correlation matrix, can you put the ensemble one on
    the corner?
  - working also on the project report, since we want to publish it
  - try to get another datasets from previous research done in the same
    topics (bacteria and fungi)

## 12 February 2020 (Wednesday)

Commits: 0

#### Worked on:

  - reset some commit
  - put the sigmoid function for the weighted average ensemble function
    (to make the function smoother) ==\> using logistic function
  - run again lstm\_emb scan oomycete that failed earlier before,
    because the mistake on the unconsistent name of the data
  - make the ensemble in the corner of Pearson correlation matrix
    functions ==\> added a little code after pivotting the data
  - started the scripts for the manual tuning
  - Upadates of the job running: there is no job finished yet
  - updated the lab book for February 2020

## 13 February 2020 (Thursday)

Commits: 0

#### Worked on:

  - put the sigmoid function for the weighted average ensemble method so
    that it will make the value smoother
  - run again the lstm embedding model that fail earlier because of a
    bug and also resource problem
  - make pearson correlation matrix ggplot with the ensemble models in
    the corner

## 14 February 2020 (Friday)

Commits: 0

#### Worked on:

  - make the manually tuning for CNN-GRU and see if the models are
    overfitting (for the models done) for all:
      - bacteria
      - oomycete
      - fungi
      - all pathogen together
  - make a script with the k-fold working
  - make a report on how well each model performs

## 17 February 2020 (Monday)

Commits: 0

#### Worked on:

  - tuning manually the model with the extra epochs
  - for fungi oomycete: training manually with more epochs, since seems
    fungi is not good enough for only 30 epochs
  - read some papers about previous project and made a losts on what
    methods and tricks to use, so that I can addapt it as well.
  - making some Rmd with different pathogens, so that we can see the
    progress of each manually tuning of the models, also to get the data
    from previous work (Transfer Learning)
  - created a heatmap and simple graphs on python for Bacteria - CNN GRU

## 18 February 2020 (Tuesday)

Commits: 0

#### Weekly minutes:

  - Focus first on getting the manually tuning done so that all of the
    models are finished and polish
  - Try to make a story on the report, nice story, picture doesnt need
    to be with long description like in a paper
  - for the line chart results of the CNN-GRU, try to smooth them
    (fourier transform)

## 19–21 February 2020

Commits: 0

#### Worked on:

  - getting the smoother transformation using fft transformation
  - took care of some job failed (because of the time limit)
  - did analysis for some jobs that has finished
  - continue the scripts on tuning manually the data
  - ran the cnn-lstm Grad cam heatmap
  - moved all the results to local directories for fungi, all, and
    bacteria

## 24–28 February 2020

#### Activities:

  - went to Germany for an interview
