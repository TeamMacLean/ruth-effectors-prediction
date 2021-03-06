
# September 2019 Lab Report <img src="figures/tsl-logo.png" align="right" width="120" />

### [Ruth Kristianingsih](https://github.com/ruthkr)

## 01–20 September 2019

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Built a skeleton for Report.
  - Started to make Rmd template and started to make the report.

##### Sub-cellular localisation

  - Getting subcellular multicompartment data from stractch (including
    read some reference papers).
  - Process:
      - From the TAIR data, get all of protein that has location and
        with the charateristics EDA.
      - Learn about the data and see if we can already group them
        together.
  - Issues:
      - Same ID proteins has more than one location, we needed to
        analyze there.
      - There are more than 150 locations, therefore it is not possible
        to have 150 labels.
  - Actions:
      - Made analysis using histogram and see how the locations of
        protein disributes.
      - Take only protein with one location.
      - Using reference to determine the location that is going to be
        used. The locations are: Chloroplast, Cytoplasm, Nucleus,
        Peroxisome, Plasma. Membrane, Vacuole, ER, Extracellular, Golgi
        apparatus, and mitochondria (labeled as reference).
      - There are three different cases in how we grouped the protein.
        (See the reports).

##### Effector and Non-effector prediction

  - Constructed four different models for effector and non-effector
    project:
      - CNN - LSTM.
      - CNN - GRU.
      - LSTM - Embedding.
      - GRU - Embedding.

When we use embedding here, we use different kind of encoding, which are
integer encoding since using embedding, we represented each amino acid
into vector.

-----

## 23 September 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Finished report before meeting.
  - Contacted Dr. Aylin from EI.

#### Weekly catch up minutes

##### Sub-cellular localisation

  - Points:
      - Ask Martin Aylin from EI to use the EI GPUs (make a draft email
        and send it to Dan first to check).
      - start the first model for subcellular compartment localisation
        (do not start from the fully connected dense, start from
        CNN-LSTM).

##### Effector-noneffector

  - Points:
      - Next time, if you know there is a service disruption, you can
        pause job or do the model check point.
      - From the results from CNN-LSTM, train manually and predict the
        test data.
      - Do grid search from other models.

## 24–27 September 2019

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Tried to get the most optimized paramaters to scan using Scikit
    Learn RandomSearchCV and GridSearchCV.

  - Created parallel scripts from each models for RandomSearchCV and
    GridSearchCV.
    
      - CNN-LSTM:
          - `0001-cnn-lstm-scan-batch16-unit16.py`.
          - `0002-cnn-lstm-scan-batch16-unit32.py`.
          - `0003-cnn-lstm-scan-batch32-unit16.py`.
          - `0004-cnn-lstm-scan-batch32-unit32.py`.
      - CNN-GRU:
          - `0005-cnn-gru-grid-batch8_gru_unit16.py`.
          - `0006-cnn-gru-grid-batch16_gru_unit32.py`.
          - `0007-cnn-gru-grid-batch8_gru_unit64.py`.
          - `0008-cnn-gru-grid-batch16_gru_unit16.py`.
      - GRU-Embedding:
          - `0010-gru-embedding_scan.py`.
      - LSTM-Embedding:
          - `0013-lstm-embedding-grid-batch8_lstm_unit16.py`.
          - `0014-lstm-embedding-grid-batch16_lstm_unit16.py`.
          - `0015-lstm-embedding-grid-batch8_lstm_unit32.py`.
          - `0016-lstm-embedding-grid-batch16_lstm_unit32.py`.

  - Read about how to get accuracy for each class in multi-class
    compartment localisation.

-----

## 30 September 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Weekly catch up minutes

##### Sub-cellular localisation

  - Points:
      - Make the scripts of subcellular multicompartment parallel.

##### Effector-noneffector

  - Points:
      - Stop the CNN-LSTM and CNN-GRU since no need to do again
        GridSearchCV.
      - Do not forget to check the memory of each job.
      - Check memory of each job.
