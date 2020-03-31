
# November 2019 Lab Report <img src="figures/tsl-logo.png" align="right" width="120" />

### [Ruth Kristianingsih](https://github.com/ruthkr)

## 01 November 2019 (Friday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Finished the confusion matrix.
  - Checked the jobs of subcellular multicompartment.
  - Make report for both.
  - Finishing lab books September, October.

-----

## 04 November 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Weekly catch up minutes

##### Sub-cellular protein project

  - The better model so far is CNN with 16 filter CNN and filter LSTM
    since it has more generality than with filter 64 by 64 or 48 by 48
    (with acc= 40 % better than random) ==\> suggests that CNN did not
    give any good representation of the model, so then we need to change
    the network ( for sequence maybe CNN only is better).

##### Effector-non effector project

  - From conffusion matrix, it suggests that the model is conservative,
    it tends to say no if they can not recognise the sequence, the
    difference between false negative and false positive is around 1/4.
  - For next step:
      - Make the heatmap from the layer of CNN after activation to see
        effector motives.

## 05 November 2019 (Tuesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Did the Visa documents scanning and uploading.
  - Tried to understand how to visualize the result of CNN filter on the
    matrix one hot encoding data:
      - Got an error `ImportError: Could not import PIL.Image. The use
        of`array\_to\_img\` requires PIL.
      - Solved error by installing `conda install pillow`.

## 06 November 2019 (Wednesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Going to talks in Autumn Symposium.

## 07 November 2019 (Thursday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Going to Croydon for the Visa application.

## 08 November 2019 (Friday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Continue with heatmap for the effector and non-effector project.
  - Thinking about the next step of the subcellular localisation (after
    the diagram that Dan showed).

-----

## 11 November 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Finishing the heatmap.
  - Read the GAN for generating sequences.

#### Weekly catch up minutes

##### Sub-cellular protein project

  - I proposed an idea to start with GAN and data augmentation for this
    project, and Dan agreed.

##### Effector-non effector project

  - I did wrong way to get the heatmap for the sequence data (looked at
    the book of Chollet).

## 12–13 November 2019

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Figure out the formula for Grad-Cam and how to implement it on
    sequence data.
  - Did the plot and tried to combine for both the original plot of
    sequence data and heatmap.
  - Install library `OpenCV` on tensorflow env, but it downgrade the
    hdf5, and the pretarined model does not want to load properly.
  - Uninstall OpenCV: `conda remove -n tensorflow opencv`.
  - `conda update -n tensorflow h5py hdf5`.

## 14 November 2019 (Thursday) – Annual Leave

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

  - Took a day off from the Annual Leave (moving to a new place).

## 15 November 2019 (Friday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Tried to fix the error on `tensorflow` conda env.

-----

## 18 November 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Making the heatmap using the OpenCV.
  - Could not get the result using the work PC, therefore using the
    personal PC to get the result.

## 19 November 2019 (Tuesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Weekly catch up meeting

##### Sub-cellular protein project

  - For the Discriminator network, Dan said it is okay to use the
    current models we have (CNN-LSTM).

##### Effector-non effector project

  - Making all of the of the combined from all of matrices for each
    datasets (training, validation, and testing).
  - Soon to wrap it up become a paper that we can publish.
  - What if we combine with the model that already exist? next week:
    purpose some models (published model).

## 20 November 2019 (Wednesday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Making the matrix heatmap calculation.
  - Could not run in the cluster (`Error: there is no keras/tensorflow
    module`).

## 21–22 November 2019

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Fixed the bugs from the matrix heatmap calculations (takes a lot of
    memory).
  - Try to create the GAN architecture.
  - Getting the value of 8 physicochemical attributes (Kawashima, Ogata
    and Kanehisa, 1998).

-----

## 25 November 2019 (Monday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Weekly catch up meeting

##### Sub-cellular protein project

  - Can be postponed first and prioritize the effector and non-effector
    one as the result is more obvious.

##### Effector and non-effector prediction

  - There are two things to do in order to see how the model behave
    specifically:
      - Shrink the datasets into different genus and run it spearately.
      - Or tweak the model into multiclass classification for different
        genus.
  - Ensemble the current model with any model, can be deep learning or
    can be other models, as long as the models resulting in the same
    classifocation with the model that we have now (and also as long as
    it classify all kind of effector).

## 26–28 November 2019

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Get the data grouping for both cases.
  - Encode them.
  - Had problem with my conda env.

## 29 November 2019 (Friday)

<!-- Commits: `#r gitlogr::get_git_commit_count(curr_date)` -->

#### Worked on

  - Train all of the datasets in 2 different models (CNN-LSTM, CNN-GRU).
