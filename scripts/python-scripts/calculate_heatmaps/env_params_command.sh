# Single run example
# CLASS_NAME="bacteria" DATA_NAME="x_test" MODEL_PATH="models/bacteria/cnn-lstm/model_1.30-0.41.hdf5" LAYER="conv1d_4" FROM_SAMPLE=0 SAMPLE_LENGTH=25 jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=None --execute heatmap_calculate_in_parts_splits_params.ipynb

export CLASS_NAME_RAW="fungi"
export LAST_SAMPLE_RAW=38
export SAMPLE_LEGTH_RAW=20
export DATA_NAME_RAW="x_test"
export MODEL_PATH_RAW="../../../data/secreted_data/saved_models/fungi/cnn_gru/sequential_1.60-0.43.hdf5"
export LAYER_RAW="conv1d_1"

for sample in $(seq 0 $SAMPLE_LEGTH_RAW $LAST_SAMPLE_RAW); do
	CLASS_NAME=$CLASS_NAME_RAW DATA_NAME=$DATA_NAME_RAW MODEL_PATH=$MODEL_PATH_RAW LAYER=$LAYER_RAW FROM_SAMPLE=$sample SAMPLE_LENGTH=$SAMPLE_LEGTH_RAW jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=None --execute heatmap_calculate_in_parts_splits_params.ipynb
done

osascript -e 'display notification "Finished x_test calculations!" with title "JupyterLab" sound name "Hero"'

export CLASS_NAME_RAW="fungi"
export LAST_SAMPLE_RAW=38
export SAMPLE_LEGTH_RAW=20
export DATA_NAME_RAW="x_val"
export MODEL_PATH_RAW="../../../data/secreted_data/saved_models/fungi/cnn_gru/sequential_1.60-0.43.hdf5"
export LAYER_RAW="conv1d_1"

for sample in $(seq 0 $SAMPLE_LEGTH_RAW $LAST_SAMPLE_RAW); do
	echo $sample
	CLASS_NAME=$CLASS_NAME_RAW DATA_NAME=$DATA_NAME_RAW MODEL_PATH=$MODEL_PATH_RAW LAYER=$LAYER_RAW FROM_SAMPLE=$sample SAMPLE_LENGTH=$SAMPLE_LEGTH_RAW jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=None --execute heatmap_calculate_in_parts_splits_params.ipynb
done

osascript -e 'display notification "Finished x_val calculations!" with title "JupyterLab" sound name "Hero"'

export CLASS_NAME_RAW="fungi"
export LAST_SAMPLE_RAW=118
export SAMPLE_LEGTH_RAW=25
export DATA_NAME_RAW="x_train"
export MODEL_PATH_RAW="../../../data/secreted_data/saved_models/fungi/cnn_gru/sequential_1.60-0.43.hdf5"
export LAYER_RAW="conv1d_1"

for sample in $(seq 0 $SAMPLE_LEGTH_RAW $LAST_SAMPLE_RAW); do
	echo $sample
	CLASS_NAME=$CLASS_NAME_RAW DATA_NAME=$DATA_NAME_RAW MODEL_PATH=$MODEL_PATH_RAW LAYER=$LAYER_RAW FROM_SAMPLE=$sample SAMPLE_LENGTH=$SAMPLE_LEGTH_RAW jupyter nbconvert --to notebook --inplace --ExecutePreprocessor.timeout=None --execute heatmap_calculate_in_parts_splits_params.ipynb
done

osascript -e 'display notification "Finished x_train calculations!" with title "JupyterLab" sound name "Hero"'
