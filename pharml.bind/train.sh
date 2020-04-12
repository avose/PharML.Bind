#!/bin/bash -l

source ../tools/enable_env_ubuntu.sh

MAP_NAME=new_train
MAP_TRAIN_PATH=../datasets/dud-e/data/map/dataset.map
MAP_TEST_PATH=../datasets/dud-e/data/map/dataset.map
MLP_LATENT=32,32
MLP_LAYERS=2,2
GNN_LAYERS=5,5
NUM_FEATURES=16,16
MODE=classification
TEMP_DIR=${MAP_NAME}
BS=2
NT=2
EPOCHS=100
LR0=0.000000001

echo "============================================================"
echo "Training: name \"${MAP_NAME}\", train ${MAP_TRAIN_PATH}, test ${MAP_TEST_PATH}"
echo "============================================================"

rm -rf ${TEMP_DIR}
mkdir -p ${TEMP_DIR}
cd ${TEMP_DIR}
# Start the training run on a single model.
python ../mldock_gnn.py \
       --map_train "../${MAP_TRAIN_PATH}" \
       --map_test  "../${MAP_TEST_PATH}" \
       --batch_size ${BS} \
       --batch_size_test ${BS} \
       --mlp_latent ${MLP_LATENT} \
       --mlp_layers ${MLP_LAYERS} \
       --gnn_layers ${GNN_LAYERS} \
       --lr_init ${LR0} \
       --use_clr True \
       --hvd False \
       --num_features ${NUM_FEATURES} \
       --data_threads ${NT} \
       --mode ${MODE} \
       --epochs ${EPOCHS} 2>&1 |& tee train.out
cd ..
echo "============================================================"
echo "Training finished."
echo "============================================================"
