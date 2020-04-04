#!/bin/bash

source ../tools/enable_env_ubuntu.sh

MAP_NAME=6vsb-fda
MAP_PATH=../datasets/covid19/data/map/dataset.map
MLP_LATENT=32,32
MLP_LAYERS=2,2
GNN_LAYERS=5,5
NUM_FEATURES=16,16
MODE=classification
EPOCHS=1
INFER_OUT=inference_output.map
TEMP_DIR=./pharml-covid19
BS=8
NT=2

for model in `ls -d ../pretrained-models/mh-gnnx5-ensemble/model_*` ; do
    echo "============================================================"
    echo "Running inference with name \"${MAP_NAME}\" using ${MAP_PATH}"
    echo "Current ensemble model: ${model}"
    echo "============================================================"
    rm -rf $TEMP_DIR
    mkdir -p ${TEMP_DIR}
    cd ${TEMP_DIR}
    # Start the inference run on a single model
    python ../mldock_gnn.py \
        --map_train "../${MAP_PATH}" \
        --map_test  "../${MAP_PATH}" \
        --batch_size ${BS} \
        --batch_size_test ${BS} \
        --mlp_latent ${MLP_LATENT} \
        --mlp_layers ${MLP_LAYERS} \
        --gnn_layers ${GNN_LAYERS} \
        --hvd False \
        --num_features ${NUM_FEATURES} \
        --data_threads ${NT} \
        --mode ${MODE} \
        --inference_only True \
        --restore="../${model}/checkpoints/model0.ckpt" \
        --inference_out ${INFER_OUT} \
        --epochs 1
    cd ..
    echo "============================================================"
    echo "Inference finished with name \"${MAP_NAME}\" using ${MAP_PATH}"
    echo "Current ensemble model: ${model}"
    echo "============================================================"
done
echo "============================================================"
echo "Done testing COVID-19 compounds against all models"
echo "============================================================"
wait
