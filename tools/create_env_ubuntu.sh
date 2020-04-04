#!/bin/bash -l

set -e

export PHARML_ENV=pharml.bind

conda create -y -n $PHARML_ENV python=3.6 cudatoolkit=10.0 cudnn
source activate $PHARML_ENV
conda install -y -c conda-forge rdkit biopython scipy dask
conda install -y -c anaconda pip
pip install tensorflow_gpu==1.15
pip install graph_nets dm-sonnet==1.25 'tensorflow-probability<0.9.0' matplotlib
pip install horovod
source deactivate
conda env list
