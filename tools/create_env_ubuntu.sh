#!/bin/bash -l

set -e

export PHARML_ENV=pharml.bind

conda create -n $PHARML_ENV -c conda-forge python=3.6 cudatoolkit=10.0 cudnn biopython scipy dask pip rdkit
conda activate $PHARML_ENV
pip install tensorflow_gpu==1.15
pip install graph_nets dm-sonnet==1.25 'tensorflow-probability<0.9.0' matplotlib
pip install horovod
source deactivate
conda env list
