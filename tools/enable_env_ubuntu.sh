#!/bin/bash -l

set -e

export PHARML_ENV=pharml.bind
export PYTHONIOENCODING=utf8

source ~/anaconda3/etc/profile.d/conda.sh
conda activate $PHARML_ENV
