#!/bin/bash -l

set -e

export PHARML_ENV=pharml.bind
export PYTHONIOENCODING=utf8

source activate $PHARML_ENV
