#!/bin/bash -l

set -e

source ../../tools/enable_env_ubuntu.sh
curl -o zinc15.csv 'http://zinc15.docking.org/substances.csv?count=all'

python csv_to_dataset.py --csv zinc15.csv --out ./data
