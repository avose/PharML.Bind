#!/bin/bash -l

set -e

source ../../tools/enable_env_ubuntu.sh
mkdir -p ./data/pdb
wget 'https://files.rcsb.org/download/6VSB.pdb'
mv 6VSB.pdb ./data/pdb

python ../../tools/create_map.py --lig_dir ../zinc/data/lig --pdbs ./data/pdb --out ./data
