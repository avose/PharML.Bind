#!/bin/bash -l

set -e

source ../../tools/enable_env_ubuntu.sh

mkdir -p obsolete
echo "****** Get 6VSB.pdb and Convert to Neightborhood Graph ******"
python ../../tools/pdb_to_nhg.py --pdbids "6VSB" --out ./data
echo "****** Create a Map File to Inference from ../bindingdb-fda/ Ligands ******"
python ../../tools/create_map.py --lig_dir ../bindingdb-fda/data/lig --pdbs "6VSB" --out ./data
