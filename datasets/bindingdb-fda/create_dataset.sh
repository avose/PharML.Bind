#!/bin/bash -l

set -e

source ../../tools/enable_env_ubuntu.sh
wget 'http://www.bindingdb.org/bind/drugatfda_BindingDBf2D.sdf'

CORES=`cat /proc/cpuinfo | grep 'processor' | grep ':' | wc -l`
mkdir -p obsolete
python ../../tools/sdf_to_dataset.py --sdf ./drugatfda_BindingDBf2D.sdf --out ./data --threads ${CORES} --lig_only
