#!/bin/bash -l

set -e

source ../../tools/enable_env_ubuntu.sh
wget 'http://www.bindingdb.org/bind/downloads/BindingDB_All_terse_2D_2019m4.sdf.zip'
unzip BindingDB_All_terse_2D_2019m4.sdf.zip

CORES=`cat /proc/cpuinfo | grep 'processor' | grep ':' | wc -l`
python sdf_to_dataset.py --sdf ./BindingDB_All_terse_2D.sdf --out ./data --threads ${CORES} --ic50 10000
