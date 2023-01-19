.PHONY: all dataset clean

CORES=$$(cat /proc/cpuinfo | grep 'processor' | grep ':' | wc -l)

all: dataset

BindingDB_All_terse_3D.sdf:
	wget 'https://www.bindingdb.org/bind/downloads/BindingDB_All_terse_3D_2023m0.sdf.zip'
	unzip BindingDB_All_terse_3D_2023m0.sdf.zip

dataset: BindingDB_All_terse_3D.sdf
	mkdir -p obsolete
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/sdf_to_dataset.py --sdf ./BindingDB_All_terse_3D.sdf --out ./data --threads ${CORES} --ic50 10000"

clean:
	rm -rf data
	rm -rf zinc15.csv
