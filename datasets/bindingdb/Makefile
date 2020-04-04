.PHONY: all dataset clean

CORES=$$(cat /proc/cpuinfo | grep 'processor' | grep ':' | wc -l)

all: dataset

BindingDB_All_terse_2D.sdf:
	wget 'http://www.bindingdb.org/bind/downloads/BindingDB_All_terse_2D_2019m4.sdf.zip'
	unzip BindingDB_All_terse_2D_2019m4.sdf.zip

dataset: BindingDB_All_terse_2D.sdf
	mkdir -p obsolete
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/sdf_to_dataset.py --sdf ./BindingDB_All_terse_2D.sdf --out ./data --threads ${CORES} --ic50 10000"

clean:
	rm -rf data
	rm -rf zinc15.csv
