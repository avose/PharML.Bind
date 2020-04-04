.PHONY: all dataset clean

CORES=$$(cat /proc/cpuinfo | grep 'processor' | grep ':' | wc -l)

all: dataset

drugatfda_BindingDBf2D.sdf:
	wget 'http://www.bindingdb.org/bind/drugatfda_BindingDBf2D.sdf'


dataset: drugatfda_BindingDBf2D.sdf
	mkdir -p obsolete
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/sdf_to_dataset.py --sdf ./drugatfda_BindingDBf2D.sdf --out ./data --threads ${CORES} --lig_only"

clean:
	rm -rf data
	rm -rf drugatfda_BindingDBf2D.sdf
