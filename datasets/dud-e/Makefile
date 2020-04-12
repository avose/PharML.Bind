.PHONY: all dataset clean

CORES=$$(cat /proc/cpuinfo | grep 'processor' | grep ':' | wc -l)

all: dataset

dud-e:
	wget -c 'http://dude.docking.org/db/subsets/all/all.tar.gz'
	mv all.tar.gz dud-e.tar.gz
	tar -zxvf dud-e.tar.gz
	mv all dud-e
	rm -rf dud-e/fgfr1

data/map/dataset.map: dud-e
	mkdir -p obsolete
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/dud-e_to_dataset.py"

dataset: data/map/dataset.map

clean:
	rm -rf data
	rm -rf dud-e.tar.gz
	rm -rf dud-e
