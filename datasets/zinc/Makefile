.PHONY: all dataset clean

all: dataset

zinc15.csv:
	curl -o zinc15.csv 'http://zinc15.docking.org/substances.csv?count=all'

dataset: zinc15.csv
	mkdir -p obsolete
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python csv_to_dataset.py --csv zinc15.csv --out ./data"

clean:
	rm -rf data
	rm -rf zinc15.csv
