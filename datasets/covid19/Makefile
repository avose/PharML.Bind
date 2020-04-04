.PHONY: all dataset clean

all: dataset

../bindingdb-fda/data/lig:
	$(MAKE) -C ../bindingdb-fda/

dataset: data/map/dataset.map ../bindingdb-fda/data/lig
	mkdir -p obsolete
	@echo "****** Get 6VSB.pdb and Convert to Neightborhood Graph ******"
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/pdb_to_nhg.py --pdbids 6VSB --out ./data"
	@echo "****** Create a Map File to Inference from ../bindingdb-fda/ Ligands ******"
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/create_map.py --lig_dir ../bindingdb-fda/data/lig --pdbs 6VSB --out ./data"

clean:
	rm -rf data
