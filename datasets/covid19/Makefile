.PHONY: all dataset dataset-bindingdb-fda dataset-bindingdb-zinc clean

all: dataset

dataset: dataset-bindingdb-fda

####

data/pdb/6vsb.cif:
	mkdir -p obsolete
	@echo "****** Get 6VSB from PDB and Convert to Neightborhood Graph ******"
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/pdb_to_nhg.py --pdbids 6VSB --out ./data"

####

../bindingdb-fda/data/lig:
	$(MAKE) -C ../bindingdb-fda/

dataset-bindingdb-fda: ../bindingdb-fda/data/lig data/pdb/6vsb.cif
	rm -rf ./data/lig
	ln -s ../../bindingdb-fda/data/lig ./data/lig
	@echo "****** Create a Map File to Inference from ../bindingdb-fda/ Ligands ******"
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/create_map.py --lig_dir ./data/lig --pdbs 6VSB --out ./data"

####

../zinc/data/lig:
	$(MAKE) -C ../zinc/

dataset-zinc: ../zinc/data/lig data/pdb/6vsb.cif
	rm -rf ./data/lig
	ln -s ../../zinc/data/lig ./data/lig
	@echo "****** Create a Map File to Inference from ../zinc/ Ligands ******"
	bash -c "source ../../tools/enable_env_ubuntu.sh ; python ../../tools/create_map.py --lig_dir ./data/lig --pdbs 6VSB --out ./data"

####

clean:
	rm -rf data
