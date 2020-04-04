#!/usr/bin/env python
############################################################


# http://zinc15.docking.org/substances.csv?count=all


############################################################


import os
import sys
import argparse
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools
from rdkit import RDLogger
sys.path.append("../../tools/")
from sdf_to_dataset import write_lig_file


############################################################


def load_csv(csv_file_name):
    print("Loading CSV.")
    # Parse the CSV file.
    rdk_lg = RDLogger.logger()
    rdk_lg.setLevel(RDLogger.CRITICAL)
    with open(csv_file_name,"r") as csvf:
        ligands = [ list(line.split(",")) for line in csvf.read().split("\n") ]
    ligands = ligands[1:]
    # Convert to mol objects.
    print("Converting ligands to mol objects.")
    valid_ligands = []
    for ndx,ligand in enumerate(ligands):
        if len(ligand) == 1 and ligand[0] == "":
            continue
        if len(ligand) != 2:
            print(ligand)
            continue
        ligand.append(Chem.MolFromSmiles(ligand[1]))
        valid_ligands.append(ligand)
        if ndx < 10:
            print(ligand)
        elif ndx == 10:
            print("...")
    print("Done creating mol objects.")
    return valid_ligands


def write_ligands(ligands,outdir="data/"):
    print("Writing %d .lig files."%(len(ligands)))
    for ligand in ligands:
        write_lig_file(ligand[2],outdir+"/lig/"+ligand[0]+".lig")
    print("Done writing %d files."%(len(ligands)))


############################################################


if __name__ == "__main__":
    # Parse command line args.
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv',     type=str,   required=True,   help='Path to CSV file.')
    parser.add_argument('--out',     type=str,   default="data",  help='Output directory.')
    args = parser.parse_args()
    # Create any needed subdirectories in advance.
    for subdir in ("lig",):
        if not os.path.exists(args.out+"/"+subdir):
            os.makedirs(args.out+"/"+subdir)
    # Read input CSV file.
    ligands = load_csv(args.csv)
    # Write out as .lig files.
    write_ligands(ligands,outdir=args.out)
    # Done.
    print("Success!")
    
    
############################################################
