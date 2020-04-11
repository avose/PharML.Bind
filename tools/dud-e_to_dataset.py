#!/usr/bin/env python
############################################################


import os
import sys
import argparse
import math
import shutil
import warnings
import collections
import numpy as np
import pandas as pd
from scipy import spatial
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools
from rdkit import RDLogger
from Bio.PDB import *
from pdb_to_nhg import convert_pdbs
from sdf_to_dataset import split_sdf


############################################################

max_lig_sz  = 200   # atoms
min_lig_sz  = 8     # atoms

############################################################

  
def write_map_file(items,out="data/"):
    mapfn = out + "/map/dataset.map"
    with open(mapfn,"w") as mapf:
        out = ""
        for item in items:
            pdbid, ligid, bind = item
            # Write two input files, nhg before lig.
            out += "2 %s %s"%("../nhg/"+str(pdbid)+".nhg", "../lig/lig"+str(ligid)+".lig")
            # Write one output, the bind state as 1 / 0.
            bind_float = 1.0 if bind else 0.0
            out += " 1 %f"%(bind_float)
            # Go ahead and add a text tag for bind / nobind.
            bind_tag = "bind" if bind else "nobind"
            out += " 1 %s\n"%(bind_tag)
        # Write the final map file.
        mapf.write(out)
        return


############################################################


def main():
    # Parse command line args.
    parser = argparse.ArgumentParser()
    parser.add_argument('--out',     type=str, default="data",  help='Output directory.')
    parser.add_argument('--inp',     type=str, default="dud-e", help='Input directory.')
    parser.add_argument('--threads', type=int, default=2,       help='Number of threads for PDB processing.')
    args = parser.parse_args()
    # Create any needed subdirectories in advance.
    for subdir in ("pdb", "nhg", "lig", "map"):
        if not os.path.exists(args.out+"/"+subdir):
            os.makedirs(args.out+"/"+subdir)
    # Copy PDBs into data dir and convert to nhg.
    pdbids = []
    for entry in os.listdir(args.inp):
        pdbids.append(entry)
        pdbf = args.inp + "/"     + entry + "/receptor.pdb"
        dest = args.out + "/pdb/" + entry + ".pdb"
        print("%s -> %s"%(pdbf,dest))
        shutil.copyfile(pdbf,dest)
    rejected = convert_pdbs(pdbids, outdir=args.out, nworkers=args.threads, download=False)
    rejected = { pdbid:ndx for ndx,pdbid in enumerate(rejected) }
    accepted = {}
    for pdbid in pdbids:
        if pdbid not in rejected:
            accepted[pdbid] = pdbid
    # Parse all the SDF ligand files.
    binds   = []
    nobinds = []
    for pdbid in pdbids:
        bind_sdf   = args.inp + "/" + pdbid + "/actives_final.sdf"
        nobind_sdf = args.inp + "/" + pdbid + "/decoys_final.sdf"
        print(bind_sdf,nobind_sdf)
        if not os.path.exists(bind_sdf) and os.path.exists(bind_sdf+".gz"):
            os.system("gunzip "+bind_sdf+".gz")
        if not os.path.exists(nobind_sdf) and os.path.exists(nobind_sdf+".gz"):
            os.system("gunzip "+nobind_sdf+".gz")
        bind, _   = split_sdf(bind_sdf,   outdir=args.out, lig_only=True, tag=pdbid+"_1_")
        nobind, _ = split_sdf(nobind_sdf, outdir=args.out, lig_only=True, tag=pdbid+"_0_")
        for key in bind:
            smiles, mol = bind[key]
            binds.append( [pdbid, key, True] )
        for key in nobind:
            smiles, mol = nobind[key]
            nobinds.append( [pdbid, key, False] )
    # Write out the final map file.
    prots_ligs = {}
    for key in bind:
        prots_ligs[key] = bind[key]
    for key in nobind:
        prots_ligs[key] = nobind[key]
    print("Writing map file for the dataset:")
    write_map_file(prots_ligs,out=args.out)
    print("  Proteins:     %d"%(len(accepted)))
    print("  Ligands:      %d"%(len(prots_ligs)))
    print("  Bind pairs:   %d"%(len(binds)))
    print("  Nobind pairs: %d"%(len(nobinds)))
    #stats = write_map_file(ligands,proteins)
    print("Success!")


if __name__ == "__main__":
    main()


############################################################
