#!/usr/bin/env python
############################################################


import os
import sys
import argparse
import math
import warnings
import collections
import dask
import numpy as np
import pandas as pd
from dask.diagnostics import ProgressBar
from scipy import spatial
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools
from rdkit import RDLogger
from Bio.PDB import *
from pdb_to_nhg import convert_pdbs


############################################################


IC50_cutoff = 10000 # nM
BATCH_SIZE  = 512   # PDBs
max_lig_sz  = 200   # atoms
min_lig_sz  = 8     # atoms


############################################################


def write_lig_file(mol,ligfn):
    if os.path.exists(ligfn):
        return
    with open(ligfn,"w") as ligf:
        atoms = mol.GetAtoms();
        # Write number of atoms
        ligf.write("%d\n"%(len(atoms)))
        # Write atomic charges
        for atom in atoms:
            ligf.write("%s "%(str(atom.GetAtomicNum())))
        ligf.write("\n")
        # Write formal charges
        for atom in atoms:
            ligf.write("%s "%(str(atom.GetFormalCharge())))
        ligf.write("\n")
        # Write bond matrix
        for ndx_a in range(0, len(atoms)):
            for ndx_b in range(0, len(atoms)):
                bond = mol.GetBondBetweenAtoms(ndx_a, ndx_b)
                if bond != None:
                    if bond.GetIsAromatic():
                        ligf.write("4")
                    else:
                        btype = str(bond.GetBondType())
                        if btype == "SINGLE":
                            btype = "1"
                        elif btype == "DOUBLE":
                            btype = "2"
                        elif btype == "TRIPLE":
                            btype = "3"
                        ligf.write(str(btype))
                else:
                    ligf.write("0")                    
                ligf.write(" ")
            ligf.write("\n")
        # Write (empty) distance matrix
        for ndx_a in range(0, len(atoms)):
            for ndx_b in range(0, len(atoms)):
                ligf.write("0 ")
            ligf.write("\n")


def write_map_file(ligands,proteins,outdir="data/"):
    mapfn = outdir + "/map/dataset.map"
    bind_stats = {}
    with open(mapfn,"w") as mapf:
        out = ""
        for pdbid in proteins:
            n_bind   = 0
            n_nobind = 0
            for ligid in proteins[pdbid]:
                # Get the ligand and convert to bind / nobind.
                lig = ligands[ligid]
                ic50_str = str(lig[2])
                bind = True if ic50_str[0] != '>' and float(ic50_str[1:]) <= IC50_cutoff else False
                n_bind   += int(bind)
                n_nobind += int(not bind)
                # Write two input files, nhg before lig.
                out += "2 %s %s"%("../nhg/"+str(pdbid)+".nhg", "../lig/lig"+str(ligid)+".lig")
                # Write one output, the bind state as 1 / 0.
                bind_float = 1.0 if bind == True else 0.0
                out += " 1 %f"%(bind_float)
                # Go ahead and add a text tag for bind / nobind.
                bind_tag = "bind" if bind == True else "nobind"
                out += " 1 %s\n"%(bind_tag)
            # Flush the output buffer once per protein.
            mapf.write(out)
            out = ""
            bind_stats[pdbid] = (n_bind, n_nobind)
    return bind_stats


############################################################


def split_sdf(sdf_file_name,outdir="data/",lig_only=False,tag=""):
    print("Loading sdf.")
    # Parse the SDF file into a Pandas dataframe.
    rdk_lg = RDLogger.logger()
    rdk_lg.setLevel(RDLogger.CRITICAL)
    df = PandasTools.LoadSDF(sdf_file_name,
                             smilesName='SMILES',
                             molColName='Molecule',
                             includeFingerprints=False)
    #print("Raw cols = ", [str(x) for x in df.columns])
    # Take a simpler path when only extracting ligands.
    if lig_only:
        df_list=['SMILES','Molecule']
        df_selected = df[df_list].copy()
        # Drop any rows with missing data.
        df_selected = df_selected.replace('',  np.nan)
        df_selected = df_selected.replace(',', np.nan)
        df_selected = df_selected.dropna()
        r_rows = len(df.index)
        s_rows = len(df_selected.index)
        print("Raw rows = ", r_rows)
        print("Sel rows = ", s_rows)
        print("Keep pct = %.2f%s"%(((float(s_rows)/float(r_rows))*100.0),'%'))
        # Build ligand dictionary.
        uligs = {}
        for lndx,row in enumerate(df_selected.values):
            latoms = len(row[1].GetAtoms())
            if latoms >= min_lig_sz and latoms <= max_lig_sz:
                uligs[ tag+str(lndx) ] = row
        # Write out .lig files and return the data dictionaries.
        print("Writing per-ligand output files.")
        for key in uligs:
            write_lig_file(uligs[key][1],outdir+"/lig/lig%s.lig"%str(key))
        return uligs, None
    # Select only the needed columns and merge the two PDB cols.
    df_list=['PDB ID(s) for Ligand-Target Complex','PDB ID(s) of Target Chain','SMILES','IC50 (nM)','Molecule']
    df_selected = df[df_list].copy()
    df_selected["PDB IDs"] = df_selected['PDB ID(s) for Ligand-Target Complex'] + ',' + df_selected['PDB ID(s) of Target Chain']
    print("Selected cols = ", [str(x) for x in df_selected.columns])
    df_selected = df_selected[ ["PDB IDs"] + df_list[2:] ]
    # Drop any rows with missing data.
    df_selected = df_selected.replace('',  np.nan)
    df_selected = df_selected.replace(',', np.nan)
    df_selected = df_selected.dropna()
    r_rows = len(df.index)
    s_rows = len(df_selected.index)
    print("Raw rows = ", r_rows)
    print("Sel rows = ", s_rows)
    print("Keep pct = %.2f%s"%(((float(s_rows)/float(r_rows))*100.0),'%'))
    # Build ligand dictionary and a protein dictionary.
    print("Building protein-ligand dictionary.")
    uligs = {}
    prots_ligs = {}
    for lndx,row in enumerate(df_selected.values):
        latoms = len(row[-1].GetAtoms())
        if latoms < min_lig_sz or latoms > max_lig_sz:
            continue
        pdbs = row[0].split(',')
        for pdb in pdbs:
            if pdb == '':
                continue
            if pdb not in prots_ligs:
                prots_ligs[pdb] = []
            prots_ligs[pdb] += [ lndx ]
        uligs[ lndx ] = row
    print("Unique proteins = ", len(prots_ligs))
    print("Writing per-ligand output files.")
    # Write out .lig files and return the data dictionaries.
    for key in uligs:
        ndx = str(key)
        lig = uligs[key]
        write_lig_file(lig[3],outdir+"/lig/lig%s.lig"%ndx)
    return uligs, prots_ligs


############################################################


def main():
    # Parse command line args.
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf',     type=str,           required=True,   help='Path to SDF file.')
    parser.add_argument('--out',     type=str,           default="data",  help='Output directory.')
    parser.add_argument('--threads', type=int,           default=16,      help='Number of threads for PDB processing.')
    parser.add_argument('--ic50',    type=float,         default=10000.0, help='IC50 max cutoff for bind in nM.')
    parser.add_argument('--nhgr',    type=float,         default=4.0,     help='Local NHG radius in A.')
    parser.add_argument('--lig_only',action='store_true',default=False,   help='Only extract ligands from SDF.')
    args = parser.parse_args()
    IC50_cutoff = args.ic50
    # Create any needed subdirectories in advance.
    for subdir in ("pdb", "nhg", "lig", "map"):
        if not os.path.exists(args.out+"/"+subdir):
            os.makedirs(args.out+"/"+subdir)
    # First parse and process the SDF file for ligands and PDB-IDs.
    ligands, proteins = split_sdf(args.sdf,outdir=args.out,lig_only=args.lig_only)
    if args.lig_only:
        print("Success!")
        return
    # Next download and process the PDB files for the proteins.
    pdbids = [ str(key) for key in proteins ]
    rejected = convert_pdbs(pdbids,outdir=args.out,nworkers=args.threads,nhgr=args.nhgr)
    # Write out a map file listing all the resultant data items.
    print("Writing map file for the dataset:")
    rejected = { pdbid:ndx for ndx,pdbid in enumerate(rejected) }
    accepted = {}
    for pdbid in proteins:
        if pdbid not in rejected:
            accepted[pdbid] = proteins[pdbid]
    stats = write_map_file(ligands,accepted)
    # Get some stats / distributions.
    bind = 0
    nobind = 0
    for pdb in stats:
        bind += stats[pdb][0]
        nobind += stats[pdb][1]
    print("  Proteins:     %d"%(len(accepted)))
    print("  Ligands:      %d"%(len(ligands)))
    print("  Bind pairs:   %d"%(bind))
    print("  Nobind pairs: %d"%(nobind))
    dist = [ stats[pdb][0]+stats[pdb][1] for pdb in stats ]
    dist.sort()
    max_count = max(dist)
    dist = collections.Counter(dist)
    print("Writing distribution of ligand counts to 'pdb_dist.dat'.")
    with open("pdb_dist.dat","w") as distf:
        for count in range(1,max_count+1):
            num = dist[count] if count in dist else 0
            distf.write("%d %d\n"%(count,num))
    print("Success!")


if __name__ == "__main__":
    main()


############################################################
