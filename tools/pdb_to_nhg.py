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


############################################################


BATCH_SIZE  = 512   # PDBs
max_ptn_sz  = 10000 # atoms
min_ptn_sz  = 500   # atoms


############################################################


NAME_TO_CHARGE = {
    'H':1,
    'D':1,   #Deuterium?
    'C':6, 
    'N':7, 
    'O':8, 
    'F':9, 
    'G':12,  #Mg
    'P':15, 
    'S':16, 
    'L':17,  #Cl
    'Z':30,  #Zn
    'E':34,  #Se
    'R':35   #Br
}


def atomic_charge(aname):
    # Map first letter of atom name from PDB to an atomic charge.
    atype = aname[0]
    if atype not in NAME_TO_CHARGE:
        print("Error: Unknown atom type: '%s' from '%s'!"%(atype,aname))
        sys.exit()
    return NAME_TO_CHARGE[atype]

    
def write_nhg_file(atoms,edges,nhgfn):
    with open(nhgfn,"w") as nhgf:
        # Write num atoms and num edges
        out = [ len(edges), len(atoms) ]
        out = np.array(out, dtype=np.int32)
        out.tofile(nhgf)
        # Write each atom as a 5-tuple: (t,b,x,y,z)
        out = []
        for atom in atoms:
            apos = atom.get_vector()
            atype = atomic_charge(atom.get_name())
            out += [ atype, 0.0, apos[0], apos[1], apos[2] ]
        out = np.array(out, dtype=np.float32)
        out.tofile(nhgf)
        # Write each edge as a 3-tuple: (d,ndx_i,ndx_j)
        out = []
        for edge in edges:
            out += [ edge[2], edge[0], edge[1] ]
        out = np.array(out, dtype=np.float32)
        out.tofile(nhgf)


############################################################


def convert_pdb(pdbid,outdir="data/",nhgr=4):
    # Print header for each protein and download if needed.
    status = "---------------- %s ----------------\n"%pdbid
    bp_pdbl = PDBList(verbose=False)
    bp_parser = MMCIFParser()
    bp_pdbl.retrieve_pdb_file(pdbid,pdir=outdir+"/pdb/",file_format="mmCif")
    pdb_fn = outdir+"/pdb/"+pdbid.lower()+".cif"
    # Parse the PDB / CIF file into a structure object.
    structure = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            structure = bp_parser.get_structure("",pdb_fn)
        except Exception as ex:
            status += "!! Failed to parse '%s', skipping 'invalid' PDB ('%s').\n"%(pdb_fn,str(ex))
            return (pdbid, False, status)
    # Traverse the PDB structure, looking for the protein chain.
    status += "Models: %d\n"%(len(structure))
    chain_A = []
    for mndx, model in enumerate(structure):
        status += "  models[%d]: %d chains\n"%(mndx,len(model))
        chain = None
        if 'A' in model:
            chain = model['A']
            status += "    Selected chain: 'A'\n"
        else:
            for c in model:
                chain = c
                break
            status += "    Selected chain: first ('%s')\n"%(chain.get_id())
        atoms = []
        for residue in chain:
            resid = residue.get_full_id()
            if resid[3][0] == "W" or (resid[3][0][0] == "H" and resid[3][0][1] == "_"):
                continue                    
            for atom in residue:
                atoms.append( atom )
        chain_A = atoms
        # !!av: For now, just consider the first model.
        break
    status += "    Selected chain len: %d\n"%(len(chain_A))
    # Convert the selected chain into a NHG.
    status += "Building neighborhood graph.\n"   
    try:
        # Get an array of the atom positions.
        positions = []
        for atom in chain_A:
            positions += list(atom.get_vector()) 
        positions = np.array(positions, dtype=np.float32)
        positions = np.reshape(positions, (len(chain_A),3))
        if positions.shape[0] > max_ptn_sz:
            raise Exception('Protein too large (%d)!'%(positions.shape[0]))
        if positions.shape[0] < min_ptn_sz:
            raise Exception('Protein too small (%d)!'%(positions.shape[0]))
        # Find pairwise distances between all atoms.
        distances = spatial.distance.pdist(positions)
        distances = spatial.distance.squareform(distances)
        # Find local neighborhoods by removing long (and self) distances.
        distances[distances == 0.0] = 1.0 + nhgr
        neighborhoods = np.nonzero(distances <= nhgr)
        # Turn the distances / indicies into an explicit edge list.
        nh_edges = [ (ndx_a, ndx_b, distances[ndx_a][ndx_b]) for ndx_a, ndx_b in zip(*neighborhoods) ]
        # Examine the local neighborhood sizes to check connectivity.
        local_nh_sizes = { ndx:0 for ndx in np.unique(neighborhoods[0]) }
        for ndx in neighborhoods[0]:
            local_nh_sizes[ndx] += 1
        if min(local_nh_sizes.values()) < 3:
            raise Exception('Too few local edges (%d)!'%(min(local_nh_sizes.values())))
        status += "Neighborhood graph info:\n"
        status += "  Total edges: %d"%(len(nh_edges))
        status += "  Local edges:\n"
        status += "    Min edges: %d\n"%(min(local_nh_sizes.values()))
        status += "    Max edges: %d\n"%(max(local_nh_sizes.values()))
    except Exception as ex:
        # Return a failure if anything went wrong.
        status += "!! Neighborhood graph failed for '%s', skipping PDB ('%s').\n"%(pdb_fn,str(ex))
        return (pdbid, False, status)
    # Write the output NHG and return success.
    status += "Writing neighborhood graph for '%s'.\n"%(pdbid)
    write_nhg_file(chain_A,nh_edges,outdir+"/nhg/%s.nhg"%(pdbid))
    status += "Done converting PDB.\n"
    return (pdbid, True, status)


def chunks(lst, size):
    # Turn the list into cunks of size 'size'.
    for i in range(0, len(lst), size):
        yield lst[i:i+size]


def convert_pdbs(pdbids,outdir="data/",nworkers=16,nhgr=4):
    print("Processing PDBs.")
    # Maintain a list of rejected PDBs for later.
    rejected_pdbs = []
    # Split all the PDBs into batches for Dask.
    batches = [ batch for batch in chunks(pdbids, BATCH_SIZE) ]
    for bndx, batch in enumerate(batches):
        print("PDB batch %d of %d (%.2f%s):"%(bndx+1,len(batches),float(bndx)/len(batches)*100.0,"%"))
        # Create lsit of Dask tasks and launch them with threads.
        results = [ dask.delayed(convert_pdb)(pdbid,outdir,nhgr) for pdbid in batch ]
        with ProgressBar(dt=0.5):
            results = dask.compute(*results, scheduler='threads', num_workers=nworkers)
        # Process the results from the batch.
        for pdbid, flag, status in results:
            print("%s"%(status),end="")
            if not flag:
                rejected_pdbs.append(pdbid)
    # Finished, so print stats and return the rejected list.
    print("--------------------------------------")
    print("Done converting PDBs:")
    print("  Rejected: %d/%d"%(len(rejected_pdbs),len(pdbids)))
    print("  Accepted: %.2f%s"%(float(len(pdbids)-len(rejected_pdbs))*100.0/float(len(pdbids)),'%'))
    return rejected_pdbs


############################################################


def main():
    # Parse command line args.
    parser = argparse.ArgumentParser()
    parser.add_argument('--out',     type=str,   default="data",  help='Output directory.')
    parser.add_argument('--threads', type=int,   default=16,      help='Number of threads for PDB processing.')
    parser.add_argument('--nhgr',    type=float, default=4.0,     help='Local NHG radius in A.')
    parser.add_argument('--pdbids',  type=str,   required=True,   help='List of PDB IDs to convert to NHGs.')
    args = parser.parse_args()
    # Create any needed subdirectories in advance.
    for subdir in ("pdb", "nhg"):
        if not os.path.exists(args.out+"/"+subdir):
            os.makedirs(args.out+"/"+subdir)
    # Next download and process the PDB files for the proteins.
    pdbids = args.pdbids.split(",")
    rejected = convert_pdbs(pdbids,outdir=args.out,nworkers=args.threads,nhgr=args.nhgr)
    rejected = { pdbid:ndx for ndx,pdbid in enumerate(rejected) }
    accepted = 0
    for pdbid in pdbids:
        if pdbid not in rejected:
            accepted += 1
    print("  Proteins: %d"%(accepted))
    print("Success!")


if __name__ == "__main__":
    main()


############################################################
