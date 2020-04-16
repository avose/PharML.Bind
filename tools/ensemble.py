#!/usr/bin/env python
############################################################


import os
import sys
import time
import argparse
import numpy as np

from chemio import read_map 


############################################################


def ensemble_maps(maps):
    # Do a sanity-check on item0.
    for m in maps:
        first_item = m[0] #[inputs, outputs, tags]
        if len(first_item[0]) != 2:
            print("Expected two inputs (NHG and LIG)!")
            sys.exit(1)
        if len(first_item[1]) != 1:
            print("Expected one output (bind/nobind)!")
            sys.exit(1)
        if len(first_item[2]) != 4:
            print("Expected four tags (actual, prediction, pred0, pred1)!")
            sys.exit(1)
    # Bin all items for each map by (p,l) pair.
    dmaps = [ {} for m in maps ]
    for mndx,m in enumerate(maps):
        for item in m:
            inputs, outputs, tags = item
            dmaps[mndx][tuple(inputs)] = [outputs,tags]
        print("  Map%d: %d items"%(mndx,len(dmaps[mndx])))
    # Make sure maps have the same items.
    for i,idmap in enumerate(dmaps):
        for ikey in idmap:
            for j,jdmap in enumerate(dmaps):
                if ikey not in jdmap:
                    print("Item from map %d not in map %d!"%(i,j))
                    sys.exit(1)
    # Each item has "[outputs,tags]", augment with ensemble info.
    ensemble = {}
    ensemble_sum = {}
    for key in dmaps[0]:
        actual = 1.0 if float(dmaps[0][key][0][0]) == 1.0 else 0.0
        binds = 0
        ensemble_sum[key] = [0.0, 0.0, actual]
        for dmap in dmaps:
            outputs, tags = dmap[key]
            prediction = 1.0 if "pred:1" in tags else 0.0
            if prediction == 1.0:
                binds += 1
            pred0 = float(tags[2].split(":")[1])
            pred1 = float(tags[3].split(":")[1])
            ensemble_sum[key][0] += pred0
            ensemble_sum[key][1] += pred1
        ensemble[key] = (key, actual, binds)
    ensemble_sum = sorted(ensemble_sum.values(), key=lambda e: e[0])
    # Get some totals / global stats.
    total_binds = 0
    total_nobinds = 0
    for key in ensemble:
        k, actual, binds = ensemble[key]
        if actual == 1.0:
            total_binds += 1
        else:
            total_nobinds += 1
    print("  actual: %.2f%% binds (%d/%d)"%(100.0*(total_binds/(total_binds+total_nobinds)),total_binds,total_nobinds))
    # Float combination
    apct = 0.0
    acnt = 0
    bcnt = 0
    with open("accuracy.txt","w") as accfile:
        for ndx,item in enumerate(ensemble_sum):
            pred0, pred1, actual = item
            acnt += 1
            if actual == 1.0:
                bcnt += 1
            if actual == 0.0 and pred0 > pred1:
                apct += 1
            elif actual == 1.0 and pred0 < pred1:
                apct += 1
            if ndx % int(len(ensemble_sum)/100) == 0:
                accfile.write("%f %f %f\n"%(100.0*float(ndx)/len(ensemble_sum),100.0*apct/acnt,bcnt))
    apct = 100.0 * apct / acnt
    print("Full ensemble float accuracy: %.2f%c"%(apct,'%'))
    # Build a dict based on bind prediction counts.
    counts = { i:[] for i in reversed(range(len(dmaps)+1)) }
    for item in sorted(ensemble.values(), key=lambda e: e[2]):
        key, actual, binds = item
        counts[binds].append(item)
    # Build some accuracy numbers.
    apct = 0.0
    acnt = 0
    for count in counts:
        for item in counts[count]:
            key, actual, binds = item
            acnt += 1
            if count == len(dmaps):
                if actual == 1.0:
                    apct += 1
            elif count != len(dmaps):
                if actual == 0.0:
                    apct += 1
    apct = 100.0 * apct / acnt
    print("Full ensemble accuracy: %.2f%c"%(apct,'%'))
    # Print ensemble info for each prediction count.
    cum_cbinds = 0
    cum_ncount = 0
    for count in counts:
        ncount = len(counts[count])
        ncpct = 100.0 * ncount / len(maps[0])
        cbinds = 0
        for item in counts[count]:
            key, actual, binds = item
            if actual == 1.0:
                cbinds += 1
        bpct = 100.0 * cbinds / ncount
        cum_cbinds += cbinds
        cum_ncount += ncount
        cbpct = 100.0 * cum_cbinds / cum_ncount
        if total_binds == 0:
            IEF = -1.0
            EF = -1.0
        else:
            IEF = 1.0 / (float(total_binds) / (total_binds+total_nobinds))
            EF = IEF * (float(cum_cbinds) / cum_ncount) 
        print("  count%d: %d  (%.1f%%)  actual_bind%%: %.1f (%.1f cum)  EF(cum): %.2f (ideal %.2f)"%(count,ncount,ncpct,bpct,cbpct,EF,IEF))
    print("  count*: %d  (100.0%%)"%(len(maps[0])))
    return


############################################################


def parse_args():
    # Parse command line args.
    parser = argparse.ArgumentParser(prog='ensemble.py: Combine inference outputs into an ensemble.', description='Combine inference outputs into an ensemble.')
    parser.add_argument('--maps', type=str, required=True, help='Colon-seperated list of .map file paths.')
    args = parser.parse_args()
    print(args)
    # Return parsed args.
    return args


def main():
    # Parse command-line args.
    args = parse_args()
    map_fns = args.maps.split(':')
    # Verbose print.
    print("Creating ensemble from %d map files:"%(len(map_fns)))
    for map_fn in map_fns:
        print("  %s"%(map_fn))
    # Read the map files.
    maps = [ read_map(map_fn) for map_fn in map_fns ]
    # Write out a map file representing the ensemble.
    ensemble_maps(maps)
    print("Success!")

    
if __name__== "__main__":
    main()


############################################################
