'''
save_closure_dics.py -- Create closure dictionaries for both polarizations,
    linear and circular, from the idxsl and idxsc dicts:
        closl[sr][tr][di],
        closc[sr][tr][di],
        clotl[tr][sr][di]
        clotc[tr][sr][di].
    Also, from the dics idxl and idxs, loaded from disk, find the baselines,
    common for both linear and circular polproducts in the bls list and
    the tribl dict of baseline triangles, based on bls.

    Pickle them and save on disk.

    Note: this script is actually the tail part of the make_idx_2187.py script.

'''

import os, sys, pickle, copy
from libvp import find_baseline_triangles, make_closure_dic, clos_to_clot

with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
with open('idx2187cI.pkl', 'rb') as finp: idxc = pickle.load(finp)

with open('idxs2187lI.pkl', 'rb') as finp: idxsl = pickle.load(finp)
with open('idxs2187cI.pkl', 'rb') as finp: idxsc = pickle.load(finp)


# =======================================================================
# ======= Creation of dictionaries with all the closures possible =======
# ======= as clos[src][tri][data_item] (closl and closc),         =======
# ======= where data_items are 'cloph', 'tau_mbd', 'tau_sbd' etc. =======
# =======================================================================
#
#
# Circular polarization data use FEWER baselines than linear pol.
# 
# This is because PolConvert, by some reason, omitted the 'Y' (i.e. 'Yj')
# station, so it is not present in the baseline list.
#
# Below we read both linear and circular baselines and select only those
# baselines that are present in both cases.
#
bls_l = set(idxl.keys())    # Linear baselines (set)
nbls_l = len(bls_l)
bls_c = set(idxc.keys())    # Circular baselines (set)
nbls_c = len(bls_c)

bls = list(bls_l & bls_c)   # Find the lin and cir sets intersection 
bls.sort()                  # Sort baselines lexicographically

#
# Exclude the 'ST' baseline: the S and T stations are too close
# to each other
#
if 'ST' in bls:
    iST = bls.index('ST')
    bls.pop(iST)

nbls = len(bls)

print("ST baseline excluded\n")
print("Found %d baselines common for linear and circular polarization:" %
      nbls)
print(bls, '\n')
print()

tribl = find_baseline_triangles(bls)

# with open('bls_2187.pkl', 'wb') as fout: pickle.dump(bls, fout)
# print("Baseline list pickled and saved on disk in bls_2187.pkl")
# print()

# with open('tribl_2187.pkl', 'wb') as fout: pickle.dump(tribl, fout)
# print("Dict of baseline triangles pickled and saved on disk in " \
#       "tribl_2187.pkl")
# print()


#
# Create dictionaries of closures, clos[sr][tr][di] and clot[tr][sr][di]
# for both polarizations, linear and circular
#

closl = make_closure_dic(idxsl, bls)
closc = make_closure_dic(idxsc, bls)

print("Created dictionary of data closures, closl, linear polarization")
print("Created dictionary of data closures, closc, circular polarization\n")

clotl = make_closure_dic(closl, tribl)
clotc = make_closure_dic(closc, tribl)

print("Created dictionary of data closures, clotl, linear polarization")
print("Created dictionary of data closures, clotc, circular polarization\n")

# with open('clos2187lI.pkl', 'wb') as fout: pickle.dump(closl, fout)
# with open('clos2187cI.pkl', 'wb') as fout: pickle.dump(closc, fout)

# with open('clot2187lI.pkl', 'wb') as fout: pickle.dump(clotl, fout)
# with open('clot2187cI.pkl', 'wb') as fout: pickle.dump(clotc, fout)

# print("Linear pol. data closure dictionary pickled and saved on disk")
# print("Circular pol. data closure dictionary pickled and saved on disk\n")

print("To load from disk:")
print('''
import pickle

with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)
with open('clos2187cI.pkl', 'rb') as finp: closc = pickle.load(finp)

with open('clot2187lI.pkl', 'rb') as finp: clotl = pickle.load(finp)
with open('clot2187cI.pkl', 'rb') as finp: clotc = pickle.load(finp)

Load bls, list of baselines common for linear and circular polarization:

with open('bls_2187.pkl', 'rb') as finp: bls = pickle.load(finp)

Load tribl, dict of baseline triangles based on bls,
     tribl[tr] --> (bl1, bl2, bl3): 

with open('tribl_2187.pkl', 'rb') as finp: tribl = pickle.load(finp)

''')
print()




