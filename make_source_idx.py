#
# make_source_idx.py
#
# Create hash-table (dictionary) idxs:
#    idxs[source][time][baseline]
# For each source and each time it shell have a list of the baselines pointing
# at the source at the time.
# Also, it shall have the index of source into asrc source array:
#    idxs[source] : <its index into asrc>
# 
#

import pickle
import numpy as np
import matplotlib.pyplot as pl

pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

#
# Unpickle it:
#
with open('idx2187lI.pkl', 'rb') as finp:
    idxl = pickle.load(finp)

with open('idx2187cI.pkl', 'rb') as finp:
    idxc = pickle.load(finp)

#
# Circular pol uses FEWER baselines than linear pol.
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
# Exclude the 'ST' baseline: the S and T stations are too close to each other
#
if 'ST' in bls:
    iST = bls.index('ST')
    bls.pop(iST)

nbls = len(bls)

#
# Find the time counts
#
tim1 = {}     # Original time points with some of them missing. 
tim = {}      # Time points. The gaps will be replaced with NaNs

tim_total = []     # List of all the time counts for all the baselines

for bl in bls:
    timbl = idxl[bl]['I']['time']
    tim_total.extend(timbl)
    tim1[bl] = np.array(timbl)  #/ 3600 # Sec -> hours

ttim = np.unique(tim_total)   # Unite all the time counts in one array
ttim0 = ttim[0]    # 'global' time start

# Set time for each baseline start at 'global' zero
for bl in bls:
    tim1[bl] = tim1[bl] - ttim0  # Set time start at 'global' zero
    print("tim1[%s] = %.2f" % (bl, tim1[bl][0]))

    
#
# Find the minimum time between the scans over all the baselines
#
min_t_scan = 1000000000

for bl in bls:
    t_scans = np.diff(tim1[bl])
    bl_min_t_scan = np.min(t_scans)
    if bl_min_t_scan < min_t_scan:
        min_t_scan = bl_min_t_scan






#
# Create array of all the sources
#
srcl = []
srcc = []
for bl in bls:
    srcl.extend(idxl[bl]['I']['source'])
    srcc.extend(idxl[bl]['I']['source'])
rcls = srcl + srcc
asrc = np.unique(rcls)     # Turn into np.array leaving only unique source names
asrc.sort()                # Sort the source names lexicographically

nsrc = len(asrc)    # Number of sources

#
# Create hash-table (dictionary) idxs:
#    idxs[source][time][baseline]
# For each source and each time it shell have a list of the baselines pointing
# at the source at the time.
# Also, it shall have the index of source into asrc source array:
#    idxs[source] : <its index into asrc>
#
idxs = {}

for bl in bls:
    srcs = idxl[bl]['I']['source']
    atms =  np.array(idxl[bl]['I']['time']) - ttim0
    nt = len(atms)

    for i in range(nt):
        sr = srcs[i]
        tm = atms[i]
        
        print(sr, tm)

        if sr in idxs.keys():
            if tm in idxs[sr].keys():
                idxs[sr][tm].append(bl)
                #print(tm)
            else:
                idxs[sr] = {tm : [bl]}
        else:
            # idxs[sr] = {}
            # idxs[sr][tm] = [bl]
            idxs[sr] = {tm : [bl]}














