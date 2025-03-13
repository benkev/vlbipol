help_text = '''
plot_scan_times_2187.py - Visualize time gaps for each baseline,
                          session 2187 (VO2187).
Usage:
%run plot_scan_times_2187.py [save]

    save is optional: save figures in pdf format.
'''

plotColorLegend =   False
plotAvailableTime = True

import sys

sf = False  # Save figure request

if len(sys.argv) > 1:
    if sys.argv[1] == '--help':
        print(help_text)
        sys.exit(0)

    if sys.argv[1] == 'save':
        sf = True
    else:
        print("Argument can only be '--help' or 'save'. Entered '%s'." %
              sys.argv[1])
        sys.exit(0)

import pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm

#
# Unpickle the index dictionary for linear pol:
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

    

tim1 = {}     # Original time points with some of them missing. 
tim = {}      # Time points. The gaps will be replaced with NaNs

for bl in bls:
    tim1[bl] = np.array(idxl[bl]['I']['time'])  #/ 3600 # Sec -> hours
    tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

sys.exit(0)

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
# Search for the maximum number of time counts, ntim,
# among all the baselines. At the end of the loop below, ntim will hold
# the length for the parameter storage in arrays.
# The parameter sets shorter than ntim contain time gaps to be filled with NaNs.
#
ntim = -1       # Will contain the maximum time count
for bl in bls:
    bl_ntim = np.max(tim1[bl]/min_t_scan)    # Max t counts for the baseline
    #print("bl_ntim = %f" % bl_ntim)
    if bl_ntim > ntim:
        ntim = np.int64(np.ceil(bl_ntim))

    print("len(tim1['%s']) = %d; Max t counts = %f" %
                                   (bl, len(tim1[bl]), bl_ntim))

print("Max time counts: %d;  min scan time: %d s." % (ntim, min_t_scan))
    
#
# Insert NaNs in the time gaps
#
for bl in bls:
    itim = np.int64(tim1[bl]/min_t_scan) # Indices of non-NaN elements into tim
    tim[bl] = np.zeros(ntim)    
    tim[bl][:] = np.nan
    tim[bl][itim] = tim1[bl]



#
# Plot available times for each baseline
#
if plotAvailableTime:
    #cols_bl = cm.rainbow(np.linspace(0, 1, nbls))
    #cols_bl = cm.nipy_spectral(np.linspace(0, 1, nbls))
    #cols_bl = cm.gist_rainbow(np.linspace(0, 1, nbls))
    cols_bl = cm.jet(np.linspace(0, 1, nbls))

    fig4, ax41 =  pl.subplots()

    fig4.text(0.22, 0.95, "Baseline Times with Missed Scans", fontsize=14)

    #sh = 1 # Just arbitrary shift to plot the lines 
    sh = np.ones(ntim) # Horizontal line with gaps
    for ib in range(nbls):
        bl = bls[ib]
        t = tim[bl]/3600
        y = nbls - ib
        yy = y*sh                    # Array of heights
        ax41.plot(t, yy, color=cols_bl[ib,:], lw=3)
        ax41.plot(t, yy, 'k.', markersize=5)
        ax41.text(-0.55, y-0.35, bl, fontsize=14)
        print("ib = %d, y = %d" % (ib, y))

    ax41.grid(True)
    ax41.set_xlabel("hours", fontsize=14)
    ax41.set_yticks([])
    ax41.set_ylim(0, nbls+1)
    ax41.set_xlim(-0.8, 24)
    fig4.tight_layout(rect=(0.00, 0.00, 0.98, 0.95))

    pl.savefig("Gaps_in_Time.pdf", format='pdf')


pl.show()

#
# Save figures on request
#
if sf:
    pl.figure(fig4)
    pl.savefig("VO2187_Baseline_times_with_missed_scans.pdf", format='pdf')





