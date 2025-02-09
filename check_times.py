import os, sys
import re, glob
import pickle
import numpy as np
import matplotlib.pyplot as pl

if len(sys.argv) < 3  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python check_times <expm>  <parname> [save], ")
    print("       where <expm> is the 4-digit experiment number (like 3819),")
    print("       and <parname> is either MBD or SBD or SNR.")
    print("       save (optional): save  figures in pdf format.")
    sys.exit(0)


expm = sys.argv[1]
if not re.match("[0-9]{4}" , expm):
    print("Only 4-digit experiment numbers allowed. Entered " + expm +
          ". Exiting.")
    sys.exit(1)
    
#
# Find the indices with pseudo-Stokes pol prods (file names with 'lI' and 'cI')
#
fidxl = glob.glob("idx" + expm + "lI.pkl")  # Linear pol with pseudo-Stokes I
fidxc = glob.glob("idx" + expm + "cI.pkl")  # Circular pol with pseudo-Stokes I

no_idxl = False
if fidxl is []:
    no_idxl = True
elif not re.match(r"idx[0-9]{4}lI\.pkl", fidxl[0]):
    no_idxl = True
if no_idxl:
    print("No linear polprod data file idx" + expm + "lI.pkl.  Exiting.")
    sys.exit(0)

no_idxc = False
if fidxc is []:
    no_idxc = True
elif not re.match(r"idx[0-9]{4}cI\.pkl", fidxc[0]):
    no_idxc = True
if no_idxc:
    print("No circular polprod data file idx" + expm + "cI.pkl.  Exiting.")
    sys.exit(0)
    
fidxl = fidxl[0] 
fidxc = fidxc[0] 

    
# sys.exit(0)

parname = sys.argv[2]
parname = parname.upper()
if parname != 'MBD' and parname != 'SBD' and parname != 'SNR':
    print("Argument can be MBD or SBD or SNR. Entered '%s'. Exiting." %
          sys.argv[2])
    sys.exit(1)

sf = False  # Request to save figures
if len(sys.argv) == 4:
    if sys.argv[3] == 'save':
        sf = True
    else:
        print("Argument can only be 'save'. Entered '%s'. Exiting." %
              sys.argv[3])
        sys.exit(1)
    
#
# Determine the parameter 'par': 'mbdelay', 'sbdelay', or 'snr'
#
if parname == 'MBD':
    par = 'mbdelay'
elif parname == 'SBD':
    par = 'sbdelay'
else:
    par = 'snr'

if parname == 'MBD' or parname == 'SBD':
    ps = "(ps)"
else: # if parname == 'SNR':
    ps = ""
    
    
#
# Unpickle the indices:
#
with open(fidxl, 'rb') as finp:  # Linear pols index file
    idxl = pickle.load(finp)

with open(fidxc, 'rb') as finp:  # Circular pols index file
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
nbls = len(bls)

print("\nbl, tl, tc:")
for bl in bls:   # Loop over the baselines
    tl = np.array(idxl[bl]['I']['time']) / 60; tl = tl - tl[0]
    tc = np.array(idxl[bl]['I']['time']) / 60; tc = tc - tc[0]
    print("%2s %4d %4d %r" % (bl, len(tl), len(tc), len(tl) == len(tc)))

print("\nbl, tl, snr_l, snr_c, eq?, par_l, par_c, eq?:")
for bl in bls:   # Loop over the baselines
    tl = np.array(idxl[bl]['I']['time']) / 60; tl = tl - tl[0]
    tc = np.array(idxl[bl]['I']['time']) / 60; tc = tc - tc[0]
    snr_l = np.array(idxl[bl]['I']['snr'])
    snr_c = np.array(idxc[bl]['I']['snr'])
    par_l = np.array(idxl[bl]['I'][par])
    par_c = np.array(idxc[bl]['I'][par])
    print("%2s %4d %4d %4d %5r %4d %4d %5r" %
          (bl, len(tl), len(snr_l), len(snr_c), len(snr_l) == len(snr_c),
           len(par_l), len(par_c), len(par_l) == len(par_c)))



