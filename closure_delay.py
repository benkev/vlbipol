help_text = '''
closure_delay.py - for  MBD or SBD.
'''
import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_time_var.py <par> [save], ")
    print("       where <par> is either MBD or SBD or SNR.")
    print("       save (optional): save  figures in pdf format.")
    sys.exit(0)
    
par = sys.argv[1]
par = par.upper()
if par != 'MBD' and par != 'SBD' and par != 'SNR':
    print("Argument can be MBD or SBD or SNR. Entered '%s'. Exiting." %
          sys.argv[1])
    sys.exit(1)

sf = False  # Save figure request
if len(sys.argv) == 3:
    if sys.argv[2] == 'save':
        sf = True
    else:
        print("Argument can only be 'save'. Entered '%s'. Exiting." %
              sys.argv[2])
        sys.exit(1)
    
    
import pickle
import numpy as np
import matplotlib.pyplot as pl
from itertools import combinations

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

#
# Unpickle it:
#
with open('idx3819l.pkl', 'rb') as finp:
    idx3819l_1 = pickle.load(finp)

# with open('idx3819c.pkl', 'rb') as finp:
#     idx3819c_1 = pickle.load(finp)

with open('idx3819cI.pkl', 'rb') as finp:
    idx3819c_1 = pickle.load(finp)

if par == 'MBD' or par == 'SBD':
    ps = "(ps)"
else: # if par == 'SNR':
    ps = ""
    
bls = list(idx3819l_1.keys())   # Baselines
bls.sort()                      # Lexigraphically sorted baselines
nbls = len(bls)

# Set of station letters stset
ststr = ''
for bl in bls: ststr = ststr + bl  # Concatenate baseline strings in ststr
stset = set(ststr)  # Leave only unique station letters in the sts set

# String of station letters ststr
nsts = len(stset)
# ststr = ''
# for st in stset: ststr = ststr + st
ststr = ''.join(sorted(stset))

#
# Find all the baseline triplets with 3 stations
#
trian = {}
ntri = 0   # Number of baseline triangles
for ab, bc, ca in combinations(bls, 3):
    stset = set(''.join(ab+bc+ca))
    trist = ''.join(sorted(stset))
    #if len(stset) == 3:
    if len(trist) == 3:
        #print(stset)
        print(trist)
        print(ab, bc, ca)
        trian[trist] = trist 
        ntri = ntri + 1   # Number of baseline triangles



tau = np.zeros((ntri,13), dtype=int) # Only first 13 times are the same

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

tim = {}
tim1 = {}
trul = 605.*np.arange(35) # Time ruler
for bl in idx3819l_1.keys():
    tim[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
    tim1[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
    #tim[bl] = np.array(idx3819l_1[bl]['I']['time']) # / 60 # Sec -> min
    tim[bl] = tim[bl] - tim[bl][0]  # Set time start at zero
    tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero
    #
    # Insert NaNs in the time gaps (ie over 605. s away)
    #
    itim = 0
    while itim < 35:
        if tim[bl][itim] != trul[itim]:
            tim[bl] = np.insert(tim[bl], itim, np.NaN)
            itim = itim + 1
        itim = itim + 1
        if itim == len(tim[bl]):
            break

# pl.figure()
# for bl in idx3819l_1.keys():
#     pl.plot(tim[bl])
#     pl.plot(tim[bl], '.', markersize=3)
    
# pl.grid(True)

# pl.show()

sys.exit()


#
# Loop over the baseline triangles (bla, blb, blc) to find closure delays
#
for bla, blb, blc in combinations(bls, 3): 
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]

    para_l = np.array(idx3819l_1[bla]['I'][parname])[istart:]*1e6 # In pseconds
    parb_l = np.array(idx3819l_1[blb]['I'][parname])[istart:]*1e6 # In pseconds
    parc_l = np.array(idx3819l_1[blc]['I'][parname])[istart:]*1e6 # In pseconds

    # tau[para_l + parb_l + parc_l)
    
    #par_c = np.array(idx3819c_1[bl]['I'][parname])[istart:]*1e6 # In pseconds

