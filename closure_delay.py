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


ntri = nbls*(nbls-1)*(nbls-2)//(1*2*3)   # Number of baseline triangles
tau = np.zeros((), dtype=int)

#
# Loop over the baseline triangles (bla, blb, blc) to find closure delays
#
for bla, blb, blc in combinations(bls, 3): 
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]

    para_l = np.array(idx3819l_1[bla]['I'][parname])[istart:]*1e6 # In pseconds
    parb_l = np.array(idx3819l_1[blb]['I'][parname])[istart:]*1e6 # In pseconds
    parc_l = np.array(idx3819l_1[blc]['I'][parname])[istart:]*1e6 # In pseconds

    print(para_l + parb_l + parc_l)
    
    #par_c = np.array(idx3819c_1[bl]['I'][parname])[istart:]*1e6 # In pseconds

