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

#
# Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or 'snr'
#
if par == 'MBD':
    parname = 'mbdelay'
elif par == 'SBD':
    parname = 'sbdelay'
else:
    parname = 'snr'

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
# To start processing from istart;  exclude bad data before istart.
#
istart = 2

tim = {}      # Time points. The gaps will be replaced with NaNs
tim1 = {}     # Original time points with some of them missing. 
trul = 605.*np.arange(35) # Time ruler
par_l = {}
par_c = {}

for bl in bls:
    tim[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
    tim1[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
    tim[bl][13:] = tim[bl][13:] + 605.
    tim1[bl][13:] = tim[bl][13:] + 605.
    
    #tim[bl] = np.array(idx3819l_1[bl]['I']['time']) # / 60 # Sec -> min
    tim[bl] = tim[bl] - tim[bl][0]  # Set time start at zero
    tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

    snr_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
    snr_a = (abs(snr_l.mean()) + abs(snr_c.mean()))/2 # Avg Lin and Cir means
    
    if par == 'MBD' or par == 'SBD':
        par_l_us = np.array(idx3819l_1[bl]['I'][parname])[istart:] # In useconds
        par_c_us = np.array(idx3819c_1[bl]['I'][parname])[istart:] # In useconds
        par_l[bl] = par_l_us*1e6           # Convert us to ps
        par_c[bl] = par_c_us*1e6           # Convert us to ps
    else: # if par == 'SNR':
        par_l[bl] = np.copy(snr_l[bl])
        par_c[bl] = np.copy(snr_c[bl])
    #
    # Insert NaNs in the time gaps (ie over 605. s away)
    # Accordingly, insert NaNs in the parameter arrays
    #
    irul = 0
    itim = 0
    itim1 = 0
    for irul in range(35): ????????????????????????????????????????????????
        dt = tim1[bl][itim1] - trul[irul]
        if dt > 0:  # The current, itim'th, value jumps over ruler
            ndt = int(dt/605) # Number of time jumps
            for idt in range(ndt):
                tim[bl] = np.insert(tim[bl], itim, np.NaN)
                par_l[bl] = np.insert(par_l[bl], itim, np.NaN)
                par_c[bl] = np.insert(par_c[bl], itim, np.NaN)
                itim = itim + 1
        else:
            itim = itim + 1
            itim1 = itim1 + 1
        if itim >= len(tim1[bl]):
            break

sh = 0
pl.figure()
for bl in bls:
    pl.plot(tim[bl] + sh)
    pl.plot(tim[bl] + sh, '.', markersize=3)
    sh = sh + 1000
    
pl.grid(True)

pl.show()

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

