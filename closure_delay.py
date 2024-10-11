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

arg_to_par = {'mbd':'mbdelay', 'sbd':'sbdelay', 'snr':'snr', 'tmbd':'tot_mbd',
              'tsbd':'tot_sbd', 'rmbd':'resid_mbd', 'rsbd':'resid_sbd'}
    
arg1 = sys.argv[1]

# arg1 = arg1.upper()
# if arg1 != 'MBD' and arg1 != 'SBD' and arg1 != 'SNR':

arg1 = arg1.lower()
if arg1 not in arg_to_par.keys():
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
import copy

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
# Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or 'snr' or
#     'tmbd' or 'tsbd' or 'rmbd' or 'rsbd'
#
parname = arg_to_par[arg1]

ps = "(ps)"
if arg1 == 'snr': ps = "(ps)"
    
# if par == 'MBD':
#     parname = 'mbdelay'
# elif par == 'SBD':
#     parname = 'sbdelay'
# else:
#     parname = 'snr'

# if par == 'MBD' or par == 'SBD':
#     ps = "(ps)"
# else: # if par == 'SNR':
#     ps = ""
    

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
# Find all the baseline triplets with 3 stations (ie station triangles)
#
# trist: a string of three station letters, like 'EMS', 'MSY', 'TVY' etc.
#
#trian = {}
trians = [] # List of station riangles: 'EVY', 'EMV', 'ESV', 'ETV' etc
tribl = {}  # Dict to translate triangle to baselines, like MSV -> MS, SV, MV
ntri = 0   # Number of baseline triangles
for ab, bc, ca in combinations(bls, 3):
    stset = set(''.join(ab+bc+ca))
    trist = ''.join(sorted(stset))
    #print("trist = ", trist)
    #if len(stset) == 3:
    if len(trist) == 3:
        #print(stset)
        #print(trist)
        #print(ab, bc, ca)
        #trian[trist] = trist
        trians.append(trist)
        tribl[trist] = (ab, bc, ca)
        ntri = ntri + 1   # Number of baseline triangles

#
# Reorder tuples of baselines in tribl into pattern 'AB', 'BC', 'AC' 
#

tribl1 = {}
for trist in tribl.keys():
    l3bl = tribl[trist]  # List of the 3 baselines
    trian = l3bl
    
    st_end = l3bl[0][1]
    if l3bl[2][0] == st_end:
        trian = (l3bl[0], l3bl[2], l3bl[1])
        tribl1[trist] = trian
        print(tribl[trist], '->', trian)

    st_end = l3bl[1][1]
    if l3bl[0][0] == st_end:
        trian = (l3bl[1], l3bl[0], l3bl[2])
        tribl1[trist] = trian
        print(tribl[trist], '->', trian)
    elif l3bl[2][0] == st_end:
        trian = (l3bl[1], l3bl[2], l3bl[0])
        tribl1[trist] = trian
        print(tribl[trist], '->', trian)

tribl = tribl1
        
#sys.exit()


#
# To start processing from istart;  exclude bad data before istart.
#
istart = 2

tim = {}      # Time points. The gaps will be replaced with NaNs
tim1 = {}     # Original time points with some of them missing. 
# trul = 605.*np.arange(35) # Time ruler
par_l = {}
par_c = {}
snr_l = {}
snr_c = {}
snr_a = {}

#
# The arrays to contain the same data as the dictionaries par_l ans par_c etc.
#
atim   = np.zeros((ntri,35), dtype=float)
apar_l = np.zeros((ntri,35), dtype=float)
apar_c = np.zeros((ntri,35), dtype=float)
asnr_l = np.zeros((ntri,35), dtype=float)
asnr_c = np.zeros((ntri,35), dtype=float)

itri = 0
for bl in bls:
    tim1[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
    tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

    snr1_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr1_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
    snr1_a = (abs(snr1_l.mean()) + abs(snr1_c.mean()))/2 # Avg Lin and Cir
    
    if parname == 'snr':
        par1_l = np.copy(snr1_l)
        par1_c = np.copy(snr1_c)
    else:
        par_l_us = np.array(idx3819l_1[bl]['I'][parname])[istart:] # In useconds
        par_c_us = np.array(idx3819c_1[bl]['I'][parname])[istart:] # In useconds
        par1_l = par_l_us*1e6           # Convert us to ps
        par1_c = par_c_us*1e6           # Convert us to ps
        
    #
    # Insert NaNs in the time gaps (ie over 605. seconds away)
    # Accordingly, insert NaNs in the parameter arrays
    #
    itim = np.int64(tim1[bl]/605) # Indices of non-NaN elements into tim and par
    tim[bl] = np.zeros(35)    
    tim[bl][:] = np.NaN
    tim[bl][itim] = tim1[bl]
    atim[itri,:] = tim[bl]
    
    snr_l[bl] = np.zeros(35)    
    snr_l[bl][:] = np.NaN
    snr_l[bl][itim] = snr1_l
    asnr_l[itri,:] = snr_l[bl]

    snr_c[bl] = np.zeros(35)    
    snr_c[bl][:] = np.NaN
    snr_c[bl][itim] = snr1_c
    asnr_c[itri,:] = snr_c[bl]

    par_l[bl] = np.zeros(35)    
    par_l[bl][:] = np.NaN
    par_l[bl][itim] = par1_l
    apar_l[itri,:] = par_l[bl]

    par_c[bl] = np.zeros(35)    
    par_c[bl][:] = np.NaN
    par_c[bl][itim] = par1_c
    apar_c[itri,:] = par_c[bl]

    itri = itri + 1

sh = 0 # Just arbitrary shift to splot the lines 
pl.figure()
for bl in bls:
    pl.plot(tim[bl] + sh)
    pl.plot(tim[bl] + sh, '.', markersize=3)
    sh = sh + 1000
    
pl.grid(True)

# sys.exit()


#
# Loop over the baseline triangles (bla, blb, blc) to find closure delays
#
tau_l = {} # Dict trist : array of closure delays (ab+bc+ca) for 35 times
tau_c = {} # Dict trist : array of closure delays (ab+bc+ca) for 35 times

#
# Arrays for tau
#
atau_l = np.zeros((ntri,35), dtype=float)
atau_c = np.zeros((ntri,35), dtype=float)

itri = 0
for trist in trians:
    print(trist, ': ', tribl[trist])
    ab, bc, ac = tribl[trist]
    tau_l[trist] = par_l[ab] + par_l[bc] - par_l[ac]
    tau_c[trist] = par_c[ab] + par_c[bc] - par_c[ac]

    atau_l[itri,:] = np.copy(tau_l[trist])
    atau_c[itri,:] = np.copy(tau_c[trist])
    
    itri = itri + 1

upar = parname.upper()

pl.figure()
for trist in tribl.keys():
    pl.plot(tim['VY']/3600, tau_l[trist], '.')
pl.grid(1)
pl.title("%s Closure Delay (Linear) vs Time" % upar)
pl.xlabel("hr")
pl.ylabel("ps")
pl.savefig("%s_Closure_Delay_(Linear)_vs_Time.pdf" % upar, format='pdf')

pl.figure()
for trist in tribl.keys():
    pl.plot(tim['VY']/3600, tau_c[trist], '.')
pl.grid(1)
pl.title("%s Closure Delay (Circular) vs Time" % upar)
pl.xlabel("hr")
pl.ylabel("ps")
pl.savefig("%s_Closure_Delay_(Circular)_vs_Time.pdf" % upar, format='pdf')


pl.figure()
pl.hist(abs(atau_l.flatten()), 100)
pl.xlabel("ps")
pl.title("%s Abs Magnitude Closure Delay Distribution (Linear)" % upar)
pl.savefig("%s_Abs_Magnitude_Closure_Delay_Distribution_(Linear).pdf" % upar,
           format='pdf')


pl.figure()
pl.hist(abs(atau_c.flatten()), 100)
pl.xlabel("ps")
pl.title("%s Abs Magnitude Closure Delay Distribution (Circular)" % upar)
pl.savefig("%s_Abs_Magnitude_Closure_Delay_Distribution_(Circular).pdf" % upar,
           format='pdf')








pl.show()








