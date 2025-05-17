'''
plot_ST_closures.py - plot closures for triangles including stations S and T,
                      session VO2187.
'''
import matplotlib
print("matplotlib.get_backend() --> ", matplotlib.get_backend())

import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print('''
    Plot closures for triangles including stations S and T,
                      session VO2187.
    ''')
    print("Usage:")
    print("python plot_ST_closures.py <par> [save], ")
    print("       <par> is one of phase, mbd, sbd, tmbd, or tsbd")
    print("       save (optional): save  figures in pdf format.")
    sys.exit(0)

arg_to_par = {'phase':'phase', 'mbd':'mbdelay', 'sbd':'sbdelay',
              'tmbd':'tot_mbd', 'tsbd':'tot_sbd'}

arg_to_clo = {'phase':'cloph', 'mbd':'tau_mbd', 'sbd':'tau_sbd',
              'tmbd':'tau_tmbd', 'tsbd':'tau_tsbd'}

parlist = list(arg_to_par.keys())
pararg = (sys.argv[1]).lower()
if pararg not in arg_to_par.keys():
    print("Argument can be one of %r. Entered '%s'. Exiting." %
          (parlist, sys.argv[1]))
    sys.exit(0)

sf = False  # Request to save figure 
if len(sys.argv) == 3:
    if sys.argv[2] == 'save':
        sf = True
    else:
        print("Argument can only be 'save'. Entered '%s'. Exiting." %
              sys.argv[2])
        sys.exit(0)
    
parname = arg_to_par[pararg]           # Like 'mbdelay' for 'mbd' etc.
upar = parname.upper()
clonm = arg_to_clo[pararg]

import os, sys, copy, pickle
import numpy as np
import matplotlib.pyplot as pl
# from matplotlib.pyplot import cm
# from libvp import find_tail_bounds, group_tails
import libvp  # Needs to reset backend for vpal sets it to non-interactive Agg!
# from libplt import plot_closure_legend
# from libvp import find_tail_bounds



parname = 'mbdelay'


# with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
# with open('idx2187cI.pkl', 'rb') as finp: idxc = pickle.load(finp)

# with open('idxs2187lI.pkl', 'rb') as finp: idxsl = pickle.load(finp)
# with open('idxs2187cI.pkl', 'rb') as finp: idxsc = pickle.load(finp)

# with open('idxf2187lI.pkl', 'rb') as finp: idxfl = pickle.load(finp)
# with open('idxf2187cI.pkl', 'rb') as finp: idxfc = pickle.load(finp)

with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)
with open('clos2187cI.pkl', 'rb') as finp: closc = pickle.load(finp)

# with open('clot2187lI.pkl', 'rb') as finp: clotl = pickle.load(finp)
# with open('clot2187cI.pkl', 'rb') as finp: clotc = pickle.load(finp)

with open('bls_2187.pkl', 'rb') as finp: bls = pickle.load(finp)
with open('tribl_2187.pkl', 'rb') as finp: tribl = pickle.load(finp)

trians = list(tribl.keys())
ntri = len(trians)

trist = []
for tr in trians:
    if 'S' in tr or 'T' in tr:
        trist.append(tr)

print("Triangles with S and T stations:")
print(trist)



sr = '1803+784'

trS = ['EGS', 'GHS', 'GIS', 'GMS', 'EHS', 'HMS', 'HIS', 'EIS', 'IMS', 'EMS']
trT = ['EGT', 'GHT', 'GIT', 'GMT', 'EHT', 'HMT', 'HIT', 'EIT', 'IMT', 'EMT']

ntr = len(trS)

for itr in range(ntr):

    trs = trS[itr]
    th = closl[sr][trs]['thour']
    cl = closl[sr][trs][clonm]

    pl.figure()

    pl.plot(th, cl, 'r.', ms=8, label=trs)

    trt = trT[itr]
    th = closl[sr][trt]['thour']
    cl = closl[sr][trt][clonm]

    pl.plot(th, cl, 'g.', ms=8, label=trt)
    pl.title("VO2187 Source %s: %s Closure for %s and %s" %
             (sr, upar, trs, trt))
    pl.grid(1)

    pl.xlabel("t (hours)")
    pl.ylabel(clonm)
    pl.legend()
    
    pl.show()
   
    if sf: pl.savefig("VO2187_Source_%s_%s_Closure_for_%s_and_%s.pdf" %
               (sr, upar, trs, trt), format='pdf')








