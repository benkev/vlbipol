import os, sys, copy, pickle
import pickle, copy
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
# from libvp import find_tail_bounds, group_tails
import libvp  # Needs to reset backend for vpal sets it to non-interactive Agg!
from libplt import plot_closure_legend
from libvp import find_tail_bounds



parname = 'mbdelay'


with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
with open('idx2187cI.pkl', 'rb') as finp: idxc = pickle.load(finp)

with open('idxs2187lI.pkl', 'rb') as finp: idxsl = pickle.load(finp)
with open('idxs2187cI.pkl', 'rb') as finp: idxsc = pickle.load(finp)

with open('idxf2187lI.pkl', 'rb') as finp: idxfl = pickle.load(finp)
with open('idxf2187cI.pkl', 'rb') as finp: idxfc = pickle.load(finp)

with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)
with open('clos2187cI.pkl', 'rb') as finp: closc = pickle.load(finp)

with open('clot2187lI.pkl', 'rb') as finp: clotl = pickle.load(finp)
with open('clot2187cI.pkl', 'rb') as finp: clotc = pickle.load(finp)

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

#
# Make dict of colors cols[triangle] --> 4-array of a color
#
acols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri)) # Array of colors
cols = {}
for ic in range(ntri):
    tr = trians[ic]
    cols[tr] = acols[ic]

upar = parname.upper()

#
# Plotting for the 1803+784 source
#
sr = '1803+784'


gs_kw1 = dict(width_ratios=[0.75, 0.25], height_ratios=[0.15, 0.425, 0.425])
fig1, ax1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                               ['dis_lin', 'his_lin'],
                               ['dis_cir', 'his_cir']],
                               gridspec_kw=gs_kw1, figsize=(8.4, 10))
fig1.tight_layout(rect=(0.00, 0.00, 0.99, 0.99))

#
# Plot color legend on top
#
ax_col = ax1['col_legend']    # Get the axis for color legend

plot_closure_legend(ax_col, cols, upar, fs=9)  
ax_col.set_title("VO2187 Source 1803+784 %s Closure" % upar, fontsize=12)
hist_colr = 'red'

axdl = ax1['dis_lin']

tr = 'EGS'
cl = closl[sr][tr]
tim = closl[sr][tr]['thour']
tau = closl[sr][tr]['tau_mbd']

axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols['EGS'])
axdl.set_title("%s, Linear" % tr, fontsize=10)
axdl.grid(1)

axdc = ax1['dis_cir']
tr = 'EGT'
tim = closl[sr][tr]['thour']
tau = closl[sr][tr]['tau_mbd']

axdc.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols['EGS'])
axdc.set_title("%s, Linear" % tr, fontsize=10)
axdc.grid(1)

# axhl = ax1['his_lin']
# axhl.hist(tcl, 51, orientation='horizontal', color='r')
# axhl.grid(1)
# axhl.set_title("%d pts" % npt, fontsize=10)
# if clonm == 'cloph': axhl.set_xlim(-1,50) 




sys.exit(0)


#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
    if tr in ['EGS', 'EGT']:
        nt = len(closc[sr][tr]['thour'])
        print("nt = %d, closc[%8s][%3s]['thour']" % (nt, sr, tr)) 
        print("N triangle closc[sr][tr]['thour']  cl['thour']  " \
              "closc[sr][tr]['tau_mbd']  cl[clonm]")
        for i in range(nt):
              print("%3d %3s  %10.4f  %10.4f  %10.4f  %10.4f" % \
                    (i, tr, closc[sr][tr]['thour'][i], cl['thour'][i],
                     closc[sr][tr]['tau_mbd'][i], cl[clonm][i]))
    axdc.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcc = np.append(tcc, cl[clonm])
    thc = np.append(thc, cl['thour'])
    npt = npt + len(cl['thour'])
axdc.grid(1)

axdc.set_title("Circular", fontsize=10)


axhc = ax1['his_cir']
axhc.hist(tcc, 51, orientation='horizontal', color='r')
axhc.grid(1)
axhc.set_title("%d pts" % npt, fontsize=10)
axhc.set_xlim(-1,50)

pl.show()

if sf: pl.savefig("VO2187_Source_1803+784_%s_Closure.pdf" % upar, format='pdf')







