'''
plot_closures_2187.py - plot closure phase and closure delay for residual and
                             total MBD or SBD, session 2187 (VO2187).
'''
import matplotlib
print("1. matplotlib.get_backend() ", matplotlib.get_backend())

import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print('''
    Plot closure phase and closure delay for residual and
                             total MBD or SBD, session 2187 (VO2187).
    ''')
    print("Usage:")
    print("python plot_closures.py <par> [save], ")
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
clonm = arg_to_clo[pararg]

# matplotlib.use('qtagg')  # Force the interactive backend

import pickle, copy
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
# from libvp import find_tail_bounds, group_tails
import libvp  # Needs to reset backend for vpal sets it to non-interactive Agg!
from libplt import plot_closure_legend
from libvp import find_tail_bounds

# pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.

# print("2. matplotlib.get_backend() ", matplotlib.get_backend())

# matplotlib.use('qtagg', force=True)  # Force reset the backend due to vpal

def denoise(tm, var):
    '''
    Remove outliers in var using its histogram
    '''
    nv = len(var)
    ni, be = np.histogram(tcl, nv)
    



#
# Unpickle it:
#
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

# bls = list(idxc.keys())    # List of baselines 
# bls.sort()                 # Sort baselines lexicographically

# #
# # Exclude the 'ST' baseline: the S and T stations are too close to each other
# #
# if 'ST' in bls:
#     iST = bls.index('ST')
#     bls.pop(iST)

nbls = len(bls)

# tribl = libvp.find_baseline_triangles(bls)

trians = list(tribl.keys())
ntri = len(trians)

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
#ttl_lin = "Fourfit Pseudo-I, %s vs Time (%d triangles)" # % (upar, ntri)

# cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    if tr in ['EGS', 'EGT']:
        nt = len(closl[sr][tr]['thour'])
        print("nt = %d, closl[%8s][%3s]['thour']" % (nt, sr, tr)) 
        print("N triangle closl[sr][tr]['thour']  cl['thour']  " \
              "closl[sr][tr]['tau_mbd']  cl[clonm]")
        for i in range(nt):
              print("%3d %3s  %10.4f  %10.4f  %10.4f  %10.4f" % \
                    (i, tr, closl[sr][tr]['thour'][i], cl['thour'][i],
                     closl[sr][tr]['tau_mbd'][i], cl[clonm][i]))
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-1,50) 

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

sys.exit(0)

#
# Plotting for the 2229+695 source
#

sr = '2229+695'

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

ax_col.set_title("VO2187 Source %8s %s Closure" % (sr, upar), fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

# cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    # if tr in ['EGS', 'EGT']:
    #     nt = len(closl[sr][tr]['thour'])
    #     print("nt = %d, closl[%8s][%3s]['thour']" % (nt, sr, tr)) 
    #     print("N triangle closl[sr][tr]['thour']  cl['thour']  " \
    #           "closl[sr][tr]['tau_mbd']  cl[clonm]")
    #     for i in range(nt):
    #           print("%3d %3s  %10.4f  %10.4f  %10.4f  %10.4f" % \
    #                 (i, tr, closl[sr][tr]['thour'][i], cl['thour'][i],
    #                  closl[sr][tr]['tau_mbd'][i], cl[clonm][i]))
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-1,60)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
    axdc.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcc = np.append(tcc, cl[clonm])
    thc = np.append(thc, cl['thour'])
    npt = npt + len(cl['thour'])
axdc.grid(1)

axdc.set_title("Circular", fontsize=10)


axhc = ax1['his_cir']
axhc.hist(tcc, 51, orientation='horizontal', color='r')
axhc.set_title("%d pts" % npt, fontsize=10)
axhc.grid(1)
axhc.set_xlim(-1,60)

pl.show()

if sf: pl.savefig("VO2187_Source_2229+695_%s_Closure.pdf" % upar, format='pdf')


# sys.exit(0)

#
# Plotting for the 1849+670 source
#

sr = '1849+670'

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

ax_col.set_title("VO2187 Source 1849+670 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-1,115)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-1,115)


pl.show()

if sf: pl.savefig("VO2187_Source_1849+670_%s_Closure.pdf" % upar, format='pdf')





#
# Plotting for the 0059+581 source
#

sr = '0059+581'

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

ax_col.set_title("VO2187 Source 0059+581 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-1,70)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-1,70)


pl.show()

if sf: pl.savefig("VO2187_Source_0059+581_%s_Closure.pdf" % upar, format='pdf')







#
# Plotting for the 0613+570 source
#

sr = '0613+570'

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

ax_col.set_title("VO2187 Source 0613+570 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-1,100)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-1,100)


pl.show()

if sf: pl.savefig("VO2187_Source_0613+570_%s_Closure.pdf" % upar, format='pdf')






#
# Plotting for the 3C418 source
#

sr = '3C418'

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

ax_col.set_title("VO2187 Source 3C418 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-0.5,15)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-0.5,15)


pl.show()

if sf: pl.savefig("VO2187_Source_3C418_%s_Closure.pdf" % upar, format='pdf')






#
# Plotting for the 0955+476 source
#

sr = '0955+476'

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

ax_col.set_title("VO2187 Source 0955+476 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-0.5, 35)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-0.5, 35)


pl.show()

if sf: pl.savefig("VO2187_Source_0955+476_%s_Closure.pdf" % upar, format='pdf')






#
# Plotting for the 3C274 source
#

sr = '3C274'

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

ax_col.set_title("VO2187 Source 3C274 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-0.1, 7)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-0.1, 7)


pl.show()

if sf: pl.savefig("VO2187_Source_3C274_%s_Closure.pdf" % upar, format='pdf')





#
# Plotting for the DA426 source
#

sr = 'DA426'

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

ax_col.set_title("VO2187 Source DA426 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
#if clonm == 'cloph': axhl.set_xlim(-1,115)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
#axhc.set_xlim(-1,115)


pl.show()

if sf: pl.savefig("VO2187_Source_DA426_%s_Closure.pdf" % upar, format='pdf')





#
# Plotting for the OJ287 source
#

sr = 'OJ287'

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

ax_col.set_title("VO2187 Source OJ287 %s Closure" % upar, fontsize=12)


hist_colr = 'red'

axdl = ax1['dis_lin']

cl = closl[sr]

npt = 0
tcl = np.empty(0, dtype=float)
thr = np.empty(0, dtype=float)
for tr in closl[sr].keys():
    cl = closl[sr][tr]
    axdl.plot(cl['thour'], cl[clonm], '.', ms=8, color=cols[tr])
    tcl = np.append(tcl, cl[clonm])
    thr = np.append(thr, cl['thour'])
    npt = npt + len(cl['thour'])
axdl.grid(1)

axdl.set_title("Linear", fontsize=10)

# ni, be = np.histogram(tcl, 51)

axhl = ax1['his_lin']
axhl.hist(tcl, 51, orientation='horizontal', color='r')
axhl.grid(1)
axhl.set_title("%d pts" % npt, fontsize=10)
if clonm == 'cloph': axhl.set_xlim(-0.1,8)

#
# Circular
#

axdc = ax1['dis_cir']

npt = 0
tcc = np.empty(0, dtype=float)
thc = np.empty(0, dtype=float)
for tr in closc[sr].keys():
    cl = closc[sr][tr]
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
axhc.set_xlim(-0.1,8)


pl.show()

if sf: pl.savefig("VO2187_Source_OJ287_%s_Closure.pdf" % upar, format='pdf')

















# tr = 'EGH'
# axdl.plot(cl[tr]['thour'], cl[tr][clonm], '.', color=cols[tr])
# tr = 'EGM'
# axdl.plot(cl[tr]['thour'], cl[tr][clonm], '.', color=cols[tr])
# tr = 'GHM'
# axdl.plot(cl[tr]['thour'], cl[tr][clonm], '.', color=cols[tr])
# tr = 'HIS'
# axdl.plot(cl[tr]['thour'], cl[tr][clonm], '.', color=cols[tr])
# tr = 'EMS'
# axdl.plot(cl[tr]['thour'], cl[tr][clonm], '.', color=cols[tr])
# tr = 'EMT'
# axdl.plot(cl[tr]['thour'], cl[tr][clonm], '.', color=cols[tr])
