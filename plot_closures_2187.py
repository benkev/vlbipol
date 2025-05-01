'''
plot_closures_2187.py - plot closure phase and closure delay for residual and
                             total MBD or SBD, session 2187 (VO2187).
'''

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

parlist = list(arg_to_par.keys())
pararg = (sys.argv[1]).lower()
if pararg not in arg_to_par.keys():
    print("Argument can be one of %r. Entered '%s'. Exiting." %
          (parlist, sys.argv[1]))
    sys.exit(0)

sf = False  # Save figure request
if len(sys.argv) == 3:
    if sys.argv[2] == 'save':
        sf = True
    else:
        print("Argument can only be 'save'. Entered '%s'. Exiting." %
              sys.argv[2])
        sys.exit(0)
    
parname = arg_to_par[pararg]           # Like 'mbdelay' for 'mbd' etc.

import matplotlib
matplotlib.use('qtagg')  # Force the interactive backend
import pickle
import copy
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib.patches as patches
from itertools import combinations
from bisect import bisect_right  # Bisection algorithm for efficient search
from group_tails import find_tail_bounds, group_tails
import libvp
import libplt

# pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.


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

bls = list(idxc.keys())    # List of baselines 
bls.sort()                 # Sort baselines lexicographically

#
# Exclude the 'ST' baseline: the S and T stations are too close to each other
#
if 'ST' in bls:
    iST = bls.index('ST')
    bls.pop(iST)

nbls = len(bls)

trians = libvp.find_baseline_triangles(bls)
ntri = len(trians)

cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))


upar = parname.upper()

if pararg == "mbd":
    yl = 200
elif pararg == "sbd":  
    yl = 2000
    
dist_ylim = (-yl,yl)
#his_xlim = (-50,6000)   # SBD


gs_kw1 = dict(width_ratios=[0.75, 0.25], height_ratios=[0.15, 0.425, 0.425])
fig1, ax1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                               ['dis_lin', 'his_lin'],
                               ['dis_cir', 'his_cir']],
                               gridspec_kw=gs_kw1, figsize=(8.4, 10),
                               layout="constrained")
# sys.exit(0)

# #
# # Plot color legend on top
# #
# ax_col = ax1['col_legend']    # Get the axis for color legend

# #libplt.plot_closure_legend(ax_col, trians, cols, upar, fs=9)  

# hist_colr = 'red'

# ax_lin = ax1['dis_lin'] # Plot distr of FourFit ps.-I param vs Time  
# ttl_lin = "Fourfit Pseudo-I, %s vs Time (%d triangles)" # % (upar, ntri)




# #pl.show()

# sys.exit(0)


# plot_closures_distr(ax_ffd, tau_l, cols, parname, ttl_ffd, yl=dist_ylim)


# ax_ffh = ax1['his_lin']  # Plot hist of FourFit I param
# ttl_ffh = "%s (%d points)" # % (upar, nfinite)

# # nclod_l = plot_closures_hist_horiz(ax_ffh, tau_l, parname, hist_colr, ttl_ffh,
# #                                  yl=dist_ylim, xl=hist_xlim)
# # nclod_l = plot_closures_hist_horiz(ax_ffh, tau_l, parname, hist_colr, ttl_ffh)


# clod0 = np.zeros((0,), dtype=np.float64)
# for tr in tau_l.keys():
#     clod0 = np.append(clod0, tau_l[tr]['tau'])
    
# ixc = np.where(abs(clod0) < yl)
# clod = np.copy(clod0[ixc])
# ax_ffh.hist(clod, 101, color=hist_colr, orientation='horizontal')
# ax_ffh.set_ylim(dist_ylim)
# # ax_ffh.set_ylim(hist_xlim)


# ax_pcd = ax1['dis_cir'] # Plot distr of PolConvert I param vs Time
# ttl_pcd = "PolConvert I, %s vs Time (%d triangles)" # % (upar, ntri)

# plot_closures_distr(ax_pcd, tau_c, cols, parname, ttl_pcd, yl=dist_ylim)

# ax_pch = ax1['his_cir']  # Plot hist of PolConvert I param
# ttl_pch = "%s (%d points)" # % (upar, nfinite)

# # nclod_c = plot_closures_hist_horiz(ax_pch, tau_c, parname, hist_colr, ttl_pch,
# #                          yl=dist_ylim, xl=hist_xlim)
# #nclod_c = plot_closures_hist_horiz(ax_pch, tau_c, parname, hist_colr, ttl_pch)

# clod0 = np.zeros((0,), dtype=np.float64)
# for tr in tau_c.keys():
#     clod0 = np.append(clod0, tau_c[tr]['tau'])
    
# ixc = np.where(abs(clod0) < yl)
# clod = np.copy(clod0[ixc])
# ax_pch.hist(clod, 101, color=hist_colr, orientation='horizontal')
# ax_pch.set_ylim(dist_ylim)
# # ax_pch.set_ylim(hist_xlim)





