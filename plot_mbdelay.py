import sys
import pickle
import numpy as np
import matplotlib.pyplot as pl

# from vpal.utility import int_to_time, time_to_int

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
print("pl.isinteractive() -> ", pl.isinteractive())

#
# Unpickle it:
#
with open('idx3819l.pkl', 'rb') as finp:
    idx3819l_1 = pickle.load(finp)

# with open('idx3819c.pkl', 'rb') as finp:
#     idx3819c_1 = pickle.load(finp)

with open('idx3819cI.pkl', 'rb') as finp:
    idx3819c_1 = pickle.load(finp)


bls = list(idx3819l_1.keys())   # Baselines
bls.sort()                      # Lexigraphically sorted baselines
nbls = len(bls)

# lmin = np.zeros(nbls, dtype=float)
# lmax = np.zeros(nbls, dtype=float)
# cmin = np.zeros(nbls, dtype=float)
# cmax = np.zeros(nbls, dtype=float)
# dmin = np.zeros(nbls, dtype=float)
# dmax = np.zeros(nbls, dtype=float)

# #
# # Find vertical ranges for plotting
# #
# ibl = 0
# for bl in bls:   # Loop over the baselines
#     mbd_l = np.array(idx3819l_1[bl]['I']['mbdelay'])[:-25]
#     mbd_c = np.array(idx3819c_1[bl]['I']['mbdelay'])[:-25]
#     lmin[ibl] = mbd_l.min()
#     lmax[ibl] = mbd_l.max()
#     cmin[ibl] = mbd_c.min()
#     cmax[ibl] = mbd_c.max()

#     dmbd = mbd_c - mbd_l
#     dmin[ibl] = dmbd.min()
#     dmax[ibl] = dmbd.max()

    
#     ibl = ibl + 1

# ymin = min(lmin.min(), cmin.min())
# ymax = max(lmax.max(), cmax.max())

# ydmin = dmin.min()
# ydmax = dmax.max()

   
# # sys.exit(0)



fig1 = pl.figure(1, figsize=(8, 12))
fig2 = pl.figure(2, figsize=(8, 12))

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 1   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    mbd_l = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:]
    mbd_c = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:]

    pl.figure(1)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, mbd_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd_c , label='Cir_I, '+bl)
    # pl.ylim(ymin,ymax)
    pl.legend(loc='upper right')

    pl.figure(2)
    pl.subplot(5, 3, ibl)
    pl.plot(tim,  abs(mbd_c - mbd_l), color='orangered', label=bl)
    # pl.ylim(ydmin,ydmax)
    pl.legend(loc='upper right')

    ibl = ibl + 1
    
# fig.tight_layout()
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))

pl.figure(1)
pl.figtext(0.02, 0.97, "3819 Pseudo-Stokes I MultiBand Delays vs Time (min), " \
           "Linear & Circular Polarization after PolConvert", fontsize=11)

pl.figure(2)
pl.figtext(0.02, 0.97, "3819 MultiBand Delay Differences vs Time (min), " \
           " between Linear & Circular Pol after PolConvert", fontsize=11)



pl.show()





    
