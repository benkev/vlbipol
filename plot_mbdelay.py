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



fig11 = pl.figure(11, figsize=(8, 12))
fig12 = pl.figure(12, figsize=(8, 12))

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 1   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    mbd_l_us = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:] # In useconds
    mbd_c_us = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:] # In useconds

    mbd_l = mbd_l_us*1e6     # Convert us to ps
    mbd_c = mbd_c_us*1e6     # Convert us to ps
    
    pl.figure(11)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, mbd_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')


    mbd0_l = mbd_l - mbd_l.mean()        # Subtract MBD means
    mbd0_c = mbd_c - mbd_c.mean()        # Subtract MBD means

    dmbd_lc = mbd0_l - mbd0_c
    
    pl.figure(12)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, dmbd_lc, color='orangered', label=bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    pl.ylim(-25, 25)

    print(bl, ": cir-lin min and max: ", dmbd_lc.min(), dmbd_lc.max()) 
    
    ibl = ibl + 1
    
# fig.tight_layout()
fig11.tight_layout(rect=(0,0,1, 0.95))
fig12.tight_layout(rect=(0,0,1, 0.95))

pl.figure(11)
pl.figtext(0.20, 0.97, "Pseudo-Stokes I MBD (ps) vs Time (min), " \
           "Lin & Cir Pol after PolConvert", fontsize=11)

pl.figure(12)
pl.figtext(0.05, 0.97, "MBD Differences (ps) vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means subtracted)", \
           fontsize=11)


pl.figure(11)
pl.savefig("MBD_Lin_I_and_Cir_I.eps", format='eps')
pl.figure(12)
pl.savefig("MBD_Lin_I_minus_Cir_I.eps", format='eps')

pl.show()





    
