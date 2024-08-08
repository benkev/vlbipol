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
#     sbd_l = np.array(idx3819l_1[bl]['I']['sbdelay'])[:-25]
#     sbd_c = np.array(idx3819c_1[bl]['I']['sbdelay'])[:-25]
#     lmin[ibl] = sbd_l.min()
#     lmax[ibl] = sbd_l.max()
#     cmin[ibl] = sbd_c.min()
#     cmax[ibl] = sbd_c.max()

#     dsbd = sbd_c - sbd_l
#     dmin[ibl] = dsbd.min()
#     dmax[ibl] = dsbd.max()

    
#     ibl = ibl + 1

# ymin = min(lmin.min(), cmin.min())
# ymax = max(lmax.max(), cmax.max())

# ydmin = dmin.min()
# ydmax = dmax.max()

   
# # sys.exit(0)



fig21 = pl.figure(21, figsize=(8, 12))
fig22 = pl.figure(22, figsize=(8, 12))

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 1   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    sbd_l_us = np.array(idx3819l_1[bl]['I']['sbdelay'])[istart:] # In useconds
    sbd_c_us = np.array(idx3819c_1[bl]['I']['sbdelay'])[istart:] # In useconds

    sbd_l = sbd_l_us*1e6     # Convert us to ps
    sbd_c = sbd_c_us*1e6     # Convert us to ps
    
    pl.figure(21)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, sbd_l , label='Lin_I, '+bl)
    pl.plot(tim, sbd_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    

    sbd0_l = sbd_l - sbd_l.mean()        # Subtract SBD means
    sbd0_c = sbd_c - sbd_c.mean()        # Subtract SBD means

    dsbd_cl = sbd0_c - sbd0_l
    
    pl.figure(22)
    pl.subplot(5, 3, ibl)
    pl.plot(tim,  dsbd_cl, color='orangered', label=bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    
    pl.ylim(-300, 300)

    print(bl, ": cir-lin min and max: ", dsbd_cl.min(), dsbd_cl.max()) 
    
    ibl = ibl + 1
    
# fig.tight_layout()
fig21.tight_layout(rect=(0,0,1, 0.95))
fig22.tight_layout(rect=(0,0,1, 0.95))

pl.figure(21)
pl.figtext(0.15, 0.97, "Pseudo-Stokes I SBD (ps) vs Time (min), " \
           "Lin & Cir Pol after PolConvert", fontsize=11)

pl.figure(22)
pl.figtext(0.05, 0.97, "SBD Differences (ps) vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means suntracted)", \
           fontsize=11)

pl.figure(21)
pl.savefig("SBD_Lin_I_and_Cir_I.eps", format='eps')
pl.figure(22)
pl.savefig("SBD_Lin_I_minus_Cir_I.eps", format='eps')


pl.show()





    
