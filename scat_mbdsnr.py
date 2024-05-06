help_text = '''
scat_mbdsnr.py: Scatterplots of the linear Stokes-I MBD vs the polconvert
                circular Stokes-I MBD to visualize the differences.
                Scatterplots of the respective SNRs.
'''
import sys
import pickle
import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl

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
# np.random.seed(19680801)
# cols = np.random.rand(nbls)  # Random colors for each baseline

# Colors for each baseline
cols = ['r', 'g', 'b', 'c', 'm', 'yellow', 'gray', 'orange', 'brown', 'olive',
        'r', 'g', 'b', 'c', 'm']
mar = ['o', 's', '^', '+', '.']

fig1 = pl.figure(1, figsize=(8, 8))
fig2 = pl.figure(2, figsize=(8, 8))

#find 
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 0   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    mbd_l = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:]
    mbd_c = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:]
    snr_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]

    ic = ibl % 5
    im = ibl // 5
    
    pl.figure(1)

    #pl.scatter(mbd_l, mbd_c, s=10, c=col) # s=area, c=colors, alpha=0.5)
    pl.plot(mbd_l, mbd_c, mar[im], color=cols[ic], markersize=5, alpha=0.4, 
            label=bl)
    pl.grid(1)
    pl.xlabel("MBD I (us), Lin. Pol.") 
    pl.ylabel("MBD I (us), Circ. Pol.") 
    pl.legend(loc='upper left')

    pl.figure(2)

    #pl.scatter(snr_l, snr_c, s=10, c=col)  # s=area, c=colors, alpha=0.5)
    pl.plot(snr_l, snr_c, mar[im], color=cols[ic], markersize=5, alpha=0.4,
            label=bl)
    pl.xlabel("SNR I, Lin. Pol.") 
    pl.ylabel("SNR I, Circ. Pol.") 
    pl.grid(1)
    pl.legend(loc='upper left')

    ibl = ibl + 1
    
# fig.tight_layout()
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))

pl.figure(1)
pl.figtext(0.2, 0.97, "ER2201 Linear Stokes-I MBD vs Polconvert "
                "Circular Stokes-I MBD ", fontsize=11)

pl.figure(2)
pl.figtext(0.2, 0.97, "ER2201 Linear Stokes-I SNR vs Polconvert "
                "Circular Stokes-I SNR ", fontsize=11)


pl.show()





    
