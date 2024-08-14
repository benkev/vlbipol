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


fig11 = pl.figure(11, figsize=(8, 12))
fig12 = pl.figure(12, figsize=(8, 12))

#
# Compute and save RMSE and Pearson's correlation coefficients for each baseline
#
# rmse: Root mean square errors (RMSE) between lin pol and cit pol curves
# r_corr: Correlation coefficients  between lin pol and cit pol curves
#
rmse = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for MBD
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for MBD


#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 0   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    mbd_l_us = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:] # In useconds
    mbd_c_us = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:] # In useconds

    mbd_l = mbd_l_us*1e6     # Convert us to ps
    mbd_c = mbd_c_us*1e6     # Convert us to ps
    
    mbd0_l = mbd_l - mbd_l.mean()        # Subtract MBD means, lin pol
    mbd0_c = mbd_c - mbd_c.mean()        # Subtract MBD means, cir pol
    dmbd = mbd0_l - mbd0_c
    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dmbd**2)/npt)
    r_corr[ibl] = sum(mbd0_l*mbd0_c)/np.sqrt(sum(mbd0_l**2)*sum(mbd0_c**2))

    pl.figure(11)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, mbd_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    ax1 = pl.gca()
    
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    pl.figure(12)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, dmbd, color='orangered', label=bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    ax2 = pl.gca()

    pl.ylim(-25, 25)

    pl.text(.03, .92, "r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .80, "RMSE: %.4f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)

    print("%s: dmbd.min = %.3f, dmbd.max = %.3f, \t rmse = %f, r_corr = %f" % \
          (bl, dmbd.min(), dmbd.max(), rmse[ibl], r_corr[ibl])) 
    
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





    
