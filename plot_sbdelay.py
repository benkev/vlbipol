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


fig1 = pl.figure(figsize=(8, 12))
fig2 = pl.figure(figsize=(8, 12))
fig3 = pl.figure(figsize=(8, 12))

#
# Compute and save RMSE and Pearson's correlation coefficients for each baseline
#
# rmse: Root mean square errors (RMSE) between lin pol and cir pol curves
# r_corr: Correlation coefficients  between lin pol and cir pol curves
#
# WRONG:
# rmse_r: RMSE reduced with respect to the mean of the average between
#         lin pol and cir pol curves
#
rmse = np.zeros(nbls, dtype=float)    # Root mean square error (RMSE) for SBD
# rmse_r = np.zeros(nbls, dtype=float)  # RMSE reduced wrt abs of average
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for SBD

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 0   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    sbd_l_us = np.array(idx3819l_1[bl]['I']['sbdelay'])[istart:] # In useconds
    sbd_c_us = np.array(idx3819c_1[bl]['I']['sbdelay'])[istart:] # In useconds

    sbd_l = sbd_l_us*1e6     # Convert us to ps
    sbd_c = sbd_c_us*1e6     # Convert us to ps
    bsbd = sbd_l - sbd_c                 # Bias
    
    sbd0_l = sbd_l - sbd_l.mean()        # Subtract SBD means
    sbd0_c = sbd_c - sbd_c.mean()        # Subtract SBD means
    dsbd = sbd0_l - sbd0_c
    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dsbd**2)/npt)
    r_corr[ibl] = sum(sbd0_l*sbd0_c)/np.sqrt(sum(sbd0_l**2)*sum(sbd0_c**2))

    pl.figure(fig1)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, sbd_l , label='Lin_I, '+bl)
    pl.plot(tim, sbd_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    ax1 = pl.gca()
    
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    pl.figure(fig2)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim,  dsbd, color='orangered', label=bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    ax2 = pl.gca()
    
    pl.ylim(-300, 300)

    pl.text(.03, .92, "r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .80, "RMSE: %.4f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)

    pl.figure(fig3)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, bsbd, color='brown', label=bl)
    pl.grid(True)
    ax3 = pl.gca()

    pl.ylim(-1700, 1400)

    pl.text(.90, .90, bl, transform=ax3.transAxes, fontsize=10)
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax3.transAxes, \
            fontsize=9)
    pl.text(.03, .90, "bias mean: %.1f" % bsbd.mean(), transform=ax3.transAxes,\
            fontsize=9)

    print("%s: dsbd.min = %7.2f, dsbd.max = %7.2f, bsbd.mean = %7.1f \t "
          "rmse = %5.2f, r_corr = %f" % \
          (bl, dsbd.min(), dsbd.max(), bsbd.mean(), rmse[ibl], r_corr[ibl])) 
    
    ibl = ibl + 1
    
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))
fig3.tight_layout(rect=(0,0,1, 0.95))

pl.figure(fig1)
pl.figtext(0.20, 0.96, "Pseudo-Stokes I SBD (ps) vs Time (min), " \
           "Lin & Cir Pol after PolConvert", fontsize=11)

pl.figure(fig2)
pl.figtext(0.08, 0.96, "SBD Residuals (ps) vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means subtracted)", \
           fontsize=11)

pl.figure(fig3)
pl.figtext(0.20, 0.96, "SBD Bias (ps) vs Time (min), " \
           " between Lin & Cir Pol after PolConvert", \
           fontsize=11)

pl.figure(fig1)
pl.savefig("SBD_Lin_I_and_Cir_I.pdf", format='pdf')
pl.figure(fig2)
pl.savefig("SBD_Lin_I_minus_Cir_I.pdf", format='pdf')
pl.figure(fig3)
pl.savefig("SBD_bias_between_Lin_I_and_Cir_I.pdf", format='pdf')


pl.show()





    
