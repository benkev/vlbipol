import sys
import pickle
import numpy as np
import matplotlib.pyplot as pl

# from vpal.utility import int_to_time, time_to_int

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

#
# Control how plain print() behaves without formatting
#
np.set_printoptions(suppress=True, precision=1)
#
# suppress : bool, optional
#     If True, always print floating point numbers using fixed point
#     notation, in which case numbers equal to zero in the current precision
#     will print as zero.  If False, then scientific notation is used when
#     absolute value of the smallest number is < 1e-4 or the ratio of the
#     maximum absolute value to the minimum is > 1e3. The default is False.
# precision : int or None, optional
#     Number of digits of precision for floating point output (default 8).
#     May be None if `floatmode` is not `fixed`, to print as many digits as
#     necessary to uniquely specify the value.


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
# rmse_r: RMSE reduced with respect to the absolute value of the mean
#         of the average between lin pol and cir pol curves
#
rmse = np.zeros(nbls, dtype=float)    # Root mean square error (RMSE) for MBD
# rmse_r = np.zeros(nbls, dtype=float)  # RMSE reduced wrt abs of average
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for MBD

# mbd_l_mean = np.zeros(nbls, dtype=float)
# mbd_c_mean = np.zeros(nbls, dtype=float)

#
# Table header
#
print("BL  avg MBD   rmse  relerr,%   avg bias   r_corr") 

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

    mbd_l = mbd_l_us*1e6           # Convert us to ps
    mbd_c = mbd_c_us*1e6           # Convert us to ps
    # mbd_l_mean[ibl] = mbd_l.mean()
    # mbd_c_mean[ibl] = mbd_c.mean()
    bmbd = mbd_l - mbd_c                 # Bias
    mbd_a = (abs(mbd_l.mean()) + abs(mbd_c.mean()))/2 # Avg Lin and Cir means
    
    mbd0_l = mbd_l - mbd_l.mean()        # Subtract MBD means, lin pol
    mbd0_c = mbd_c - mbd_c.mean()        # Subtract MBD means, cir pol
    dmbd = mbd0_l - mbd0_c               # Residuals

    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dmbd**2)/npt)
#    rmse_r[ibl] = rmse[ibl]/abs(mbd_a.mean()) # RMSE reduced wrt abs of average
    r_corr[ibl] = sum(mbd0_l*mbd0_c)/np.sqrt(sum(mbd0_l**2)*sum(mbd0_c**2))

    pl.figure(fig1)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, mbd_l , label='Lin_I, mean: %.1f' % mbd_l.mean())
    pl.plot(tim, mbd_c , label='Cir_I, mean: %.1f' % mbd_c.mean())
    pl.grid(True)
    pl.legend(loc='upper left', prop={'size': 9})
    ax1 = pl.gca()
    
    pl.text(.88, .02, bl, transform=ax1.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    pl.figure(fig2)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, dmbd, color='red')
    pl.grid(True)
    ax2 = pl.gca()

    pl.ylim(-25, 25)

    pl.text(.88, .90, bl, transform=ax2.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .92, "r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .80, "RMSE: %.4f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)

    pl.figure(fig3)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, bmbd, color='brown') #, label=bl)
    pl.grid(True)
    ax3 = pl.gca()

    pl.ylim(-200, 250)

    pl.text(.88, .90, bl, transform=ax3.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax3.transAxes, \
            fontsize=9)
    pl.text(.03, .90, "bias mean: %.1f" % bmbd.mean(), transform=ax3.transAxes,\
            fontsize=9)

    #
    # Table
    #
    rel_err = 100*abs(rmse[ibl]/mbd_a)
    print("%s  %7.1f   %4.2f    %4.2f     %6.1f     %8.6f" % \
          (bl, mbd_a, rmse[ibl], rel_err, bmbd.mean(), r_corr[ibl])) 
    
    ibl = ibl + 1
    
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))
fig3.tight_layout(rect=(0,0,1, 0.95))

pl.figure(fig1)
pl.figtext(0.20, 0.96, "Pseudo-Stokes I MBD (ps) vs Time (min), " \
           "Lin & Cir Pol after PolConvert", fontsize=11)

pl.figure(fig2)
pl.figtext(0.08, 0.96, "MBD Residuals (ps) vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means subtracted)", \
           fontsize=11)

pl.figure(fig3)
pl.figtext(0.20, 0.96, "MBD Bias (ps) vs Time (min), " \
           " between Lin & Cir Pol after PolConvert", \
           fontsize=11)


pl.figure(fig1)
pl.savefig("MBD_Lin_I_and_Cir_I.pdf", format='pdf')
pl.figure(fig2)
pl.savefig("MBD_Lin_I_minus_Cir_I.pdf", format='pdf')
pl.figure(fig3)
pl.savefig("MBD_bias_between_Lin_I_and_Cir_I.pdf", format='pdf')

pl.show()

#
# Reset to default the control how plain print() behaves without formatting
#
np.set_printoptions(suppress=False, precision=8)
#
# suppress : bool, optional
#     If True, always print floating point numbers using fixed point
#     notation, in which case numbers equal to zero in the current precision
#     will print as zero.  If False, then scientific notation is used when
#     absolute value of the smallest number is < 1e-4 or the ratio of the
#     maximum absolute value to the minimum is > 1e3. The default is False.
# precision : int or None, optional
#     Number of digits of precision for floating point output (default 8).
#     May be None if `floatmode` is not `fixed`, to print as many digits as
#     necessary to uniquely specify the value.





    
