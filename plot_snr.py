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


fig31 = pl.figure(31, figsize=(8, 12))
fig32 = pl.figure(32, figsize=(8, 12))

#
# Compute and save RMSE and Pearson's correlation coefficients for each baseline
#
# rmse: Root mean square errors (RMSE) between lin pol and cir pol curves
# rmse_r: RMSE reduced with respect to the mean of the average between
#         lin pol and cir pol curves
# r_corr: Correlation coefficients  between lin pol and cir pol curves
#
rmse = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for SNR
rmse_r = np.zeros(nbls, dtype=float)  # RMSE reduced wrt abs of average
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for SNR


#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 0   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    snr_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:] # In useconds
    snr_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:] # In useconds
    snr_a = (snr_l + snr_c)/2      # Average of the lin and cir curves

    snr0_l = snr_l - snr_l.mean()        # Subtract SNR means
    snr0_c = snr_c - snr_c.mean()        # Subtract SNR means
    dsnr = snr0_l - snr0_c
    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dsnr**2)/npt)
    rmse_r[ibl] = rmse[ibl]/abs(snr_a.mean()) # RMSE reduced wrt abs of average
    r_corr[ibl] = sum(snr0_l*snr0_c)/np.sqrt(sum(snr0_l**2)*sum(snr0_c**2))


    
    pl.figure(31)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, snr_l , label='Lin_I, '+bl)
    pl.plot(tim, snr_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    ax1 = pl.gca()
        
    pl.ylim(0, 6000)
   
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    pl.figure(32)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim,  dsnr, color='orangered', label=bl)
    pl.grid(True)
    pl.legend(loc='upper right')
    ax2 = pl.gca()
    
    pl.ylim(-100, 100)

    pl.text(.03, .92, "r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .80, "RMSE: %.4f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .68, "RMSE_r: %.5f" % rmse_r[ibl], transform=ax2.transAxes, \
            fontsize=9)

    print("%s: dsnr.min = %.3f, dsnr.max = %.3f, \t rmse = %.4f, "\
          "rmse_r = %8f, r_corr = %f" % \
          (bl, dsnr.min(), dsnr.max(), rmse[ibl], rmse_r[ibl], r_corr[ibl])) 
    
    ibl = ibl + 1
    
# fig.tight_layout()
fig31.tight_layout(rect=(0,0,1, 0.95))
fig32.tight_layout(rect=(0,0,1, 0.95))

pl.figure(31)
pl.figtext(0.20, 0.96, "Pseudo-Stokes I SNR vs Time (min), " \
           "Lin & Cir Pol after PolConvert", fontsize=11)

pl.figure(32)
pl.figtext(0.08, 0.96, "SNR Residuals vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means subtracted)", \
           fontsize=11)

pl.figure(31)
pl.savefig("SNR_Lin_I_and_Cir_I.eps", format='eps')
pl.figure(32)
pl.savefig("SNR_Lin_I_minus_Cir_I.eps", format='eps')




pl.show()





    
