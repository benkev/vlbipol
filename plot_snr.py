import sys
import pickle
import numpy as np
import matplotlib.pyplot as pl

# from vpal.utility import int_to_time, time_to_int

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

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
rmse = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for SNR
# rmse_r = np.zeros(nbls, dtype=float)  # RMSE reduced wrt abs of average
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for SNR


#
# Table header
#
print("BL  avg SNR   rmse  relerr,%   avg bias    r_corr") 

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
    bsnr = snr_l - snr_c                 # Bias
    snr_a = (abs(snr_l.mean()) + abs(snr_c.mean()))/2 # Avg Lin and Cir means

    snr0_l = snr_l - snr_l.mean()        # Subtract SNR means
    snr0_c = snr_c - snr_c.mean()        # Subtract SNR means
    dsnr = snr0_l - snr0_c
    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dsnr**2)/npt)
    r_corr[ibl] = sum(snr0_l*snr0_c)/np.sqrt(sum(snr0_l**2)*sum(snr0_c**2))


    
    pl.figure(fig1)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, snr_l , label='Lin_I, mean: %.1f' % snr_l.mean())
    pl.plot(tim, snr_c , label='Cir_I, mean: %.1f' % snr_c.mean())
    pl.grid(True)
    pl.legend(loc='upper left', prop={'size': 9})
    ax1 = pl.gca()
        
    pl.ylim(0, 6000)
   
    pl.text(.88, .02, bl, transform=ax1.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    pl.figure(fig2)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim,  dsnr, color='red')
    pl.grid(True)
    ax2 = pl.gca()
    
    pl.ylim(-100, 100)

    pl.text(.88, .90, bl, transform=ax2.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .92, "r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .80, "RMSE: %.4f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)

    pl.figure(fig3)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, bsnr, color='brown')
    pl.grid(True)
    ax3 = pl.gca()

    pl.ylim(-150, 350)

    pl.text(.88, .90, bl, transform=ax3.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax3.transAxes, \
            fontsize=9)
    pl.text(.03, .90, "bias mean: %.1f" % bsnr.mean(), transform=ax3.transAxes,\
            fontsize=9)

    #
    # Table
    #
    rel_err = 100*abs(rmse[ibl]/snr_a)
    print("%s  %7.1f   %4.1f   %4.2f     %6.1f    %8.6f" % \
          (bl, snr_a, rmse[ibl], rel_err, bsnr.mean(), r_corr[ibl])) 
    
    
    ibl = ibl + 1

    
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))
fig3.tight_layout(rect=(0,0,1, 0.95))

pl.figure(fig1)
pl.figtext(0.20, 0.96, "Pseudo-Stokes I SNR vs Time (min), " \
           "Lin & Cir Pol after PolConvert", fontsize=11)

pl.figure(fig2)
pl.figtext(0.08, 0.96, "SNR Residuals vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means subtracted)", \
           fontsize=11)

pl.figure(fig3)
pl.figtext(0.20, 0.96, "SNR Bias vs Time (min), " \
           " between Lin & Cir Pol after PolConvert", \
           fontsize=11)

pl.figure(fig1)
pl.savefig("SNR_Lin_I_and_Cir_I.pdf", format='pdf')
pl.figure(fig2)
pl.savefig("SNR_Lin_I_minus_Cir_I.pdf", format='pdf')
pl.figure(fig3)
pl.savefig("SNR_bias_between_Lin_I_and_Cir_I.pdf", format='pdf')


pl.show()





    
