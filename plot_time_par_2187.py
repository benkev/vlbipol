help_text = '''
plot_time_par.py - Plot temporal variations of MBD, SBD, or SNR.
'''
import sys, glob, re

expm = '2187'     # Experiment number: VO2187
snr_floor = 30    # Data with SNR < snr_floor will be discarded.
n_sig = 6         # Data deviating beyond n_sig*std will be discarded.


if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_time_var.py <par> [save], ")
    print("       where <par> is either MBD or SBD or SNR.")
    print("       save (optional): save  figures in pdf format.")
    sys.exit(0)
    
par = sys.argv[1]
par = par.upper()
if par != 'MBD' and par != 'SBD' and par != 'SNR':
    print("Argument can be MBD or SBD or SNR. Entered '%s'. Exiting." %
          sys.argv[1])
    sys.exit(1)

sf = False  # Save figure request
if len(sys.argv) == 3:
    if sys.argv[2] == 'save':
        sf = True
    else:
        print("Argument can only be 'save'. Entered '%s'. Exiting." %
              sys.argv[2])
        sys.exit(1)
    
#
# Find the indices with pseudo-Stokes pol prods (file names with 'lI' and 'cI')
#
fidxl = glob.glob("idx" + expm + "lI.pkl")  # Linear pol with pseudo-Stokes I
fidxc = glob.glob("idx" + expm + "cI.pkl")  # Circular pol with pseudo-Stokes I

no_idxl = False
if fidxl is []:
    no_idxl = True
elif not re.match(r"idx[0-9]{4}lI\.pkl", fidxl[0]):
    no_idxl = True
if no_idxl:
    print("No linear polprod data file idx" + expm + "lI.pkl.  Exiting.")
    sys.exit(0)

no_idxc = False
if fidxc is []:
    no_idxc = True
elif not re.match(r"idx[0-9]{4}cI\.pkl", fidxc[0]):
    no_idxc = True
if no_idxc:
    print("No circular polprod data file idx" + expm + "cI.pkl.  Exiting.")
    sys.exit(0)
    
fidxl = fidxl[0] 
fidxc = fidxc[0] 

    
    
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

#
# Unpickle it:
#
# with open('idx3819l.pkl', 'rb') as finp:
#     idx3819l_1 = pickle.load(finp)

# with open('idx3819cI.pkl', 'rb') as finp:
#     idx3819c_1 = pickle.load(finp)

with open(fidxl, 'rb') as finp:  # Linear pols index file
    idxl = pickle.load(finp)

with open(fidxc, 'rb') as finp:  # Circular pols index file
    idxc = pickle.load(finp)

    

if par == 'MBD' or par == 'SBD':
    ps = "(ps)"
else: # if par == 'SNR':
    ps = ""
    
#
# Circular pol uses FEWER baselines than linear pol.
# 
# This is because PolConvert, by some reason, omitted the 'Y' (i.e. 'Yj')
# station, so it is not present in the baseline list.
#
# Below we read both linear and circular baselines and select only those
# baselines that are present in both cases.
#
bls_l = set(idxl.keys())    # Linear baselines (set)
nbls_l = len(bls_l)
bls_c = set(idxc.keys())    # Circular baselines (set)
nbls_c = len(bls_c)

bls = list(bls_l & bls_c)   # Find the lin and cir sets intersection 
bls.sort()                  # Sort baselines lexicographically

#
# Exclude the 'ST' baseline: the S and T stations are too close to each other
#
if 'ST' in bls:
    iST = bls.index('ST')
    bls.pop(iST)

nbls = len(bls)


# sys.exit(0)


fig1 = pl.figure(figsize=(8, 12))
fig2 = pl.figure(figsize=(8, 12))
fig3 = pl.figure(figsize=(8, 12))

#
# Compute and save RMSE and Pearson's correlation coefficients for each baseline
#
# rmse: Root mean square errors (RMSE) between lin pol and cir pol curves
# r_corr: Correlation coefficients  between lin pol and cir pol curves
#
rmse = np.zeros(nbls, dtype=float)    # Root mean square error (RMSE) for MBD
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for MBD

#
# Table header
#
if par == 'SNR':
    print("BL  avg SNR    rmse   relerr    avg bias   r_corr")
else: # if par == 'MBD' or par == 'SBD': 
    print("BL  avg SNR  avg %s     rmse   relerr,%%   avg bias   r_corr" % par)
#
# Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or 'snr'
#
if par == 'MBD':
    parname = 'mbdelay'
elif par == 'SBD':
    parname = 'sbdelay'
else:
    parname = 'snr'
    
    
ibl = 0   # Baseline
for bl in bls:   # Loop over the baselines
    snr_l = np.array(idxl[bl]['I']['snr'])
    snr_c = np.array(idxc[bl]['I']['snr'])
    tim = np.array(idxl[bl]['I']['time'])
    tim0 = tim[0]            # Save the starttime in case snr[0] <= snr_floor

    #
    # Exclude data with SNR <= snr_floor
    #
    isnr_floorl = np.where(snr_l <= snr_floor)[0]
    isnr_floorc = np.where(snr_c <= snr_floor)[0]
    # Merge lin & cir indices where snr <= snr_floor
    isnr_floor = np.concatenate((isnr_floorl, isnr_floorc))
    isnr_floor = np.unique(isnr_floor)
    isnr_floor.sort()                 # Sort the array in-place

    tim = np.delete(tim, isnr_floor) # Exclude data with SNR <= 30
    snr_l = np.delete(snr_l, isnr_floor) # Exclude data with SNR <= 30
    snr_c = np.delete(snr_c, isnr_floor) # Exclude data with SNR <= 30

    tim = tim - tim0    # Count time in minutes from the session start
    tim = tim / 3600    # Change timestamps from seconds to hours
   
    # Average of Lin and Cir means of SNR
    snr_a = (abs(snr_l.mean()) + abs(snr_c.mean()))/2
    
    if par == 'MBD' or par == 'SBD':
        par_l_us = np.array(idxl[bl]['I'][parname]) # In useconds
        par_c_us = np.array(idxc[bl]['I'][parname]) # In useconds
        par_l_us =  np.delete(par_l_us, isnr_floor)  # Excl data with SNR <= 30
        par_c_us =  np.delete(par_c_us, isnr_floor)  # Excl data with SNR <= 30
        par_l = par_l_us*1e6           # Convert us to ps
        par_c = par_c_us*1e6           # Convert us to ps
    else: # if par == 'SNR':
        par_l = np.copy(snr_l)
        par_c = np.copy(snr_c)
   
    #
    # Exclude points beyond +-n_sig*std
    #
    par0_l = par_l - par_l.mean()        # Subtract MBD means, lin pol
    par0_c = par_c - par_c.mean()        # Subtract MBD means, cir pol

    sl = n_sig*np.std(par0_l)
    sc = n_sig*np.std(par0_c)
    isl = np.where((par0_l > sl) | (par0_l < -sl))[0]
    isc = np.where((par0_c > sc) | (par0_c < -sc))[0]    
    isg = np.concatenate((isl, isc))
    isg = np.unique(isg)
    isg.sort()                 # Sort the array in-place

    tim = np.delete(tim, isg)       # Exclude data with |param| > n_sig*std 
    par0_l = np.delete(par0_l, isg) # Exclude data with |param| > n_sig*std 
    par0_c = np.delete(par0_c, isg) # Exclude data with |param| > n_sig*std 
    par_l = np.delete(par_l, isg)   # Exclude data with |param| > n_sig*std 
    par_c = np.delete(par_c, isg)   # Exclude data with |param| > n_sig*std 

    bpar = par_l - par_c                 # Bias between Lin and Cir
    
    # Average of Lin and Cir means of the parameter (MBD or SDD)
    par_a = (abs(par_l.mean()) + abs(par_c.mean()))/2 # Avg Lin and Cir means
    
    dpar = bpar - bpar.mean()                   # Residuals

    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dpar**2)/npt)
    r_corr[ibl] = sum(par0_l*par0_c)/np.sqrt(sum(par0_l**2)*sum(par0_c**2))


    # sys.exit(0)


    
    pl.figure(fig1)
    pl.subplot(7, 3, ibl+1)
    # pl.plot(tim, par_l, label='Lin_I, mean: %.1f' % par_l.mean())
    # pl.plot(tim, par_c, label='Cir_I, mean: %.1f' % par_c.mean())
    # pl.plot(tim, par_l, 'b', label='Lin_I, mean: %.1f' % par_l.mean())
    pl.plot(tim, par_l, 'b.', markersize=2)
    #pl.plot(tim, par_c, 'g', label='Cir_I, mean: %.1f' % par_c.mean())
    pl.plot(tim, par_c, 'g.', markersize=2)
    
    pl.grid(True)
    #pl.figtext("minutes")
    #pl.ylabel("ps")
    # pl.legend(loc='upper left', prop={'size': 9})
    ax1 = pl.gca()

    # pl.plot(tim[isl], par_l[isl], 'c.', markersize=8)
    # pl.plot(tim[isc], par_c[isc], 'm.', markersize=8)

    # if par == 'SNR':
    #     pl.ylim(0, 6000)
    
    pl.text(.88, .02, bl, transform=ax1.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .04, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    # sys.exit(0)

    
    pl.figure(fig2)
    pl.subplot(7, 3, ibl+1)
    # pl.plot(tim, dpar, color='red')
    pl.plot(tim, dpar, '.', markersize=2, color='red')
    pl.grid(True)
    ax2 = pl.gca()

    # if  par == 'MBD':
    #     #pl.ylim(-25, 25)
    #     pl.ylim(-1000, 1000)
    # elif par == 'SBD':
    #     pl.ylim(-300, 300)
    # else: # if par == 'SNR':
    #     pl.ylim(-100, 100)

    pl.text(.88, 1.01, bl, transform=ax2.transAxes, fontsize=10, weight="bold")
    # pl.text(.03, .92, r"r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes,\
    #         fontsize=9)
    pl.text(.03, 1.01, r"RMSE: %.2f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .82, r"$\overline{\mathrm{SNR}}$:  %.1f" % \
            snr_a, transform=ax2.transAxes, fontsize=9)

    pl.figure(fig3)
    pl.subplot(7, 3, ibl+1)
    #pl.plot(tim, bpar, color='brown')
    pl.plot(tim, bpar, '.', markersize=3, color='brown')
    pl.grid(True)
    ax3 = pl.gca()

    # if  par == 'MBD':
    #     #pl.ylim(-200, 250)
    #     pl.ylim(-1000, 1000)
    # elif par == 'SBD':
    #     pl.ylim(-1700, 1400)
    # else: # if par == 'SNR':
    #     pl.ylim(-150, 350)

    pl.text(.88, 1.01, bl, transform=ax3.transAxes, fontsize=10, weight="bold")
    # pl.text(.84, .89, bl, transform=ax3.transAxes, fontsize=10, weight="bold")
    # pl.text(.03, .02, r"r_corr: %.6f" % r_corr[ibl], transform=ax3.transAxes,
    #         fontsize=9)
    #pl.text(.03, .84, r"$\overline{\mathrm{bias}}$:  %.1f" % \
    pl.text(.03, 1.01, r"$\overline{\mathrm{bias}}$:  %.1f" % \
            bpar.mean(), transform=ax3.transAxes, fontsize=9)

    #
    # Table
    #
    rel_err = 100*abs(rmse[ibl]/par_a) # Relative error
    if par == 'SNR':
        print("%s  %7.1f   %5.1f   %5.2f     %5.1f     %8.6f" % \
              (bl, par_a, rmse[ibl], rel_err, bpar.mean(), r_corr[ibl])) 
    else: # if par == 'MBD' or par == 'SBD': 
        print("%s  %7.1f  %7.1f  %7.1f  %7.2f    %7.1f    %8.6f" % \
              (bl, snr_a, par_a, rmse[ibl], rel_err, bpar.mean(), r_corr[ibl])) 

    ibl = ibl + 1
    
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))
fig3.tight_layout(rect=(0,0,1, 0.95))

pl.figure(fig1)
pl.figtext(0.03, 0.96, "VO%s: Pseudo-Stokes I %s %s vs Time (hours), " \
           "Lin (blue) & Cir (green) Pol after PolConvert" % (expm, par, ps), \
           fontsize=11)
pl.figtext(0.75, 0.10, "SNR > %d" % snr_floor, fontsize=16)
pl.figtext(0.75, 0.07, r"|%s| < %d$\sigma$" % (par, n_sig), fontsize=16)


pl.figure(fig2)
pl.figtext(0.08, 0.96, "VO%s: %s Residuals %s vs Time (hours), " \
           " between Lin & Cir Pol after PolConvert" \
           % (expm, par, ps), fontsize=11)
pl.figtext(0.75, 0.10, "SNR > %d" % snr_floor, fontsize=16)
pl.figtext(0.75, 0.07, r"|%s| < %d$\sigma$" % (par, n_sig), fontsize=16)


pl.figure(fig3)
pl.figtext(0.08, 0.96, "VO%s: %s Bias %s vs Time (hours), " \
           " between Lin & Cir Pol after PolConvert" % (expm, par, ps), \
           fontsize=11)
pl.figtext(0.75, 0.10, r"SNR > %d" % snr_floor, fontsize=16)
pl.figtext(0.75, 0.07, r"|%s| < %d$\sigma$" % (par, n_sig), fontsize=16)

if sf:
    pl.figure(fig1)
    pl.savefig("VO%s_%s_Lin_I_and_Cir_I_SNR_floor_%d.pdf" % \
               (expm, par, snr_floor), format='pdf')
    pl.figure(fig2)
    pl.savefig("VO%s_%s_Lin_I_minus_Cir_I_SNR_floor_%d.pdf" % \
               (expm, par, snr_floor), format='pdf')
    pl.figure(fig3)
    pl.savefig("VO%s_%s_bias_between_Lin_I_and_Cir_I_SNR_floor_%d.pdf" % \
               (expm, par, snr_floor), format='pdf')

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
#




    
