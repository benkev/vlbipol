help_text = '''
plot_time_var.py - Plot temporal variations of MBD, SBD, or SNR.
'''
import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_time_var.py <par> [save], ")
    print("       where <par> is either MBD or SBD or SNR.")
    print("       save (optional): save pdf of figures.")
    sys.exit(1)
    
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
with open('idx3819l.pkl', 'rb') as finp:
    idx3819l_1 = pickle.load(finp)

# with open('idx3819c.pkl', 'rb') as finp:
#     idx3819c_1 = pickle.load(finp)

with open('idx3819cI.pkl', 'rb') as finp:
    idx3819c_1 = pickle.load(finp)

if par == 'MBD' or par == 'SBD':
    ps = "(ps)"
else: # if par == 'SNR':
    ps = ""
    
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
rmse = np.zeros(nbls, dtype=float)    # Root mean square error (RMSE) for MBD
r_corr = np.zeros(nbls, dtype=float)  # Pearson's correlation for MBD

#
# Table header
#
print("BL  avg %s   rmse  relerr,%   avg bias   r_corr" % (par, ps))

if par == 'MBD':
    parname = 'mbdelay'
elif par == 'SBD':
    parname = 'sbdelay'
else:
    parname = 'snr'
    
    
#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 0   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]

    if par == 'MBD' or par == 'SBD':
        par_l_us = np.array(idx3819l_1[bl]['I'][parname])[istart:] # In useconds
        par_c_us = np.array(idx3819c_1[bl]['I'][parname])[istart:] # In useconds
        par_l = par_l_us*1e6           # Convert us to ps
        par_c = par_c_us*1e6           # Convert us to ps
    else: # if par == 'SNR':
        par_l = np.array(idx3819l_1[bl]['I'][parname])[istart:]
        par_c = np.array(idx3819c_1[bl]['I'][parname])[istart:]

    bpar = par_l - par_c                 # Bias
    par_a = (abs(par_l.mean()) + abs(par_c.mean()))/2 # Avg Lin and Cir means
    
    par0_l = par_l - par_l.mean()        # Subtract MBD means, lin pol
    par0_c = par_c - par_c.mean()        # Subtract MBD means, cir pol
    dpar = par0_l - par0_c               # Residuals

    
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    npt = len(tim)   # Number of points for current baseline
    rmse[ibl] = np.sqrt(np.sum(dpar**2)/npt)
    r_corr[ibl] = sum(par0_l*par0_c)/np.sqrt(sum(par0_l**2)*sum(par0_c**2))

    pl.figure(fig1)
    pl.subplot(5, 3, ibl+1)
    # pl.plot(tim, par_l, label='Lin_I, mean: %.1f' % par_l.mean())
    # pl.plot(tim, par_c, label='Cir_I, mean: %.1f' % par_c.mean())
    pl.plot(tim, par_l, 'b', label='Lin_I, mean: %.1f' % par_l.mean())
    pl.plot(tim, par_l, 'k.', markersize=3)
    pl.plot(tim, par_c, 'g', label='Cir_I, mean: %.1f' % par_c.mean())
    pl.plot(tim, par_c, 'k.', markersize=3)
    pl.grid(True)
    pl.legend(loc='upper left', prop={'size': 9})
    ax1 = pl.gca()

    if par == 'SNR':
        pl.ylim(0, 6000)
    
    pl.text(.88, .02, bl, transform=ax1.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax1.transAxes, \
            fontsize=9)

    pl.figure(fig2)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, dpar, color='red')
    pl.plot(tim, dpar, 'k.', markersize=3)
    pl.grid(True)
    ax2 = pl.gca()

    if  par == 'MBD':
        pl.ylim(-25, 25)
    elif par == 'SBD':
        pl.ylim(-300, 300)
    else: # if par == 'SNR':
        pl.ylim(-100, 100)

    pl.text(.88, .90, bl, transform=ax2.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .92, "r_corr: %.6f" % r_corr[ibl], transform=ax2.transAxes, \
            fontsize=9)
    pl.text(.03, .80, "RMSE: %.4f" % rmse[ibl], transform=ax2.transAxes, \
            fontsize=9)

    pl.figure(fig3)
    pl.subplot(5, 3, ibl+1)
    pl.plot(tim, bpar, color='brown')
    pl.plot(tim, bpar, 'k.', markersize=3)
    pl.grid(True)
    ax3 = pl.gca()

    if  par == 'MBD':
        pl.ylim(-200, 250)
    elif par == 'SBD':
        pl.ylim(-1700, 1400)
    else: # if par == 'SNR':
        pl.ylim(-150, 350)

    pl.text(.88, .90, bl, transform=ax3.transAxes, fontsize=10, weight="bold")
    pl.text(.03, .02, "r_corr: %.6f" % r_corr[ibl], transform=ax3.transAxes, \
            fontsize=9)
    pl.text(.03, .90, "bias mean: %.1f" % bpar.mean(), transform=ax3.transAxes,\
            fontsize=9)

    #
    # Table
    #
    rel_err = 100*abs(rmse[ibl]/par_a)
    print("%s  %7.1f   %4.2f    %4.2f     %6.1f     %8.6f" % \
          (bl, par_a, rmse[ibl], rel_err, bpar.mean(), r_corr[ibl])) 
    
    ibl = ibl + 1
    
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))
fig3.tight_layout(rect=(0,0,1, 0.95))

pl.figure(fig1)
pl.figtext(0.20, 0.96, "Pseudo-Stokes I %s %s vs Time (min), " \
           "Lin & Cir Pol after PolConvert" % (par, ps), fontsize=11)

pl.figure(fig2)
pl.figtext(0.08, 0.96, "%s Residuals %s vs Time (min), " \
           " between Lin & Cir Pol after PolConvert (means subtracted)" \
           % (par, ps), fontsize=11)

pl.figure(fig3)
pl.figtext(0.20, 0.96, "%s Bias %s vs Time (min), " \
           " between Lin & Cir Pol after PolConvert" % (par, ps), \
           fontsize=11)

if sf:
    pl.figure(fig1)
    pl.savefig("%s_Lin_I_and_Cir_I.pdf" % par, format='pdf')
    pl.figure(fig2)
    pl.savefig("%s_Lin_I_minus_Cir_I.pdf % par", format='pdf')
    pl.figure(fig3)
    pl.savefig("%s_bias_between_Lin_I_and_Cir_I.pdf" % par, format='pdf')

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




    
