help_text = '''

plot_st_difhist.py: Plot histograms of residuals between the pseudo-Stokes
                    parameter values before and after PolConvert,
                    for the baselines involving individual stations
                    and for the baselines of all the stations.
'''

import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_st_difhist.py <par> [save], ")
    print("       where <par> is either MBD or SBD or SNR.")
    print("       save (optional): save figures in pdf format.")
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
    
import pickle
import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm, chi2

from group_tails import group_tails

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

if par == 'MBD' or par == 'SBD':
    ps = "(ps)"
else: # if par == 'SNR':
    ps = ""
    


bls = list(idx3819l_1.keys())   # List of baselines
bls.sort()                      # Lexigraphically sorted baselines
nbls = len(bls)

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2


# Set of station letters stset
ststr = ''
for bl in bls: ststr = ststr + bl  # Concatenate baseline strings in ststr
stset = set(ststr)  # Leave only unique station letters in the sts set

# String of station letters ststr
nsts = len(stset)
# ststr = ''
# for st in stset: ststr = ststr + st
ststr = ''.join(sorted(stset))

#
# Gather par data into nsts station bins 
#
stpar = {} # Dict for stationwise par data: stpar['X'] 
stbls = {} # Dict for baselines including a station and their point numbers
stsnr = {} # Dict for average SNR for a station

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
strmse = {}   # Root mean square error (RMSE) for par for a station
str_corr = {} # Pearson's correlation for par for a station

if par == 'MBD':
    parname = 'mbdelay'
elif par == 'SBD':
    parname = 'sbdelay'
else:
    parname = 'snr'
    

for sta in ststr:
    
    tim = np.empty(0, dtype=float)   # Time for a station
    par_l = np.empty(0, dtype=float) # Lin par for a station
    par_c = np.empty(0, dtype=float) # Cir par for a station
    par0_l = np.empty(0, dtype=float) # Lin par-par.mean() for a station
    par0_c = np.empty(0, dtype=float) # Cir par-par.mean() for a station
    
    snr_a = np.empty(0, dtype=float)  # SNR avg for a station
       
    bsl = []  # List of baselines that include a station "sta"
    bsnpl = []  # List of numbers of points in the baselines with station "sta"

    ndat_st = 0 # Number of points for baselines with a station "sta"
    for bl in bls:   # Loop over the baselines
        if sta in bl:
            timbl = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
            timbl0 = timbl - timbl[0]

            snrbl_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
            snrbl_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
            snrbl_a = (abs(snrbl_l.mean()) + abs(snrbl_c.mean()))/2


            if par == 'MBD' or par == 'SBD':
                parbl_l_us = np.array(idx3819l_1[bl]['I'][parname])[istart:]
                parbl_c_us = np.array(idx3819c_1[bl]['I'][parname])[istart:]
                parbl_l = parbl_l_us*1e6           # Convert us to ps
                parbl_c = parbl_c_us*1e6           # Convert us to ps
            else: # if par == 'SNR':
                parbl_l = np.copy(snrbl_l)
                parbl_c = np.copy(snrbl_c)

            # bparbl = parbl_l - parbl_c                 # Bias
            parbl_a = (abs(parbl_l.mean()) + abs(parbl_c.mean()))/2

            parbl0_l = parbl_l - parbl_l.mean()  # Subtract par means, lin pol
            parbl0_c = parbl_c - parbl_c.mean()  # Subtract par means, cir pol

            tim = np.append(tim, timbl0)
            par_l = np.append(par_l, parbl_l)
            par_c = np.append(par_c, parbl_c)
            par0_l = np.append(par0_l, parbl0_l) # Mean subtracted
            par0_c = np.append(par0_c, parbl0_c) # Mean subtracted

            snr_a = np.append(snr_a, snrbl_a)
            
            ntim = len(timbl)
            ndat_st = ndat_st + ntim
            
            bsl.append(bl)
            bsnpl.append(ntim)
            
    print("'", sta, "': ", ndat_st) 
    #
    # Residuals Lin-Cir for baselines with a particular station sta
    #
    dpar = par0_l - par0_c
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    strmse[sta] = np.sqrt(np.sum(dpar**2)/ndat_st)
    str_corr[sta] = sum(par0_l*par0_c)/np.sqrt(sum(par0_l**2)*sum(par0_c**2))
    stsnr[sta] = snr_a.mean() 
    stpar[sta] = dpar
    stbls[sta] = [bsl, bsnpl]
    
#
# Set default array print format: as fixed-point only and as short as possible
#
np.set_printoptions(suppress=True, precision=1)

nbin_ini = 21   # Initial number of histogram bins (before tail grouping)
    
fig1 = pl.figure(figsize=(8, 10))
    
#
# Plot histograms of par residuals for the baselines including station "sta"
#
ist = 0   # Baseline number starting from 0
for sta in ststr:
    iplt = ist + 1  # Subplot number
    pl.figure(fig1)
    pl.subplot(3, 2, iplt)
    pl.hist(stpar[sta], nbin_ini, color='green')
    
    if  par == 'MBD':
        pl.xlabel("ps")
    elif par == 'SBD':
        pl.xlabel("ps")
    else: # if par == 'SNR':
        pass

    pl.grid(1)    
    ax = pl.gca()
    pl.text(.03, .92, "Station: "+sta, transform=ax.transAxes, fontsize=12)
    pl.text(.03, .84, "Bls: ", transform=ax.transAxes, fontsize=9)
    pl.text(.12, .84, ', '.join(stbls[sta][0]), transform=ax.transAxes, \
            fontsize=9)
    # ist = ist + 1

    #
    # Testing the H0 hypothesis of stpar[sta] normal distribution: FAILS!
    #
    ni, bedges = np.histogram(stpar[sta], nbin_ini) # 21 bin

    N = np.sum(ni)
    binwd = bedges[1] - bedges[0]             # Bin width
    xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
    hmean = np.sum(xi*ni)/N               # Sample mean
    sig2 = np.sum(xi**2*ni/N - hmean**2)  # Sample variance sigma^2
    sig = np.sqrt(sig2)                   # Standard deviation sigma
    #
    # Fit a normal distribution to the histogram and to the stpar[sta] data
    #
    zi = (xi - hmean)/sig                 # Standardized xi
    fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
    fni = binwd*N*fnorm              # Theoretical frequencies
    mu, stdev = norm.fit(stpar[sta]) # Fit a normal dist. to the stpar[sta] data
    #
    # 
    #
    idm = np.where(abs(stpar[sta]) < stdev)[0]    # stpar[sta] within +-stdev
    pmstd = len(idm)/N*100  # Percent of stpar[sta] within +-stdev
    pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
    print("%s dmdb: mu = %f,    std = %f" % (sta, mu, stdev))
    print("%s dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % \
          (sta, pmstd, pmstd_norm))

    #
    # Group left-tail and right-tail bins with sparse data.
    #
    ni_ini = np.copy(ni)    # Save old observed freqs with sparse tails
    fni_ini = np.copy(fni)  # Save old theoretical freqs with sparse tails

    ni, fni = group_tails(ni_ini, fni_ini) 
    nbin = len(ni)
    
    #
    # Pearson's X^2
    #
    chi2obs = np.sum((ni - fni)**2/fni)

    #
    # Critical value for chi^2 at p=0.95 confidence level
    #
    deg_fr = nbin - 2 - 1   # 2 params of normal distr. estimated, mu and sigma
    chi2cr = chi2.isf(0.05, df=deg_fr) # The same as chi2.ppf(0.95, df=deg_fr)
    q_chi2 = chi2obs/chi2cr  # Quotient

    # print("Station %s:" % sta)
    # print('Original binning with sparse tails (%d bins):' % nbin_ini)
    # print('ni:  ', ni_ini)
    # print('fni: ', fni_ini)
    # print('Sparse tails grouped: (%d bins)' % nbin)
    # print('ni:  ', ni)
    # print('fni: ', fni)
    # print()
    print("%s nbin = %d, chi2obs = %.1f, chi2cr = %.1f chi2obs/chi2cr = %.1f" %
          (sta, nbin, chi2obs, chi2cr, q_chi2))
    
    #
    # Smooth normal approximations 
    #
    if  par == 'MBD':
        x1 = np.linspace(-11, 11, 101)
    elif par == 'SBD':
        x1 = np.linspace(-200, 200, 101)
    else: # if par == 'SNR':
        x1 = np.linspace(-100, 100, 101)

    f2 = norm.pdf(x1, mu, stdev)*binwd*N

    pl.plot(x1, f2, 'b')

    stdh = 0.7*pl.ylim()[1]  # Height of std line
    pl.plot([-stdev, -stdev], [0, stdh], 'r-.')   # , lw=0.8)
    pl.plot([stdev, stdev], [0, stdh], 'r-.')     # , lw=0.8)

    pl.plot(xi, fni_ini, 'r.')

   
    ax = pl.gca()
    
    pl.text(.6, .92, "Within $\pm$std: %5.2f%%" % pmstd, \
            transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .85, "For Normal: 68.27%", transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .78, "$\mu$=%.4f, $\sigma$=%5.2f" % (mu, stdev), \
            transform=ax.transAxes, fontsize=9)
    pl.text(.67, .71, "N=%3d; bins: %2d" % (N, nbin), \
            transform=ax.transAxes, fontsize=9)
    pl.text(.67, .64, "df=%2d-2-1=%2d" % (nbin, deg_fr), \
            transform=ax.transAxes, fontsize=9)
    if chi2obs < 1000:
        pl.text(.67, .57, "$\chi^2$=%6.2f" % chi2obs, transform=ax.transAxes, \
                fontsize=9)
    else:
        pl.text(.67, .57, "$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
                fontsize=9)
    if q_chi2 <= 1:
        pl.text(.66, .50, "$\chi^2 \leq \chi^2_{cr}$=%.2f" % chi2cr, \
                transform=ax.transAxes, fontsize=11, c='red')        
    elif q_chi2 < 5:
        pl.text(.67, .50, "$\chi^2 > \chi^2_{cr}$=%.2f" % chi2cr, \
                transform=ax.transAxes, fontsize=9)
    else:
        pl.text(.67, .50, "$\chi^2 \gg \chi^2_{cr}$=%.2f" % chi2cr, \
                transform=ax.transAxes, fontsize=9)

    # if  par == 'MBD' or par == 'SBD':
    #     pl.text(.67, .41, "$\overline{SNR}$: %.1f" % stsnr[sta], \
    #             transform=ax.transAxes, fontsize=9)
    #     pl.text(.67, .36, "r_corr: %.6f" % str_corr[sta], \
    #             transform=ax.transAxes, fontsize=9)
    # else: # if  par == 'SNR':
    #     pl.text(.67, .42, "r_corr: %.6f" % str_corr[sta], \
    #             transform=ax.transAxes, fontsize=9)
    pl.text(.67, .41, "$\overline{\mathrm{SNR}}$: %.1f" % stsnr[sta], \
            transform=ax.transAxes, fontsize=9)
    pl.text(.67, .36, "r_corr: %.6f" % str_corr[sta], \
            transform=ax.transAxes, fontsize=9)

    #
    # X ticks
    #
    if  par == 'MBD':
        pxtc = -20 + 5*np.arange(9, dtype=float) # X ticks positions
    
        print("pxtc     = ", pxtc)

        i_mstd = 4    # Position of -stdev tick
        i_pstd = 6    # Position of +stdev tick
        pl.xlabel("ps")
        #hw = 12            # Histogram width
        hw = 23            # Histogram width
        
    elif par == 'SBD':
        pxtc = -300 + 100*np.arange(7, dtype=float) # X ticks positions
        i_mstd = 3    # Position of -stdev
        i_pstd = 5    # Position of +stdev
        pl.xlabel("ps")
        hw = 300           # Histogram width

    else: # if par == 'SNR':
        pxtc = -100 + 50*np.arange(7, dtype=float) # X ticks positions
        i_mstd = 2    # Position of -stdev
        i_pstd = 4    # Position of +stdev
        hw = 120           # Histogram width

    pxtc = np.insert(pxtc, i_mstd, -stdev)
    pxtc = np.insert(pxtc, i_pstd, +stdev)
    xtc = list(np.int64(pxtc))
    xtc[i_mstd] = r"$-\sigma$"    # Tick label for -stdev
    xtc[i_pstd] = r"$+\sigma$"    # Tick label for +stdev

    pl.xticks(pxtc, xtc)
    pl.xlim(-hw,+hw)

    ist = ist + 1


fig1.text(0.2, 0.96, \
          "Distributions of %s Residuals Lin_I-Cir_I for Stations" % par, \
          fontsize=12)
fig1.tight_layout(rect=(0,0,1, 0.95))



# ================= HIST FOR ALL STATIONS ===================================

nbin_ini = 21

#
# Get and plot par for all the baselines 
#
# rmse: Root mean square errors (RMSE) between lin pol and cir pol curves
# r_corr: Correlation coefficients  between lin pol and cir pol curves
#
# WRONG:
# rmse_r: RMSE reduced with respect to the absolute value of the mean
#         of the average between lin pol and cir pol curves

tim = np.empty(0, dtype=float)   # Time
par_l = np.empty(0, dtype=float) # Lin PAR
par_c = np.empty(0, dtype=float) # Cir par
snr_a = np.empty(0, dtype=float) # Lin and Cir SNR average
par0_l = np.empty(0, dtype=float) # Lin par, mean subtracted
par0_c = np.empty(0, dtype=float) # Cir par, mean subtracted

rmse = np.empty(0, dtype=float)   # Root mean square error (RMSE) for par
# rmse_r = np.empty(0, dtype=float) # RMSE reduced wrt abs of average
r_corr = np.empty(0, dtype=float) # Pearson's correlation for par

for bl in bls:   # Loop over the baselines
    timbl = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    timbl0 = timbl - timbl[0]

    snrbl_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snrbl_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
    snrbl_a = (abs(snrbl_l.mean()) + abs(snrbl_c.mean()))/2


    if par == 'MBD' or par == 'SBD':
        parbl_l_us = np.array(idx3819l_1[bl]['I'][parname])[istart:]
        parbl_c_us = np.array(idx3819c_1[bl]['I'][parname])[istart:]
        parbl_l = parbl_l_us*1e6           # Convert us to ps
        parbl_c = parbl_c_us*1e6           # Convert us to ps
    else: # if par == 'SNR':
        parbl_l = np.copy(snrbl_l)
        parbl_c = np.copy(snrbl_c)

    # bparbl = parbl_l - parbl_c                 # Bias
    parbl_a = (abs(parbl_l.mean()) + abs(parbl_c.mean()))/2

    parbl0_l = parbl_l - parbl_l.mean()  # Subtract par means, lin pol
    parbl0_c = parbl_c - parbl_c.mean()  # Subtract par means, cir pol

    tim = np.append(tim, timbl0)
    par_l = np.append(par_l, parbl_l)
    par_c = np.append(par_c, parbl_c)
    par0_l = np.append(par0_l, parbl0_l) # Mean subtracted
    par0_c = np.append(par0_c, parbl0_c) # Mean subtracted
    snr_a = np.append(snr_a, snrbl_a) # Averages of all SNRs, Lin and Cir

#
# Residuals Lin-Cir for all baselines
#
dpar = par0_l - par0_c
#
# Root mean square error (RMSE) and Pearson's correlation coefficient
#
ndat = len(tim)
rmse = np.sqrt(np.sum(dpar**2)/ndat)
r_corr = sum(par0_l*par0_c)/np.sqrt(sum(par0_l**2)*sum(par0_c**2))
snr_avg = snr_a.mean()

print("All baselines: abs(par_l).mean() = %.2f (ps),\t "
          "abs(par_c).mean() = %.2f (ps)" % \
          (abs(par_l).mean(), abs(par_c).mean()))
print("All baselines: dpar min and max: ", dpar.min(), dpar.max())

fig2 = pl.figure()

pl.figure(fig2);
pl.hist(dpar, nbin_ini, color = "g", ec="k"); pl.grid(1)

if  par == 'MBD':
    pl.xlabel("ps")
elif par == 'SBD':
    pl.xlabel("ps")
else: # if par == 'SNR':
    pass

        
fig2.text(0.15, 0.95, "Distribution of %s Residuals Lin_I-Cir_I " \
          "for All Baselines" % par, \
          fontsize=12)
fig2.tight_layout(rect=(0,0,1, 0.95))


#
# Testing the H0 hypothesis or dpar normal distribution: FAILS!
#
ni, bedges = np.histogram(dpar, nbin_ini) # 21 bin

# ni = ni[7:15]
# bedges = bedges[7:16]

N = np.sum(ni)
binwd = bedges[1] - bedges[0]             # Bin width
xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
hmean = np.sum(xi*ni)/N               # Sample mean
sig2 = np.sum(xi**2*ni/N - hmean**2)  # Sample variance sigma^2
sig = np.sqrt(sig2)                   # Standard deviation sigma
#
# Fit a normal distribution to the histogram and to the whole dpar data
#
zi = (xi - hmean)/sig                 # Standardized xi
fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
fni = binwd*N*fnorm              # Theoretical frequencies
mu, stdev = norm.fit(dpar)  # Fit a normal distribution to the WHOLE dpar data
#
# 
#
idm = np.where(abs(dpar) < stdev)[0]    # Count of dpar within +-stdev 
# idm1 = np.where(abs(dpar) >= stdev)[0]  # Count of dpar outside of +-stdev
pmstd = len(idm)/N*100  # Percent of dpar within +-stdev
pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
# print("Histogram:  hmean = %f, sig = %f" % (hmean, sig))
print("dmdb: mu = %f,    stdev = %f" % (mu, stdev))
print("dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % (pmstd, pmstd_norm))


#
# Group left-tail and right-tail bins with sparse data.
#

ni_ini = np.copy(ni)     # Save old freqs with sparse tails
fni_ini = np.copy(fni)    # Save old freqs with sparse tails

ni, fni = group_tails(ni_ini, fni_ini) 

nbin = len(ni)

#
# Pearson's X^2
#
chi2obs = np.sum((ni - fni)**2/fni) # !!!!!!!! HUGE !!!!!!!!

#
# Critical value for chi^2 at p=0.95 confidence level
#
deg_fr = nbin - 2 - 1    # 2 params of normal distr. estimated, mu and sigma
chi2cr = chi2.isf(0.05, df=deg_fr)
q_chi2 = chi2obs/chi2cr  # Quotient

# print('All stations:')
# print('Original binning with sparse tails (%d bins):' % nbin_ini)
# print('ni:  ', ni_ini)
# print('fni: ', fni_ini)
# print('Sparse tails grouped: (%d bins)' % nbin)
# print('ni:  ', ni)
# print('fni: ', fni)
# print("chi2obs/chi2cr = %f" % q_chi2)
# print()
print("%s nbin = %d, chi2obs = %.1f, chi2cr = %.1f chi2obs/chi2cr = %.1f" %
      (sta, nbin, chi2obs, chi2cr, q_chi2))

#
# Smooth normal approximations 
#
if  par == 'MBD':
    x1 = np.linspace(-11, 11, 101)
elif par == 'SBD':
    x1 = np.linspace(-200, 200, 101)
else: # if par == 'SNR':
    x1 = np.linspace(-100, 100, 101)
    
f2 = norm.pdf(x1, mu, stdev)*binwd*N


pl.figure(fig2)

pl.plot(x1, f2, 'b')

stdh = 0.85*pl.ylim()[1]  # Height of std line
pl.plot([-stdev, -stdev], [0, stdh], 'r-.')
pl.plot([stdev, stdev], [0, stdh], 'r-.')

pl.plot(xi, fni_ini, 'ro')

ax = pl.gca()
pl.text(.04, .95, "Within $\pm$std: %5.2f%%" % pmstd, transform=ax.transAxes, \
        fontsize=10)
pl.text(.04, .90, "For Normal: 68.27%", transform=ax.transAxes, \
        fontsize=10)
pl.text(.75, .95, "$\mu$=%.4f, $\sigma$=%5.2f" % (mu, stdev), \
        transform=ax.transAxes, fontsize=10)
pl.text(.75, .90, "N = %3d; bins: %2d" % (N, nbin), \
        transform=ax.transAxes, fontsize=10)
pl.text(.75, .85, "df = %2d-2-1=%2d" % (nbin, deg_fr), \
        transform=ax.transAxes, fontsize=10)
if chi2obs < 1000:
    pl.text(.75, .80, "$\chi^2$=%6.2f" % chi2obs, transform=ax.transAxes, \
            fontsize=10)
else:
    pl.text(.75, .80, "$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
            fontsize=10)

if q_chi2 <= 1:
    pl.text(.75, .74, "$\chi^2 \leq \chi^2_{cr}$=%.2f" % chi2cr, \
            transform=ax.transAxes, fontsize=11, c='red')        
elif q_chi2 < 5:
    pl.text(.75, .75, "$\chi^2 > \chi^2_{cr}$ = %.2f" % chi2cr, \
            transform=ax.transAxes, fontsize=10)
else:
    pl.text(.75, .75, "$\chi^2 \gg \chi^2_{cr}$ = %.2f" % chi2cr, \
            transform=ax.transAxes, fontsize=10)
    
#        ---------------   Include avg(SNR) in the SNR plots ? ----------
# if  par == 'MBD' or par == 'SBD':
#     pl.text(.75, .68, "$\overline{\mathrm{SNR}}$: %.1f" % snr_avg, \
#             transform=ax.transAxes, fontsize=10)
#     pl.text(.75, .64, "r_corr: %.6f" % r_corr, transform=ax.transAxes, \
#             fontsize=10)
# else: # if  par == 'SNR':
#     pl.text(.75, .70, "r_corr: %.6f" % r_corr, transform=ax.transAxes, \
#             fontsize=10)

pl.text(.75, .68, "$\overline{\mathrm{SNR}}$: %.1f" % snr_avg, \
        transform=ax.transAxes, fontsize=10)
pl.text(.75, .64, "r_corr: %.6f" % r_corr, transform=ax.transAxes, \
        fontsize=10)


#
# X ticks
#
if  par == 'MBD':
    pxtc = -20 + 5*np.arange(9, dtype=float) # X ticks positions
    i_mstd = 4    # Position of -stdev tick
    i_pstd = 6    # Position of +stdev tick
    pl.xlabel("ps")
    #hw = 12            # Histogram width
    hw = 23            # Histogram width

elif par == 'SBD':
    pxtc = -300 + 100*np.arange(7, dtype=float) # X ticks positions
    i_mstd = 3    # Position of -stdev
    i_pstd = 5    # Position of +stdev
    pl.xlabel("ps")
    hw = 300           # Histogram width
    
else: # if par == 'SNR':
    pxtc = -100 + 50*np.arange(7, dtype=float) # X ticks positions
    i_mstd = 2    # Position of -stdev
    i_pstd = 3    # Position of +stdev
    hw = 120           # Histogram width

    
pxtc = np.insert(pxtc, i_mstd, -stdev)
pxtc = np.insert(pxtc, i_pstd, +stdev)
xtc = list(np.int64(pxtc))
xtc[i_mstd] = r"$-\sigma$"    # Tick label for -stdev
xtc[i_pstd] = r"$+\sigma$"    # Tick label for +stdev


pl.xticks(pxtc, xtc)
pl.xlim(-hw,+hw)

pl.show()

if sf:
    pl.figure(fig1)
    pl.savefig("Distr_%s_Lin_I-Cir_I_Diff_Stations.pdf" % par, format='pdf')
    pl.figure(fig2)
    pl.savefig("Distr_%s_Lin_I-Cir_I_Diff.pdf" % par, format='pdf')


#
# Restore default array print format
#
np.set_printoptions(suppress=False, precision=8)







    
