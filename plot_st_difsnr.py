help_text = '''

plot_st_difsnr.py: Plot histograms of differences between the snr values
                    before and after PolConvert, for individual stations
                    and for all the stations.
'''

import sys
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



bls = list(idx3819l_1.keys())   # List of baselines
bls.sort()                      # Lexigraphically sorted baselines
nbls = len(bls)

#
# To start plotting from istart;  exclude bad data before istart.
#
# istart = 2
istart = 0


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
# Gather SNR sata into nsts station bins 
#
stsnr = {} # Dict for stationwise SNR data: stsnr['X'] 
stbls = {} # Dict for baselines including a station and their point numbers

for sta in ststr:
    
    tim = np.empty(0, dtype=float)   # Time for a particular station
    snr_l = np.empty(0, dtype=float) # Lin SNR for a particular station
    snr_c = np.empty(0, dtype=float) # Cir SNR for a particular station
    bsl = []  # List of baselines that include a particular station "sta"
    bsnpl = []  # List of numbers of points in the baselines with station "sta"

    ndat_st = 0 # Number of points for baselines with a station "sta"
    for bl in bls:   # Loop over the baselines
        if sta in bl:
            tim0 = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
            tim0 = tim0 - tim0[0]

            snr0_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
            snr0_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
            
            #
            # Subtract SNR means
            #
            snr0_l = snr0_l - snr0_l.mean()
            snr0_c = snr0_c - snr0_c.mean()

            tim = np.append(tim, tim0)
            snr_l = np.append(snr_l, snr0_l)
            snr_c= np.append(snr_c, snr0_c)
            
            ntim = len(tim0)
            ndat_st = ndat_st + ntim
            bsl.append(bl)
            bsnpl.append(ntim)
    print("'", sta, "': ", ndat_st) 
    #
    # Differences Lin-Cir for baselines with a particular station sta
    #
    dsnr = snr_l - snr_c
    stsnr[sta] = dsnr
    stbls[sta] = [bsl, bsnpl]


nbin_ini = 21   # Histogram bins
    
fig1 = pl.figure(figsize=(8, 10))
    
#
# Plot histograms of SNR differences for the baselines including station "sta"
#
hw = 40  # Histogram width: +- hw

ist = 0   # Baseline number starting from 0
for sta in ststr:
    iplt = ist + 1  # Subplot number
    pl.figure(fig1)
    pl.subplot(3, 2, iplt)
    pl.hist(stsnr[sta], nbin_ini, color='green')
    pl.grid(1)    
    ax = pl.gca()
    pl.text(.03, .92, "Station: "+sta, transform=ax.transAxes, fontsize=12)
    pl.text(.03, .84, "Bls: ", transform=ax.transAxes, fontsize=9)
    pl.text(.12, .84, ', '.join(stbls[sta][0]), transform=ax.transAxes, \
            fontsize=9)
    ist = ist + 1

    #
    # Testing the H0 hypothesis of stsnr[sta] normal distribution: FAILS!
    #
    ni, bedges = np.histogram(stsnr[sta], nbin_ini) # 21 bin

    N = np.sum(ni)
    binwd = bedges[1] - bedges[0]             # Bin width
    xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
    hmean = np.sum(xi*ni)/N               # Sample mean
    sig2 = np.sum(xi**2*ni/N - hmean**2)  # Sample variance sigma^2
    sig = np.sqrt(sig2)                   # Standard deviation sigma
    #
    # Fit a normal distribution to the histogram and to the stsnr[sta] data
    #
    zi = (xi - hmean)/sig                 # Standardized xi
    fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
    fni = binwd*N*fnorm              # Theoretical frequencies
    mu, stdev = norm.fit(stsnr[sta]) # Fit a normal dist. to the stsnr[sta] data
    #
    # 
    #
    idm = np.where(abs(stsnr[sta]) < stdev)[0]    # stsnt[sta] within +-stdev
    pmstd = len(idm)/N*100  # Percent of stsnr[sta] within +-stdev
    pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
    print("%s dsnr: mu = %f,    std = %f" % (sta, mu, stdev))
    print("%s dsnr: %5.2f%% within +-std;  %5.2f%% for normal" % \
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
    chi2obs = np.sum((ni - fni)**2/fni) # !!!!!!!! HUGE !!!!!!!!

    #
    # Critical value for chi^2 at p=0.95 confidence level
    #
    deg_fr = nbin - 2 - 1   # 2 params of normal distr. estimated, mu and sigma
    chi2cr = chi2.isf(0.05, df=deg_fr)
    q_chi2 = chi2obs/chi2cr  # Quotient
   
    print("Station %s:" % sta)
    print('Original binning with sparse tails (%d bins):' % nbin_ini)
    print('ni:  ', ni_ini)
    print('fni: ', fni_ini)
    print('Sparse tails grouped: (%d bins)' % nbin)
    print('ni:  ', ni)
    print('fni: ', fni)
    print("chi2obs/chi2cr = %f" % q_chi2)
    print()
    
    #
    # Smooth normal approximations 
    #
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
    if q_chi2 < 10:
        pl.text(.67, .50, "$\chi^2 > \chi^2_{cr}$=%.1f" % chi2cr, \
                transform=ax.transAxes, fontsize=9)
    else:
        pl.text(.67, .50, "$\chi^2 \gg \chi^2_{cr}$=%.1f" % chi2cr, \
                transform=ax.transAxes, fontsize=9)
    #
    # X ticks
    #
    pxtc = -100 + 50*np.arange(5, dtype=float)
    pxtc = np.insert(pxtc, (2, 3), (-stdev, stdev))

    xtc = list(np.int64(pxtc))
    xtc[2] = r"$-\sigma$"
    xtc[4] = r"$+\sigma$"

    pl.xticks(pxtc, xtc)
    pl.xlim(-120,+120)

fig1.text(0.2, 0.97, "Differences SNR Lin_I-Cir_I Distributions for Stations", \
          fontsize=12)
fig1.tight_layout(rect=(0,0,1, 0.95))



# ================= HIST FOR ALL STATIONS ===================================

nbin_ini = 21


    
#
# Get and plot SNR for all the baselines 
#
rmse_snr = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for SNR

dsnr = []  # Differences of SNR for all baselines

ibl = 0   # Baseline number starting from 0
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    snr_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
    dsnr_bl = np.zeros_like(tim)  # Differences of SNR for current baseline

    #
    # Subtract SNR means
    #
    snr0_l = snr_l - snr_l.mean()
    snr0_c = snr_c - snr_c.mean()
    
    #
    # Root mean square error (RMSE)
    #
    dsnr_bl = snr0_l - snr0_c
    dsnr.extend(dsnr_bl) # Add SNR differences to list
    npt = len(tim)   # Number of points for current baseline
    rmse_snr[ibl] = np.sqrt(np.sum(dsnr_bl**2)/nbls)
    
    ibl = ibl + 1

dsnr = np.array(dsnr, dtype=float)

fig5 = pl.figure()

pl.figure(fig5);
pl.hist(dsnr, nbin_ini, color = "g", ec="k"); pl.grid(1)
pl.xlim(-120,120)

fig5.text(0.15, 0.95, "Differences SNR Lin_I-Cir_I Distributions " \
          "for All Baselines", \
          fontsize=12)
fig5.tight_layout(rect=(0,0,1, 0.95))

#
# Testing the H0 hypothesis or dsnr normal distribution: FAILS!
#
ni, bedges = np.histogram(dsnr, nbin_ini) # 21 bin

# ni = ni[7:15]
# bedges = bedges[7:16]

N = np.sum(ni)
binwd = bedges[1] - bedges[0]             # Bin width
xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
hmean = np.sum(xi*ni)/N               # Sample mean
sig2 = np.sum(xi**2*ni/N - hmean**2)  # Sample variance sigma^2
sig = np.sqrt(sig2)                   # Standard deviation sigma
#
# Fit a normal distribution to the histogram and to the whole dsnr data
#
zi = (xi - hmean)/sig                 # Standardized xi
fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
fni = binwd*N*fnorm              # Theoretical frequencies
mu, stdev = norm.fit(dsnr)  # Fit a normal distribution to the WHOLE dsnr data
#
# 
#
idm = np.where(abs(dsnr) < stdev)[0]    # Count of dsnr within +-stdev 
# idm1 = np.where(abs(dsnr) >= stdev)[0]  # Count of dsnr outside of +-stdev
pmstd = len(idm)/N*100  # Percent of dsnr within +-stdev
pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
# print("Histogram:  hmean = %f, sig = %f" % (hmean, sig))
print("dsnr: mu = %f,    stdev = %f" % (mu, stdev))
print("dsnr: %5.2f%% within +-std;  %5.2f%% for normal" % (pmstd, pmstd_norm))


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
chi2obs = np.sum((ni - fni)**2/fni)

#
# Critical value for chi^2 at p=0.95 confidence level
#
deg_fr = nbin - 2 - 1    # 2 params of normal distr. estimated, mu and sigma
chi2cr = chi2.isf(0.05, df=deg_fr)
q_chi2 = chi2obs/chi2cr  # Quotient

print('All stations:')
print('Original binning with sparse tails (%d bins):' % nbin_ini)
print('ni:  ', ni_ini)
print('fni: ', fni_ini)
print('Sparse tails grouped: (%d bins)' % nbin)
print('ni:  ', ni)
print('fni: ', fni)
print()

#
# Smooth normal approximations 
#
x1 = np.linspace(-100, 100, 101)
#f1 = norm.pdf(x1, hmean, sig)*binwd*N
f2 = norm.pdf(x1, mu, stdev)*binwd*N


pl.figure(fig5)

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
        transform=ax.transAxes, fontsize=9)
pl.text(.75, .85, "df = %2d-2-1=%2d" % (nbin, deg_fr), \
        transform=ax.transAxes, fontsize=9)
if chi2obs < 1000:
    pl.text(.75, .80, "$\chi^2$=%6.2f" % chi2obs, transform=ax.transAxes, \
            fontsize=10)
else:
    pl.text(.75, .80, "$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
            fontsize=10)
if q_chi2 < 10:
    pl.text(.75, .75, "$\chi^2 > \chi^2_{cr}$ = %.1f" % chi2cr, \
            transform=ax.transAxes, fontsize=9)
else:
    pl.text(.75, .75, "$\chi^2 \gg \chi^2_{cr}$ = %.1f" % chi2cr, \
            transform=ax.transAxes, fontsize=9)


#
# X ticks
#
pxtc = -100 + 50*np.arange(5, dtype=float)
pxtc = np.insert(pxtc, (2, 3), (-stdev, stdev))

xtc = list(np.int64(pxtc))
xtc[2] = r"$-\sigma$"
xtc[4] = r"$+\sigma$"


pl.xticks(pxtc, xtc)

pl.show()


pl.figure(fig1)
pl.savefig("Distr_SNR_Lin_I-Cir_I_Diff_Stations.eps", format='eps')
pl.figure(fig5)
pl.savefig("Distr_SNR_Lin_I-Cir_I_Diff.eps", format='eps')





    
