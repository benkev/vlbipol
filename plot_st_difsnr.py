help_text = '''

plot_st_difsnr.py: Plot histograms of residuals between the snr values
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
# Save figures indicator
#
sf = False  # Do not save

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
# Gather SNR data into nsts station bins 
#
stsnr = {} # Dict for stationwise SNR data: stsnr['X'] 
stbls = {} # Dict for baselines including a station and their point numbers

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
strmse = {}   # Root mean square error (RMSE) for SNR for a station
# strmse_r = {} # RMSE reduced wrt abs of average for a station
str_corr = {} # Pearson's correlation for SNR for a station

for sta in ststr:
    
    tim = np.empty(0, dtype=float)   # Time for a particular station
    snr_l = np.empty(0, dtype=float) # Lin SNR for a particular station
    snr_c = np.empty(0, dtype=float) # Cir SNR for a particular station
    snr0_l = np.empty(0, dtype=float) # Lin SNR for a particular station
    snr0_c = np.empty(0, dtype=float) # Cir SNR for a particular station
       
    bsl = []  # List of baselines that include a particular station "sta"
    bsnpl = []  # List of numbers of points in the baselines with station "sta"

    ndat_st = 0 # Number of points for baselines with a station "sta"
    for bl in bls:   # Loop over the baselines
        if sta in bl:
            timbl = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
            timbl0 = timbl - timbl[0]

            snrbl_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
            snrbl_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
            
            #
            # Subtract SNR means
            #
            snrbl0_l = snrbl_l - snrbl_l.mean()        # Subtract SNR means
            snrbl0_c = snrbl_c - snrbl_c.mean()        # Subtract SNR means

            tim = np.append(tim, timbl0)
            snr_l = np.append(snr_l, snrbl_l)
            snr_c = np.append(snr_c, snrbl_c)
            snr0_l = np.append(snr0_l, snrbl0_l) # Zero subtracted
            snr0_c = np.append(snr0_c, snrbl0_c) # Zero subtracted
            
            ntim = len(timbl)
            ndat_st = ndat_st + ntim
            
            bsl.append(bl)
            bsnpl.append(ntim)
            
    print("'", sta, "': ", ndat_st) 
    #
    # Residuals Lin-Cir for baselines with a particular station sta
    #
    dsnr = snr0_l - snr0_c
    #
    # Root mean square error (RMSE) and Pearson's correlation coefficient
    #
    strmse[sta] = np.sqrt(np.sum(dsnr**2)/ndat_st)
    # snr_a = (snr_l + snr_c)/2      # Average of the lin and cir curves
    # strmse_r[sta] = strmse[sta]/abs(snr_a.mean()) # RMSE reduced
    str_corr[sta] = sum(snr0_l*snr0_c)/np.sqrt(sum(snr0_l**2)*sum(snr0_c**2))

    stsnr[sta] = dsnr
    stbls[sta] = [bsl, bsnpl]
    
#
# Set default array print format: as fixed-point only and as short as possible
#
np.set_printoptions(suppress=True, precision=1)

nbin_ini = 21   # Initial number of histogram bins (before tail grouping)
    
fig1 = pl.figure(figsize=(8, 10))
    
#
# Plot histograms of SNR residuals for the baselines including station "sta"
#

ist = 0   # Baseline number starting from 0
for sta in ststr:
    iplt = ist + 1  # Subplot number
    pl.figure(fig1)
    pl.subplot(3, 2, iplt)
    pl.hist(stsnr[sta], nbin_ini, color='green')
    pl.xlim(-120, 120)
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

    # ni = ni[7:15]
    # bedges = bedges[7:16]

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
    idm = np.where(abs(stsnr[sta]) < stdev)[0]    # stsnr[sta] within +-stdev
    pmstd = len(idm)/N*100  # Percent of stsnr[sta] within +-stdev
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
    chi2cr = chi2.isf(0.05, df=deg_fr)
    q_chi2 = chi2obs/chi2cr  # Quotient

    print("Station %s:" % sta)
    print('Original binning with sparse tails (%d bins):' % nbin_ini)
    print('ni:  ', ni_ini)
    print('fni: ', fni_ini)
    print('Sparse tails grouped: (%d bins)' % nbin)
    print('ni:  ', ni)
    print('fni: ', fni)
    print("chi2obs/chi2cr = %f" % (chi2obs/chi2cr))
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

    # pl.xlim(-120,120)
    
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

    pl.text(.67, .43, "RMSE: %.4f" % strmse[sta], transform=ax.transAxes, \
            fontsize=9)
    # pl.text(.67, .36, "RMSE_r: %.5f" % strmse_r[sta], transform=ax.transAxes,\
    #         fontsize=9)
    pl.text(.67, .36, "r_corr: %.6f" % str_corr[sta], transform=ax.transAxes, \
            fontsize=9)


    #
    # X ticks
    #
    pxtc = -100 + 50*np.arange(7, dtype=float)
    pxtc = np.insert(pxtc, (2, 3), (-stdev, stdev))
    
    xtc = list(np.int64(pxtc))
    xtc[2] = r"$-\sigma$"
    xtc[4] = r"$+\sigma$"

    pl.xticks(pxtc, xtc)
    pl.xlim(-120, +120)


fig1.text(0.2, 0.96, "Distributions of SNR Residuals Lin_I-Cir_I for Stations",
          fontsize=12)
fig1.tight_layout(rect=(0,0,1, 0.95))



# ================= HIST FOR ALL STATIONS ===================================

nbin_ini = 21

#
# Get and plot SNR for all the baselines 
#
# rmse: Root mean square errors (RMSE) between lin pol and cir pol curves
# r_corr: Correlation coefficients  between lin pol and cir pol curves
#
# WRONG:
# rmse_r: RMSE reduced with respect to the absolute value of the mean
#         of the average between lin pol and cir pol curves
#

tim = np.empty(0, dtype=float)   # Time
snr_l = np.empty(0, dtype=float) # Lin SNR
snr_c = np.empty(0, dtype=float) # Cir SNR
snr0_l = np.empty(0, dtype=float) # Lin SNR, mean subtracted
snr0_c = np.empty(0, dtype=float) # Cir SNR, mean subtracted

rmse = np.empty(0, dtype=float)   # Root mean square error (RMSE) for SNR
# rmse_r = np.empty(0, dtype=float) # RMSE reduced wrt abs of average
r_corr = np.empty(0, dtype=float) # Pearson's correlation for SNR

for bl in bls:   # Loop over the baselines
    timbl = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    timbl0 = timbl - timbl[0]

    snrbl_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snrbl_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]

    #
    # Subtract SNR means
    #
    snrbl0_l = snrbl_l - snrbl_l.mean()        # Subtract SNR means
    snrbl0_c = snrbl_c - snrbl_c.mean()        # Subtract SNR means

    tim = np.append(tim, timbl0)
    snr_l = np.append(snr_l, snrbl_l)
    snr_c = np.append(snr_c, snrbl_c)
    snr0_l = np.append(snr0_l, snrbl0_l) # Zero subtracted
    snr0_c = np.append(snr0_c, snrbl0_c) # Zero subtracted


#
# Residuals Lin-Cir for all baselines
#
dsnr = snr0_l - snr0_c
#
# Root mean square error (RMSE) and Pearson's correlation coefficient
#
ndat = len(tim)
rmse = np.sqrt(np.sum(dsnr**2)/ndat)
# snr_a = (snr_l + snr_c)/2       # Average of the lin and cir curves
# rmse_r = rmse/abs(snr_a.mean()) # RMSE reduced wrt abs average
r_corr = sum(snr0_l*snr0_c)/np.sqrt(sum(snr0_l**2)*sum(snr0_c**2))


print("All baselines: abs(snr_l).mean() = %.2f,\t "
          "abs(snr_c).mean() = %.2f" % \
          (abs(snr_l).mean(), abs(snr_c).mean()))
print("All baselines: dsnr min and max: ", dsnr.min(), dsnr.max())

fig5 = pl.figure()

pl.figure(fig5);
pl.hist(dsnr, nbin_ini, color = "g", ec="k"); pl.grid(1)
pl.xlim(-120, 120)

fig5.text(0.15, 0.95, "Distribution of SNR Residuals Lin_I-Cir_I " \
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
print("chi2obs/chi2cr = %f" % q_chi2)
print()

#
# Smooth normal approximations 
#
x1 = np.linspace(-100, 100, 101)
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

pl.text(.75, .70, "RMSE: %.4f" % rmse, transform=ax.transAxes, \
        fontsize=10)
# pl.text(.75, .65, "RMSE_r: %.5f" % rmse_r, transform=ax.transAxes, \
#         fontsize=10)
pl.text(.75, .65, "r_corr: %.6f" % r_corr, transform=ax.transAxes, \
        fontsize=10)

#
# X ticks
#
pxtc = -100 + 50*np.arange(5, dtype=float)
pxtc = np.insert(pxtc, (2, 3), (-stdev, stdev))

xtc = list(np.int64(pxtc))
xtc[2] = r"$-\sigma$"
xtc[4] = r"$+\sigma$"

pl.xticks(pxtc, xtc)


#
# Restore default array print format
#
np.set_printoptions(suppress=False, precision=8)

pl.show()


if sf:
    pl.figure(fig1)
    pl.savefig("Distr_SNR_Lin_I-Cir_I_Diff_Stations.pdf", format='pdf')
    pl.figure(fig5)
    pl.savefig("Distr_SNR_Lin_I-Cir_I_Diff.pdf", format='pdf')








    
