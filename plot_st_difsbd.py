help_text = '''

plot_st_diffsbd.py: Plot histograms of differences between the sbdelay values
                    before and after PolConvert, for individual stations
                    and for all the stations.
'''

import sys
import pickle
import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm

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
# Gather SBD data into nsts station bins 
#
stsbd = {} # Dict for stationwise SBD data: stsbd['X'] 
stbls = {} # Dict for baselines including a station and their point numbers

for sta in ststr:
    
    tim = np.empty(0, dtype=float)   # Time for a particular station
    sbd_l = np.empty(0, dtype=float) # Lin SBD for a particular station
    sbd_c = np.empty(0, dtype=float) # Cir SBD for a particular station
    bsl = []  # List of baselines that include a particular station "sta"
    bsnpl = []  # List of numbers of points in the baselines with station "sta"

    ndat_st = 0 # Number of points for baselines with a station "sta"
    for bl in bls:   # Loop over the baselines
        if sta in bl:
            tim0 = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
            tim0 = tim0 - tim0[0]

            sbd0_l = np.array(idx3819l_1[bl]['I']['sbdelay'])[istart:]
            sbd0_c = np.array(idx3819c_1[bl]['I']['sbdelay'])[istart:]
            
            #
            # Subtract SBD means
            #
            sbd0_l = sbd0_l - sbd0_l.mean()
            sbd0_c = sbd0_c - sbd0_c.mean()

            tim = np.append(tim, tim0)
            sbd_l = np.append(sbd_l, sbd0_l)
            sbd_c= np.append(sbd_c, sbd0_c)
            
            ntim = len(tim0)
            ndat_st = ndat_st + ntim
            bsl.append(bl)
            bsnpl.append(ntim)
    print("'", sta, "': ", ndat_st) 
    #
    # Differences Lin-Cir for baselines with a particular station sta
    #
    dsbd = sbd_l - sbd_c
    stsbd[sta] = dsbd*1e6     # Convert us to ps
    stbls[sta] = [bsl, bsnpl]


nbin = 21   # Histogram bins
    
fig1 = pl.figure(figsize=(8, 10))
    
#
# Plot SBD histograms for the baselines including station "sta"
#
hw = 12  # Histogram width: +- hw

ist = 0   # Baseline number starting from 0
for sta in ststr:
    iplt = ist + 1  # Subplot number
    pl.figure(fig1)
    pl.subplot(3, 2, iplt)
    pl.hist(stsbd[sta], nbin, color='green')
    pl.xlabel("ps")
    pl.xlim(-300, 300)
    pl.grid(1)    
    ax = pl.gca()
    pl.text(.03, .92, "Station: "+sta, transform=ax.transAxes, fontsize=12)
    pl.text(.03, .84, "Bls: ", transform=ax.transAxes, fontsize=9)
    pl.text(.12, .84, ', '.join(stbls[sta][0]), transform=ax.transAxes, \
            fontsize=9)
    ist = ist + 1

    #
    # Testing the H0 hypothesis of stsbd[sta] normal distribution: FAILS!
    #
    ni, bedges = np.histogram(stsbd[sta], nbin) # 21 bin

    # ni = ni[7:15]
    # bedges = bedges[7:16]

    N = np.sum(ni)
    binwd = bedges[1] - bedges[0]             # Bin width
    xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
    hmean = np.sum(xi*ni)/N               # Sample mean
    sig2 = np.sum(xi**2*ni/N - hmean**2)  # Sample variance sigma^2
    sig = np.sqrt(sig2)                   # Standard deviation sigma
    #
    # Fit a normal distribution to the histogram and to the stsbd[sta] data
    #
    zi = (xi - hmean)/sig                 # Standardized xi
    fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
    fni = binwd*N*fnorm              # Theoretical frequencies
    mu, stdev = norm.fit(stsbd[sta]) # Fit a normal dist. to the stsbd[sta] data
    #
    # 
    #
    idm = np.where(abs(stsbd[sta]) < stdev)[0]    # stsbd[sta] within +-stdev
    pmstd = len(idm)/N*100  # Percent of stsbd[sta] within +-stdev
    pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
    print("%s dmdb: mu = %f,    std = %f" % (sta, mu, stdev))
    print("%s dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % \
          (sta, pmstd, pmstd_norm))

    #
    # Pearson's X^2
    #
    chi2obs = np.sum((ni - fni)**2/fni) # !!!!!!!! HUGE !!!!!!!!

    #
    # Smooth normal approximations 
    #
    x1 = np.linspace(-10, 10, 101)
    f2 = norm.pdf(x1, mu, stdev)*binwd*N

    pl.plot(x1, f2, 'k')

    stdh = 0.7*pl.ylim()[1]  # Height of std line
    pl.plot([-stdev, -stdev], [0, stdh], 'r-.')   # , lw=0.8)
    pl.plot([stdev, stdev], [0, stdh], 'r-.')     # , lw=0.8)
    pl.plot(xi, fni, 'b-')
    pl.plot(xi, fni, 'r.')

    # pl.xlim(-12,12)
    
    ax = pl.gca()
    pl.text(.6, .92, "Within $\pm$std: %5.2f%%" % pmstd, \
            transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .85, "For Normal: 68.27%", transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .78, "$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .71, "$\mu$=%.4f, $\sigma$=%5.2f" % (mu, stdev), \
            transform=ax.transAxes, fontsize=9)
    #
    # X ticks
    pxtc = -20 + 5*np.arange(9, dtype=float)
    pxtc = np.insert(pxtc, 4, -stdev)
    pxtc = np.insert(pxtc, 6, stdev)

    xtc = list(np.int64(pxtc))
    xtc[4] = r"$-\sigma$"
    xtc[6] = r"$+\sigma$"

    # Cut the ends to +-hw
    ipx = np.where(abs(pxtc) < 12)[0]
    pxtc = pxtc[ipx]
    xtc1 = []
    for i in ipx:
        xtc1.append(xtc[i])
    xtc = xtc1

    # pl.xticks(pxtc, xtc)

    # pl.xlim(-hw,+hw)
    pl.xlim(-300, +300)

fig1.text(0.2, 0.97, "Differences SBD Lin_I-Cir_I Distributions for Stations", \
          fontsize=12)
fig1.tight_layout(rect=(0,0,1, 0.95))



# sys.exit(0)

    
#
# Get and plot SBD for all the baselines 
#
rmse_sbd = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for SBD

dsbd = []  # Differences of SBD for all baselines
sbd_all_l = []  # SBD for all baselines, lin pol
sbd_all_c = []  # SBD for all baselines, cir pol

ibl = 0   # Baseline number starting from 0
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    sbd_l = np.array(idx3819l_1[bl]['I']['sbdelay'])[istart:]
    sbd_c = np.array(idx3819c_1[bl]['I']['sbdelay'])[istart:]

    sbd_all_l.extend(sbd_l) # Add SBD differences to list, lin pol
    sbd_all_c.extend(sbd_c) # Add SBD differences to list, cir pol
    
    dsbd_bl = np.zeros_like(tim)  # Differences of SBD for current baseline

    #
    # Subtract SBD means
    #
    sbd0_l = sbd_l - sbd_l.mean()
    sbd0_c = sbd_c - sbd_c.mean()
    
    print("'%s': abs(1e6*sbd0_l).mean() = %.2f,\t "
          "abs(1e6*sbd0_c).mean() = %.2f"% \
          (bl, abs(1e6*sbd0_l).mean(), abs(1e6*sbd0_c).mean()))
    
    #
    # Root mean square error (RMSE)
    #
    dsbd_bl = sbd0_l - sbd0_c
    dsbd.extend(dsbd_bl) # Add SBD differences to list
    npt = len(tim)   # Number of points for current baseline
    rmse_sbd[ibl] = np.sqrt(np.sum(dsbd_bl**2)/nbls)
    
    ibl = ibl + 1

dsbd = np.array(dsbd, dtype=float)*1e6 # Convert SBD from micro- to picoseconds
sbd_all_l = np.array(sbd_all_l, dtype=float)*1e6
sbd_all_c = np.array(sbd_all_c, dtype=float)*1e6

print("All baselines: abs(sbd_all_l).mean() = %.2f (ps),\t "
          "abs(sbd_all_c).mean() = %.2f (ps)" % \
          (abs(sbd_all_l).mean(), abs(sbd_all_c).mean()))
print("All baselines: dsbd min and max: ", dsbd.min(), dsbd.max())


fig5 = pl.figure()

pl.figure(fig5);
pl.hist(dsbd, nbin, color = "g", ec="k"); pl.grid(1)
pl.xlabel("ps")
# pl.xlim(-21, 21)
fig5.text(0.15, 0.95, "Differences SBD Lin_I-Cir_I Distributions " \
          "for All Baselines", \
          fontsize=12)
fig5.tight_layout(rect=(0,0,1, 0.95))


#
# Testing the H0 hypothesis or dsbd normal distribution: FAILS!
#
ni, bedges = np.histogram(dsbd, nbin) # 21 bin

# ni = ni[7:15]
# bedges = bedges[7:16]

N = np.sum(ni)
binwd = bedges[1] - bedges[0]             # Bin width
xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
hmean = np.sum(xi*ni)/N               # Sample mean
sig2 = np.sum(xi**2*ni/N - hmean**2)  # Sample variance sigma^2
sig = np.sqrt(sig2)                   # Standard deviation sigma
#
# Fit a normal distribution to the histogram and to the whole dsbd data
#
zi = (xi - hmean)/sig                 # Standardized xi
fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
fni = binwd*N*fnorm              # Theoretical frequencies
mu, stdev = norm.fit(dsbd)  # Fit a normal distribution to the WHOLE dsbd data
#
# 
#
idm = np.where(abs(dsbd) < stdev)[0]    # Count of dsbd within +-stdev 
# idm1 = np.where(abs(dsbd) >= stdev)[0]  # Count of dsbd outside of +-stdev
pmstd = len(idm)/N*100  # Percent of dsbd within +-stdev
pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
# print("Histogram:  hmean = %f, sig = %f" % (hmean, sig))
print("dmdb: mu = %f,    stdev = %f" % (mu, stdev))
print("dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % (pmstd, pmstd_norm))

#
# Pearson's X^2
#
chi2obs = np.sum((ni - fni)**2/fni) # !!!!!!!! HUGE !!!!!!!!
# scipy.stats.chi2.ppf(1-alpha, 19) = 30.14

# mu, stdev = norm.fit(dsbd)  # Fit a normal distribution to the WHOLE dsbd data

#
# Smooth normal approximations 
#
x1 = np.linspace(-10, 10, 101)
#f1 = norm.pdf(x1, hmean, sig)*binwd*N
f2 = norm.pdf(x1, mu, stdev)*binwd*N


pl.figure(fig5)

# pl.plot(x1, f1, 'm')
pl.plot(x1, f2, 'k')


pl.figure(fig5);

stdh = 0.85*pl.ylim()[1]  # Height of std line
pl.plot([-stdev, -stdev], [0, stdh], 'r-.')
pl.plot([stdev, stdev], [0, stdh], 'r-.')
pl.plot(xi, fni, 'b-')
pl.plot(xi, fni, 'ro')

ax = pl.gca()
pl.text(.04, .95, "Within $\pm$std: %5.2f%%" % pmstd, transform=ax.transAxes, \
        fontsize=10)
pl.text(.04, .90, "For Normal: 68.27%", transform=ax.transAxes, \
        fontsize=10)
pl.text(.75, .95, "$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
        fontsize=10)
pl.text(.75, .90, "$\mu$=%.4f, $\sigma$=%5.2f" % (mu, stdev), \
        transform=ax.transAxes, fontsize=10)
#
# X ticks
# pxtc = -20 + 5*np.arange(9, dtype=float)
# pxtc = np.insert(pxtc, 4, -stdev)
# pxtc = np.insert(pxtc, 6, stdev)

# xtc = list(np.int64(pxtc))
# xtc[4] = r"$-\sigma$"
# xtc[6] = r"$+\sigma$"

# pl.xticks(pxtc, xtc)

pl.show()





    
