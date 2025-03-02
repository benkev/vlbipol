help_text = '''

plot_st_difhist_2187_grp.py: Plot histograms of residuals between
                    the pseudo-Stokes
                    parameter values before and after PolConvert,
                    for the baselines involving individual stations
                    and for the baselines of all the stations.

                    This script differs from plot_st_difhist_2187.py in that
                    the statistics is applied to the histograms after their
                    sparse tails have been grouped.
'''

import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_st_difhist_2187.py <parname> [save], ")
    print("       where <parname> is either MBD or SBD or SNR.")
    print("       save (optional): save figures in pdf format.")
    sys.exit(0)
    
parname = sys.argv[1]
parname = parname.upper()
if parname != 'MBD' and parname != 'SBD' and parname != 'SNR':
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
import matplotlib.ticker as mticker
from scipy.stats import norm, chi2

from group_tails import find_tail_bounds, group_tails

pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
print("pl.isinteractive() -> ", pl.isinteractive())

snr_floor = 30    # Data with SNR < snr_floor will be discarded.
n_sig = 6         # Data deviating beyond n_sig*std will be discarded.
nbin_ini = 41     # Initial number of histogram bins (before tail grouping)
    

#
# Unpickle it:
#
with open('idx2187lI.pkl', 'rb') as finp:
    idxl = pickle.load(finp)

with open('idx2187cI.pkl', 'rb') as finp:
    idxc = pickle.load(finp)

if parname == 'MBD' or parname == 'SBD':
    ps = "(ps)"
else: # if parname == 'SNR':
    ps = ""
    


# bls = list(idxl.keys())   # List of baselines
# bls.sort()                      # Lexigraphically sorted baselines
# nbls = len(bls)

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
sttim = {} # Dict for stationwise time data: sttim['X'] 

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

if parname == 'MBD':
    par = 'mbdelay'
elif parname == 'SBD':
    par = 'sbdelay'
else:
    par = 'snr'
    

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

    #
    # Merge the data of the baselines containing the station sta
    #
    for bl in bls:   # Loop over the baselines

        if sta in bl:  

            snrbl_l = np.array(idxl[bl]['I']['snr'])
            snrbl_c = np.array(idxc[bl]['I']['snr'])
            timbl = np.array(idxl[bl]['I']['time']) / 3600 # Seconds to hours
            timbl0 = timbl[0]  # Save the starttime in case snr[0] <= snr_floor
            
            #
            # Exclude data with SNR <= snr_floor
            #
            isnr_floorl = np.where(snrbl_l <= snr_floor)[0]
            isnr_floorc = np.where(snrbl_c <= snr_floor)[0]
            # Merge lin & cir indices where snr <= snr_floor
            isnr_floor = np.concatenate((isnr_floorl, isnr_floorc)) 
            isnr_floor = np.unique(isnr_floor)
            isnr_floor.sort()                 # Sort the array in-place

            timbl = np.delete(timbl, isnr_floor)     # Exclude low SNR data
            timbl = timbl - timbl0  # Count time from the session start
            snrbl_l = np.delete(snrbl_l, isnr_floor) # Exclude low SNR data
            snrbl_c = np.delete(snrbl_c, isnr_floor) # Exclude low SNR data

            # Average of Lin and Cir means of SNR
            snrbl_a = (abs(snrbl_l.mean()) + abs(snrbl_c.mean()))/2


            if parname == 'MBD' or parname == 'SBD':
                parbl_l = np.array(idxl[bl]['I'][par])*1e6 # us to ps
                parbl_c = np.array(idxc[bl]['I'][par])*1e6 # us to ps
                parbl_l =  np.delete(parbl_l, isnr_floor) # Exclude low SNR data
                parbl_c =  np.delete(parbl_c, isnr_floor) # Exclude low SNR data
            else: # if parname == 'SNR':
                parbl_l = np.copy(snrbl_l)
                parbl_c = np.copy(snrbl_c)

            #
            # Exclude points beyond +-n_sig*std
            #
            parbl0_l = parbl_l - parbl_l.mean()  # Subtract par means, lin pol
            parbl0_c = parbl_c - parbl_c.mean()  # Subtract par means, cir pol

            sl = n_sig*np.std(parbl0_l)
            sc = n_sig*np.std(parbl0_c)
            isl = np.where((parbl0_l > sl) | (parbl0_l < -sl))[0]
            isc = np.where((parbl0_c > sc) | (parbl0_c < -sc))[0]    
            isg = np.concatenate((isl, isc))
            isg = np.unique(isg)
            isg.sort()            # Indices where |param| > n_sig*st
            
            # Exclude data with |param| > n_sig*std
            timbl = np.delete(timbl, isg)
            parbl0_l = np.delete(parbl0_l, isg)
            parbl0_c = np.delete(parbl0_c, isg)
            parbl_l = np.delete(parbl_l, isg)
            parbl_c = np.delete(parbl_c, isg)
    
            # Average of Lin and Cir means of the parameter (MBD or SDD)
            parbl_a = (abs(parbl_l.mean()) + abs(parbl_c.mean()))/2


            tim = np.append(tim, timbl)
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
    sttim[sta] = tim
    stbls[sta] = [bsl, bsnpl]
    
for sta in ststr:
    print("'%s': len(stpar['%s'])) = %d" % (sta, sta, len(stpar[sta])))
    print("   : np.where(np.isnan(stpar['%s'])) = " % sta,
                np.where(np.isnan(stpar[sta]))[0])

#
# Set default array print format: as fixed-point only and as short as possible
#
#np.set_printoptions(suppress=True, precision=1)
np.set_printoptions(precision=6, legacy='1.25')

fig1 = pl.figure(figsize=(8, 10))
    
#
# Plot histograms of par residuals for the baselines including station "sta"
#
ist = 0   # Baseline number starting from 0
for sta in ststr:

    ni_ini, bedges = np.histogram(stpar[sta], nbin_ini) # 21 bin

    N = np.sum(ni_ini)
    binwd = bedges[1] - bedges[0]             # Bin width
    #
    # Histogram initial parameters with sparse tails
    #
    xi_ini = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals
    hmean_ini = np.sum(xi_ini*ni_ini)/N               # Sample mean
    #sig2_ini = np.sum(xi_ini**2*ni_ini/N - hmean_ini**2)  # W R O N G ! ! !
    # Sample variance: sig^2 = M(X^2) - (M(X))^2
    sig2_ini = np.sum(xi_ini**2*ni_ini/N) - hmean_ini**2  # Sample variance
    sig_ini = np.sqrt(sig2_ini)                   # Standard deviation sigma
    #
    # Fit a normal distribution to the histogram
    #
    zi_ini = (xi_ini - hmean_ini)/sig_ini             # Standardized xi
    # Standard normal PDF
    fnorm_ini = (1/(sig_ini*np.sqrt(2*np.pi)))*np.exp(-zi_ini**2/2)
    fni_ini = binwd*N*fnorm_ini              # Theoretical frequencies

    #
    # Group left-tail and right-tail bins with sparse data.
    #

    l_idx, r_idx = find_tail_bounds(ni_ini)
    ni =  group_tails(ni_ini,  (l_idx, r_idx)) 
    fni = group_tails(fni_ini, (l_idx, r_idx)) 
    xi =  np.copy(xi_ini[l_idx:r_idx]) 
    nbin = len(ni)
    
    #
    # Histogram parameters with the grouped tails
    #
    hmean = np.sum(xi*ni)/N               # Sample mean
    sig2 = np.sum(xi**2*ni/N) - hmean**2  # Sample variance sigma^2
    sig = np.sqrt(sig2)                   # Standard deviation sigma
    #
    # Fit a normal distribution to the histogram
    #
    zi = (xi - hmean)/sig                 # Standardized xi
    fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
    fni = binwd*N*fnorm              # Theoretical frequencies
    # Fit a normal distribution to the stpar[sta] data
    mu, stdev = norm.fit(stpar[sta])
    #
    # 
    #
    idm = np.where(abs(stpar[sta]) < stdev)[0]    # stpar[sta] within +-stdev
    pmstd = len(idm)/N*100  # Percent of stpar[sta] within +-stdev
    pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
    print(r"%s dmdb: mu = %f,    std = %f" % (sta, mu, stdev))
    print(r"%s dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % \
          (sta, pmstd, pmstd_norm))

    print(r"st. %s: zi[0], zi[%d] = %f, %f" % (sta, len(zi)-1, zi[0], zi[-1]))
    
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

    print("Station %s:" % sta)
    print('Original binning with sparse tails (%d bins):' % nbin_ini)
    print('ni_ini:  ', ni_ini)
    #print('fni: ', fni_ini)
    print('Sparse tails grouped: (%d bins, [%d:%d])' % (nbin, l_idx, r_idx))
    print('ni:  ', ni)
    #print('fni: ', fni)
    print()
    print("%s nbin = %d, chi2obs = %.1f, chi2cr = %.1f chi2obs/chi2cr = %.1f" %
          (sta, nbin, chi2obs, chi2cr, q_chi2))

    #
    # Plot the grouped histogram
    #
    
    iplt = ist + 1  # Subplot number
    pl.figure(fig1)
    pl.subplot(4, 2, iplt)

    # bws = binwd*np.ones(nbin) # Bin widths
        
    pl.bar(xi, ni, width=0.6*binwd, edgecolor='black', color='g',
           align='center')

    #
    # Slightly raise the histogram over the zero level
    #
    hbot, htop = pl.ylim()
    yrng = htop - hbot
    pl.ylim(hbot-0.015*yrng, htop)
    
    if  parname == 'MBD':
        pl.xlabel("ps")
    elif parname == 'SBD':
        pl.xlabel("ps")
    else: # if parname == 'SNR':
        pass

    pl.grid(1)    
    ax = pl.gca()
    pl.text(.03, .92, "Station: "+sta, transform=ax.transAxes, fontsize=12)
    pl.text(.03, .84, "Bls: ", transform=ax.transAxes, fontsize=9)
    pl.text(.12, .84, ', '.join(stbls[sta][0]), transform=ax.transAxes, \
            fontsize=9)

    #
    # Smooth normal approximations 
    #
    # if  parname == 'MBD':
    #     x1 = np.linspace(-11, 11, 101)
    # elif parname == 'SBD':
    #     x1 = np.linspace(-200, 200, 101)
    # else: # if parname == 'SNR':
    #     x1 = np.linspace(-100, 100, 101)

    ax0 = np.abs(xi[0]); ax1 = np.abs(xi[-1])
    xna = ax0 if ax0 > ax1 else ax1
    
    x1 = np.linspace(-xna, +xna, 1001)
    
    f2 = norm.pdf(x1, mu, stdev)*binwd*N

    pl.plot(x1, f2, 'r')    # Plot smooth normal approximations 

    stdh = 0.6*pl.ylim()[1]  # Height of std line
    pl.plot([-stdev, -stdev], [0, stdh], 'r-.')   # , lw=0.8)
    pl.plot([stdev, stdev], [0, stdh], 'r-.')     # , lw=0.8)
    
    pl.plot(xi, fni, 'b.')

   
    ax = pl.gca()

    pc = r'\%'
    
    pl.text(.6, .92, r"Within $\pm$std: %5.2f%s" % (pmstd, pc), \
            transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .85, r"For Normal: 68.27%s" % pc, transform=ax.transAxes, \
            fontsize=9)
    pl.text(.6, .78, r"$\mu$=%.4f, $\sigma$=%5.2f" % (mu, stdev), \
            transform=ax.transAxes, fontsize=9)
    pl.text(.67, .71, r"N=%3d; bins: %2d" % (N, nbin), \
            transform=ax.transAxes, fontsize=9)
    pl.text(.67, .64, r"df=%2d-2-1=%2d" % (nbin, deg_fr), \
            transform=ax.transAxes, fontsize=9)
    if chi2obs < 1000:
        pl.text(.67, .57, r"$\chi^2$=%6.2f" % chi2obs, transform=ax.transAxes, \
                fontsize=9)
    else:
        pl.text(.67, .57, r"$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
                fontsize=9)
    if q_chi2 <= 1:
        pl.text(.66, .50, r"$\chi^2 \leq \chi^2_{cr}$=%.2f" % chi2cr, \
                transform=ax.transAxes, fontsize=11, c='red')        
    elif q_chi2 < 5:
        pl.text(.67, .50, r"$\chi^2 > \chi^2_{cr}$=%.2f" % chi2cr, \
                transform=ax.transAxes, fontsize=9)
    else:
        pl.text(.67, .50, r"$\chi^2 \gg \chi^2_{cr}$=%.2f" % chi2cr, \
                transform=ax.transAxes, fontsize=9)

    # if  parname == 'MBD' or parname == 'SBD':
    #     pl.text(.67, .41, r"$\overline{SNR}$: %.1f" % stsnr[sta], \
    #             transform=ax.transAxes, fontsize=9)
    #     pl.text(.67, .36, r"r_corr: %.6f" % str_corr[sta], \
    #             transform=ax.transAxes, fontsize=9)
    # else: # if  parname == 'SNR':
    #     pl.text(.67, .42, r"r_corr: %.6f" % str_corr[sta], \
    #             transform=ax.transAxes, fontsize=9)
    pl.text(.67, .41, r"$\overline{\mathrm{SNR}}$: %.1f" % stsnr[sta], \
            transform=ax.transAxes, fontsize=9)
    pl.text(.67, .36, r"r_corr: %.6f" % str_corr[sta], \
            transform=ax.transAxes, fontsize=9)

    #
    # X ticks
    #
    if  parname == 'MBD':
        pxtc = -20 + 5*np.arange(9, dtype=float) # X ticks positions
    
        print("pxtc     = ", pxtc)

        i_mstd = 4    # Position of -stdev tick
        i_pstd = 6    # Position of +stdev tick
        pl.xlabel("ps")
        #hw = 12            # Histogram width
        hw = 23            # Histogram width
        
    elif parname == 'SBD':
        pxtc = -300 + 100*np.arange(7, dtype=float) # X ticks positions
        i_mstd = 3    # Position of -stdev
        i_pstd = 5    # Position of +stdev
        pl.xlabel("ps")
        hw = 300           # Histogram width

    else: # if parname == 'SNR':
        pxtc = -100 + 50*np.arange(7, dtype=float) # X ticks positions
        i_mstd = 2    # Position of -stdev
        i_pstd = 4    # Position of +stdev
        hw = 120           # Histogram width

    pxtc = np.insert(pxtc, i_mstd, -stdev)
    pxtc = np.insert(pxtc, i_pstd, +stdev)
    xtc = list(np.int64(pxtc))
    xtc[i_mstd] = r"$-\sigma$"    # Tick label for -stdev
    xtc[i_pstd] = r"$+\sigma$"    # Tick label for +stdev

    #pl.xticks(pxtc, xtc)
    #pl.xlim(-hw,+hw)
    

    # ax.set_xticks(list(ax.get_xticks()) + [-stdev, +stdev])
    # xtls = [t.get_text() for t in ax.get_xticklabels()]
    # xtls[-2:] = [r'$-\sigma$', r'$+\sigma$']
    # ax.set_xticklabels(xtls)
    yl = ax.get_ylim(); yrng = yl[1] - yl[0]; yp = - 0.1*yrng
    xl = ax.get_xlim(); xrng = xl[1] - xl[0];
    pl.text(-stdev-0.05*xrng, stdh+0.05*yrng, r'$-\sigma$',
            fontsize=12, color='red')
    pl.text(+stdev-0.03*xrng, stdh+0.05*yrng, r'$+\sigma$',
            fontsize=12, color='red')



    # xp = -stdev - 0.05*yrng
    # ax.text(xp, yp, r'$-\sigma$', fontsize=12, color='red')
    # xp = stdev - 0.05*yrng
    # ax.text(xp, yp, r'$+\sigma$', fontsize=12, color='red')
    print("bedges: %f ~ %f, xl: %f ~ %f" % (bedges[0], bedges[-1],
                                            xl[0], xl[1]))

    ist = ist + 1

    # sys.exit(0)
    

fig1.text(0.2, 0.96, \
          "VO2187: Distributions of %s Residuals Lin_I-Cir_I for Stations"
          % parname, fontsize=12)
fig1.text(0.65, 0.15, r"$\mathrm{SNR} > %d$" % snr_floor, fontsize=16)
fig1.text(0.65, 0.10, r"$|\mathrm{%s}| < 6\sigma$" % parname, fontsize=16)

fig1.tight_layout(rect=(0,0,1, 0.95))



# ================= HIST FOR ALL STATIONS ===================================

# nbin_ini = 21

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

# for bl in bls:   # Loop over the baselines

#     snrbl_l = np.array(idxl[bl]['I']['snr'])
#     snrbl_c = np.array(idxc[bl]['I']['snr'])

#     #
#     # Exclude data with SNR <= snr_floor
#     #
#     isnr_floorl = np.where(snrbl_l <= snr_floor)[0]
#     isnr_floorc = np.where(snrbl_c <= snr_floor)[0]
#     # Merge lin & cir indices where snr <= snr_floor
#     isnr_floor = np.concatenate((isnr_floorl, isnr_floorc)) 
#     isnr_floor = np.unique(isnr_floor)
#     isnr_floor.sort()           # Indices where SNR <= snr_floor

#     timbl = np.array(idxl[bl]['I']['time']) / 3600 # Seconds to hours
#     timbl = np.delete(timbl, isnr_floor)     # Exclude low SNR data
#     timbl0 = timbl - timbl[0]  # Count time from the session start

#     snrbl_l = np.delete(snrbl_l, isnr_floor) # Exclude low SNR data
#     snrbl_c = np.delete(snrbl_c, isnr_floor) # Exclude low SNR data

#     # Average of Lin and Cir means
#     snrbl_a = (abs(snrbl_l.mean()) + abs(snrbl_c.mean()))/2


#     if parname == 'MBD' or parname == 'SBD':
#         parbl_l = np.array(idxl[bl]['I'][par])*1e6 # us to ps
#         parbl_c = np.array(idxc[bl]['I'][par])*1e6 # us to ps
#         parbl_l =  np.delete(parbl_l, isnr_floor) # Exclude low SNR data
#         parbl_c =  np.delete(parbl_c, isnr_floor) # Exclude low SNR data
#     else: # if parname == 'SNR':
#         parbl_l = np.copy(snrbl_l)
#         parbl_c = np.copy(snrbl_c)


#     # Average of Lin and Cir means of the parameter (MBD or SDD)
#     parbl_a = (abs(parbl_l.mean()) + abs(parbl_c.mean()))/2

#     parbl0_l = parbl_l - parbl_l.mean()  # Subtract par means, lin pol
#     parbl0_c = parbl_c - parbl_c.mean()  # Subtract par means, cir pol

#     tim = np.append(tim, timbl0)
#     par_l = np.append(par_l, parbl_l)
#     par_c = np.append(par_c, parbl_c)
#     par0_l = np.append(par0_l, parbl0_l) # Mean subtracted
#     par0_c = np.append(par0_c, parbl0_c) # Mean subtracted
#     snr_a = np.append(snr_a, snrbl_a)


for bl in bls:   # Loop over the baselines
    
    snrbl_l = np.array(idxl[bl]['I']['snr'])
    snrbl_c = np.array(idxc[bl]['I']['snr'])
    timbl = np.array(idxl[bl]['I']['time']) / 3600 # Seconds to hours
    timbl0 = timbl[0]  # Save the starttime in case snr[0] <= snr_floor

    #
    # Exclude data with SNR <= snr_floor
    #
    isnr_floorl = np.where(snrbl_l <= snr_floor)[0]
    isnr_floorc = np.where(snrbl_c <= snr_floor)[0]
    # Merge lin & cir indices where snr <= snr_floor
    isnr_floor = np.concatenate((isnr_floorl, isnr_floorc)) 
    isnr_floor = np.unique(isnr_floor)
    isnr_floor.sort()                 # Sort the array in-place

    timbl = np.delete(timbl, isnr_floor)     # Exclude low SNR data
    timbl = timbl - timbl0  # Count time from the session start
    snrbl_l = np.delete(snrbl_l, isnr_floor) # Exclude low SNR data
    snrbl_c = np.delete(snrbl_c, isnr_floor) # Exclude low SNR data

    # Average of Lin and Cir means of SNR
    snrbl_a = (abs(snrbl_l.mean()) + abs(snrbl_c.mean()))/2


    if parname == 'MBD' or parname == 'SBD':
        parbl_l = np.array(idxl[bl]['I'][par])*1e6 # us to ps
        parbl_c = np.array(idxc[bl]['I'][par])*1e6 # us to ps
        parbl_l =  np.delete(parbl_l, isnr_floor) # Exclude low SNR data
        parbl_c =  np.delete(parbl_c, isnr_floor) # Exclude low SNR data
    else: # if parname == 'SNR':
        parbl_l = np.copy(snrbl_l)
        parbl_c = np.copy(snrbl_c)

    #
    # Exclude points beyond +-n_sig*std
    #
    parbl0_l = parbl_l - parbl_l.mean()  # Subtract par means, lin pol
    parbl0_c = parbl_c - parbl_c.mean()  # Subtract par means, cir pol

    sl = n_sig*np.std(parbl0_l)
    sc = n_sig*np.std(parbl0_c)
    isl = np.where((parbl0_l > sl) | (parbl0_l < -sl))[0]
    isc = np.where((parbl0_c > sc) | (parbl0_c < -sc))[0]    
    isg = np.concatenate((isl, isc))
    isg = np.unique(isg)
    isg.sort()            # Indices where |param| > n_sig*st

    # Exclude data with |param| > n_sig*std
    timbl = np.delete(timbl, isg)
    parbl0_l = np.delete(parbl0_l, isg)
    parbl0_c = np.delete(parbl0_c, isg)
    parbl_l = np.delete(parbl_l, isg)
    parbl_c = np.delete(parbl_c, isg)

    # Average of Lin and Cir means of the parameter (MBD or SDD)
    parbl_a = (abs(parbl_l.mean()) + abs(parbl_c.mean()))/2


    tim = np.append(tim, timbl)
    par_l = np.append(par_l, parbl_l)
    par_c = np.append(par_c, parbl_c)
    par0_l = np.append(par0_l, parbl0_l) # Mean subtracted
    par0_c = np.append(par0_c, parbl0_c) # Mean subtracted
    snr_a = np.append(snr_a, snrbl_a)





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

print(r"All baselines: abs(par_l).mean() = %.2f (ps),\t "
          "abs(par_c).mean() = %.2f (ps)" % \
          (abs(par_l).mean(), abs(par_c).mean()))
print(r"All baselines: dpar min and max: ", dpar.min(), dpar.max())

fig2 = pl.figure()

pl.figure(fig2);


# pl.hist(dpar, nbin_ini, color = "g", ec="k"); pl.grid(1)

# #
# # Slightly raise the histogram over the zero level
# #
# hbot, htop = pl.ylim()
# yrng = htop - hbot
# pl.ylim(hbot-0.015*yrng, htop)


if  parname == 'MBD':
    pl.xlabel("ps")
elif parname == 'SBD':
    pl.xlabel("ps")
else: # if parname == 'SNR':
    pass

        
fig2.text(0.15, 0.95, "VO2187: Distribution of %s Residuals Lin_I-Cir_I " \
          "for All Baselines" % parname, \
          fontsize=12)
fig2.tight_layout(rect=(0,0,1, 0.95))


# ni, bedges = np.histogram(dpar, nbin_ini) # 21 bin

# N = np.sum(ni)
# binwd = bedges[1] - bedges[0]             # Bin width
# xi = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals    
# hmean = np.sum(xi*ni)/N               # Sample mean
# sig2 = np.sum(xi**2*ni/N) - hmean**2  # Sample variance sigma^2
# sig = np.sqrt(sig2)                   # Standard deviation sigma
# #
# # Fit a normal distribution to the histogram and to the whole dpar data
# #
# zi = (xi - hmean)/sig                 # Standardized xi
# fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
# fni = binwd*N*fnorm              # Theoretical frequencies
# mu, stdev = norm.fit(dpar)  # Fit a normal distribution to the WHOLE dpar data
# #
# # 
# #
# idm = np.where(abs(dpar) < stdev)[0]    # Count of dpar within +-stdev 
# # idm1 = np.where(abs(dpar) >= stdev)[0]  # Count of dpar outside of +-stdev
# pmstd = len(idm)/N*100  # Percent of dpar within +-stdev
# pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
# # print(r"Histogram:  hmean = %f, sig = %f" % (hmean, sig))
# print(r"dmdb: mu = %f,    stdev = %f" % (mu, stdev))
# print(r"dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % (pmstd, pmstd_norm))


# #
# # Group left-tail and right-tail bins with sparse data.
# #

# ni_ini = np.copy(ni)     # Save old freqs with sparse tails
# fni_ini = np.copy(fni)    # Save old freqs with sparse tails

# ni, fni = group_tails(ni_ini, fni_ini) 

# nbin = len(ni)





ni_ini, bedges = np.histogram(dpar, nbin_ini) # 21 bin

N = np.sum(ni)
binwd = bedges[1] - bedges[0]             # Bin width

#
# Histogram initial parameters with sparse tails
#
xi_ini = (bedges[1:] + bedges[:-1])/2          # Middles of the intervals
hmean_ini = np.sum(xi_ini*ni_ini)/N               # Sample mean
# Sample variance: sig^2 = M(X^2) - (M(X))^2
sig2_ini = np.sum(xi_ini**2*ni_ini/N) - hmean_ini**2  # Sample variance
sig_ini = np.sqrt(sig2_ini)                   # Standard deviation sigma
#
# Fit a normal distribution to the histogram
#
zi_ini = (xi_ini - hmean_ini)/sig_ini             # Standardized xi
# Standard normal PDF
fnorm_ini = (1/(sig_ini*np.sqrt(2*np.pi)))*np.exp(-zi_ini**2/2)
fni_ini = binwd*N*fnorm_ini              # Theoretical frequencies

#
# Group left-tail and right-tail bins with sparse data.
#

l_idx, r_idx = find_tail_bounds(ni_ini)
ni =  group_tails(ni_ini,  (l_idx, r_idx)) 
fni = group_tails(fni_ini, (l_idx, r_idx)) 
xi =  np.copy(xi_ini[l_idx:r_idx]) 
nbin = len(ni)

#
# Histogram parameters with the grouped tails
#
hmean = np.sum(xi*ni)/N               # Sample mean
# Sample variance: sig^2 = M(X^2) - (M(X))^2
sig2 = np.sum(xi**2*ni/N) - hmean**2  # Sample variance sigma^2
sig = np.sqrt(sig2)                   # Standard deviation sigma
#
# Fit a normal distribution to the histogram
#
zi = (xi - hmean)/sig                 # Standardized xi
fnorm = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-zi**2/2)   # Standard normal PDF
fni = binwd*N*fnorm              # Theoretical frequencies
# Fit a normal distribution to the stpar[sta] data
mu, stdev = norm.fit(stpar[sta])
#
# 
#
idm = np.where(abs(stpar[sta]) < stdev)[0]    # stpar[sta] within +-stdev
pmstd = len(idm)/N*100  # Percent of stpar[sta] within +-stdev
pmstd_norm = 68.27         # Percent of normal data within +-std: 68.27%
print(r"%s dmdb: mu = %f,    std = %f" % (sta, mu, stdev))
print(r"%s dmdb: %5.2f%% within +-std;  %5.2f%% for normal" % \
      (sta, pmstd, pmstd_norm))

print(r"st. %s: zi[0], zi[%d] = %f, %f" % (sta, len(zi)-1, zi[0], zi[-1]))


#
# Plot the grouped histogram
#

# bws = binwd*np.ones(nbin) # Bin widths

pl.bar(xi, ni, width=0.6*binwd, edgecolor='black', color='g', align='center')

#
# Slightly raise the histogram over the zero level
#
hbot, htop = pl.ylim()
yrng = htop - hbot
pl.ylim(hbot-0.015*yrng, htop)

if  parname == 'MBD':
    pl.xlabel("ps")
elif parname == 'SBD':
    pl.xlabel("ps")
else: # if parname == 'SNR':
    pass


pl.grid(1)    







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

print('All stations:')
# print('Original binning with sparse tails (%d bins):' % nbin_ini)
print('ni_ini:  ', ni_ini)
# print('fni: ', fni_ini)
# print('Sparse tails grouped: (%d bins)' % nbin)
print('Sparse tails grouped: (%d bins, [%d:%d])' % (nbin, l_idx, r_idx))
print('ni:  ', ni)
# print('fni: ', fni)
# print("chi2obs/chi2cr = %f" % q_chi2)
# print()
print("%s nbin = %d, chi2obs = %.1f, chi2cr = %.1f chi2obs/chi2cr = %.1f" %
      (sta, nbin, chi2obs, chi2cr, q_chi2))

#
# Smooth normal approximations 
#
ax0 = np.abs(xi[0]); ax1 = np.abs(xi[-1])
xna = ax0 if ax0 > ax1 else ax1

x1 = np.linspace(-xna, +xna, 1001)

f2 = norm.pdf(x1, mu, stdev)*binwd*N


pl.figure(fig2)

pl.plot(x1, f2, 'r')

stdh = 0.85*pl.ylim()[1]  # Height of std line
pl.plot([-stdev, -stdev], [0, stdh], 'r-.')
pl.plot([stdev, stdev], [0, stdh], 'r-.')

pl.plot(xi, fni, 'bo')

ax = pl.gca()
pl.text(.04, .95, r"Within $\pm$std: %5.2f\%%" % pmstd, transform=ax.transAxes, \
        fontsize=10)
pl.text(.04, .90, r"For Normal: 68.27\%%", transform=ax.transAxes, \
        fontsize=10)
pl.text(.75, .95, r"$\mu$=%.4f, $\sigma$=%5.2f" % (mu, stdev), \
        transform=ax.transAxes, fontsize=10)
pl.text(.75, .90, r"N = %3d; bins: %2d" % (N, nbin), \
        transform=ax.transAxes, fontsize=10)
pl.text(.75, .85, r"df = %2d-2-1=%2d" % (nbin, deg_fr), \
        transform=ax.transAxes, fontsize=10)
if chi2obs < 1000:
    pl.text(.75, .80, r"$\chi^2$=%6.2f" % chi2obs, transform=ax.transAxes, \
            fontsize=10)
else:
    pl.text(.75, .80, r"$\chi^2$=%9.2e" % chi2obs, transform=ax.transAxes, \
            fontsize=10)

if q_chi2 <= 1:
    pl.text(.75, .74, r"$\chi^2 \leq \chi^2_{cr}$=%.2f" % chi2cr, \
            transform=ax.transAxes, fontsize=11, c='red')        
elif q_chi2 < 5:
    pl.text(.75, .75, r"$\chi^2 > \chi^2_{cr}$ = %.2f" % chi2cr, \
            transform=ax.transAxes, fontsize=10)
else:
    pl.text(.75, .75, r"$\chi^2 \gg \chi^2_{cr}$ = %.2f" % chi2cr, \
            transform=ax.transAxes, fontsize=10)
    
#        ---------------   Include avg(SNR) in the SNR plots ? ----------
# if  parname == 'MBD' or parname == 'SBD':
#     pl.text(.75, .68, r"$\overline{\mathrm{SNR}}$: %.1f" % snr_avg, \
#             transform=ax.transAxes, fontsize=10)
#     pl.text(.75, .64, r"r_corr: %.6f" % r_corr, transform=ax.transAxes, \
#             fontsize=10)
# else: # if  parname == 'SNR':
#     pl.text(.75, .70, r"r_corr: %.6f" % r_corr, transform=ax.transAxes, \
#             fontsize=10)

pl.text(.75, .68, r"$\overline{\mathrm{SNR}}$: %.1f" % snr_avg, \
        transform=ax.transAxes, fontsize=10)
pl.text(.75, .64, r"r_corr: %.6f" % r_corr, transform=ax.transAxes, \
        fontsize=10)


#
# X ticks
#
if  parname == 'MBD':
    pxtc = -20 + 5*np.arange(9, dtype=float) # X ticks positions
    i_mstd = 4    # Position of -stdev tick
    i_pstd = 6    # Position of +stdev tick
    pl.xlabel("ps")
    #hw = 12            # Histogram width
    hw = 23            # Histogram width

elif parname == 'SBD':
    pxtc = -300 + 100*np.arange(7, dtype=float) # X ticks positions
    i_mstd = 3    # Position of -stdev
    i_pstd = 5    # Position of +stdev
    pl.xlabel("ps")
    hw = 300           # Histogram width
    
else: # if parname == 'SNR':
    pxtc = -100 + 50*np.arange(7, dtype=float) # X ticks positions
    i_mstd = 2    # Position of -stdev
    i_pstd = 3    # Position of +stdev
    hw = 120           # Histogram width

    
pxtc = np.insert(pxtc, i_mstd, -stdev)
pxtc = np.insert(pxtc, i_pstd, +stdev)
xtc = list(np.int64(pxtc))
xtc[i_mstd] = r"$-\sigma$"    # Tick label for -stdev
xtc[i_pstd] = r"$+\sigma$"    # Tick label for +stdev

# ax.set_xticks(list(ax.get_xticks()) + [-stdev, +stdev])
# xtls = [t.get_text() for t in ax.get_xticklabels()]
# xtls[-2:] = [r'$-\sigma$', r'$+\sigma$']
# ax.set_xticklabels(xtls)

# ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
# ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

yl = ax.get_ylim(); yrng = yl[1] - yl[0]; yp = - 0.1*yrng
xl = ax.get_xlim(); xrng = xl[1] - xl[0];
pl.text(-stdev-0.05*xrng, stdh+0.05*yrng, r'$-\sigma$',
        fontsize=12, color='red')
pl.text(+stdev-0.03*xrng, stdh+0.05*yrng, r'$+\sigma$',
        fontsize=12, color='red')



# xp = -stdev - 0.05*yrng
# ax.text(xp, yp, r'$-\sigma$', fontsize=12, color='red')
# xp = stdev - 0.05*yrng
# ax.text(xp, yp, r'$+\sigma$', fontsize=12, color='red')
print("bedges: %f ~ %f, xl: %f ~ %f" % (bedges[0], bedges[-1],
                                        xl[0], xl[1]))

#pl.xticks(pxtc, xtc)
#pl.xlim(-hw,+hw)

pl.show()

if sf:
    pl.figure(fig1)
    pl.savefig("VO2187_Distr_%s_Lin_I-Cir_I_Diff_Stations.pdf" % parname,
               format='pdf')
    pl.figure(fig2)
    pl.savefig("VO2187_Distr_%s_Lin_I-Cir_I_Diff.pdf" % parname, format='pdf')


#
# Restore default array print format
#
# np.set_printoptions(suppress=False, precision=8)







    
