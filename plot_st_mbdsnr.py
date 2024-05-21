import sys
import pickle
import numpy as np
import matplotlib.pyplot as pl

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
# Separate MBD and SNR sata into nsts station bins 
#
stmbd = {} # Dict for stationwise MBD data: stmbd['X'] 
stsnr = {} # Dict for stationwise SNR data: stmbd['X'] 
stbls = {} # Dict for baselines including a station and their point numbers

for sta in ststr:
    
    tim = np.empty(0, dtype=float)   # Time for a particular station
    mbd_l = np.empty(0, dtype=float) # Lin MBD for a particular station
    mbd_c = np.empty(0, dtype=float) # Cir MBD for a particular station
    snr_l = np.empty(0, dtype=float) # Lin SNR for a particular station
    snr_c = np.empty(0, dtype=float) # Cir SNR for a particular station
    bsl = []  # List of baselines that include a particular station "sta"
    bsnpl = []  # List of numbers of points in the baselines with station "sta"

    ndat_st = 0 # Number of points for baselines with a station "sta"
    for bl in bls:   # Loop over the baselines
        if sta in bl:
            tim0 = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
            tim0 = tim0 - tim0[0]

            mbd0_l = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:]
            mbd0_c = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:]
            snr0_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
            snr0_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
            
            #
            # Subtract MBD and SNR means
            #
            mbd0_l = mbd0_l - mbd0_l.mean()
            mbd0_c = mbd0_c - mbd0_c.mean()
            snr0_l = snr0_l - snr0_l.mean()
            snr0_c = snr0_c - snr0_c.mean()

            tim = np.append(tim, tim0)
            mbd_l = np.append(mbd_l, mbd0_l)
            mbd_c= np.append(mbd_c, mbd0_c)
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
    dmbd = mbd_l - mbd_c
    dsnr = snr_l - snr_c
    stmbd[sta] = dmbd*1e6     # Convert us to ps
    stsnr[sta] = dsnr
    stbls[sta] = [bsl, bsnpl]

fig1 = pl.figure(figsize=(8, 10))
    
#
# Plot MBD histograms for the baselines including station "sta"
#
ist = 0   # Baseline number starting from 0
for sta in ststr:
    iplt = ist + 1  # Subplot number
    pl.figure(fig1)
    pl.subplot(3, 2, iplt)
    pl.hist(stmbd[sta], 31, color='green')
    pl.xlabel("ps")
    pl.xlim(-21, 21)
    pl.grid(1)    
    ax = pl.gca()
    pl.text(.04, .92, "Station: "+sta, transform=ax.transAxes, fontsize=12)
    pl.text(.04, .84, "Bls: ", transform=ax.transAxes, fontsize=10)
    pl.text(.14, .84, ', '.join(stbls[sta][0]), transform=ax.transAxes, \
            fontsize=10)
    # for bl in stbls[sta][0]:
    #     pl.text(.10, .84, bl+" ", transform=ax.transAxes, fontsize=12)
#    pl.text(.1, .9, sta,  ha='left', va='top', transform=ax.transAxes)
    ist = ist + 1

fig1.text(0.3, 0.97, "MBD Lin_I-Cir_I Distributions for Stations", fontsize=12)
fig1.tight_layout(rect=(0,0,1, 0.95))

#fig1.text(0.07, 0.9, "E", fontsize=16)



sys.exit(0)

    

rmse_mbd = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for MBD
rmse_snr = np.zeros(nbls, dtype=float)  # Root mean square error (RMSE) for SNR

# fig2 = pl.figure(figsize=(8, 12))
# fig3 = pl.figure(figsize=(8, 12))
# fig4 = pl.figure(figsize=(8, 12))

dmbd = []  # Differences of MBD for all baselines
dsnr = []  # Differences of SNR for all baselines

ibl = 0   # Baseline number starting from 0
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    mbd_l = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:]
    mbd_c = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:]
    snr_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
    dmbd_bl = np.zeros_like(tim)  # Differences of MBD for current baseline
    dsnr_bl = np.zeros_like(tim)  # Differences of SNR for current baseline
    isplt = ibl + 1  # Subplot number

    #
    # Plot MBD
    #
    pl.figure(fig1)
    pl.subplot(5, 3, isplt)
    pl.plot(tim, mbd_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    #
    # Subtract MBD means
    #
    mbd0_l = mbd_l - mbd_l.mean()
    mbd0_c = mbd_c - mbd_c.mean()
    #
    # Subtract SNR means
    #
    snr0_l = snr_l - snr_l.mean()
    snr0_c = snr_c - snr_c.mean()
    
    #
    # Root mean square error (RMSE)
    #
    dmbd_bl = mbd0_l - mbd0_c
    dsnr_bl = snr0_l - snr0_c
    dmbd.extend(dmbd_bl) # Add MBD differences to list
    dsnr.extend(dsnr_bl) # Add SNR differences to list
    npt = len(tim)   # Number of points for current baseline
    rmse_mbd[ibl] = np.sqrt(np.sum(dmbd_bl**2)/nbls)
    rmse_snr[ibl] = np.sqrt(np.sum(dsnr_bl**2)/nbls)
    
    pl.figure(fig2)
    pl.subplot(5, 3, isplt)
    pl.plot(tim, mbd0_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd0_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    #
    # Plot SNR
    #
    pl.figure(fig3)
    pl.subplot(5, 3, isplt)
    pl.plot(tim, snr_l , label='Lin_I, '+bl)
    pl.plot(tim, snr_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    pl.figure(fig4)
    pl.subplot(5, 3, isplt)
    pl.plot(tim, snr0_l , label='Lin_I, '+bl)
    pl.plot(tim, snr0_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    ibl = ibl + 1

dmbd = np.array(dmbd, dtype=float)
dsnr = np.array(dsnr, dtype=float)
    
# fig.tight_layout()
fig1.tight_layout(rect=(0,0,1, 0.95))
fig2.tight_layout(rect=(0,0,1, 0.95))
fig3.tight_layout(rect=(0,0,1, 0.95))
fig4.tight_layout(rect=(0,0,1, 0.95))

fig1.text(0.02, 0.97, "Lin I MBD and Cir I MBD after PolConvert", \
             fontsize=11)

fig2.text(0.02, 0.97, "Mean Subtracted Lin I MBD and Cir I MBD after PConv",
          fontsize=11)

fig3.text(0.02, 0.97, "Lin I SNR and Cir I SNR after PolConvert", \
             fontsize=11)

fig4.text(0.02, 0.97, "Mean Subtracted Lin I SNR and Cir I SNR after PConv",
          fontsize=11)

fig5 = pl.figure()
fig6 = pl.figure()

pl.figure(fig5); pl.hist(dmbd*1e6, 51); pl.grid(1)
fig5.text(0.1, 0.9, "Differences of MBDs")
pl.figure(fig6); pl.hist(dsnr, 51); pl.grid(1)
fig6.text(0.1, 0.9, "Differences of SNRs")

# fig7 = pl.figure()
# fig8 = pl.figure()
# pl.figure(fig7); pl.hist(rmse_mbd, 14); pl.grid(1)
# fig7.text(0.1, 0.9, "Differences of MBD for 14 Baselines")
# pl.figure(fig8); pl.hist(rmse_snr, 14); pl.grid(1)
# fig8.text(0.1, 0.9, "Differences of SNR for 14 Baselines")




pl.show()





    
