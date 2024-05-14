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

fig1 = pl.figure(figsize=(8, 12))
fig2 = pl.figure(figsize=(8, 12))
fig3 = pl.figure(figsize=(8, 12))
fig4 = pl.figure(figsize=(8, 12))

#
# To start plotting from istart;  exclude bad data before istart.
#
istart = 2

ibl = 1   # Baseline
for bl in bls:   # Loop over the baselines
    tim = np.array(idx3819l_1[bl]['I']['time'])[istart:] / 60
    tim = tim - tim[0]
    mbd_l = np.array(idx3819l_1[bl]['I']['mbdelay'])[istart:]
    mbd_c = np.array(idx3819c_1[bl]['I']['mbdelay'])[istart:]
    snr_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]

    #
    # Plot MBD
    #
    pl.figure(fig1)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, mbd_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    #
    # Subtract MBD means
    #
    mbd0_l = mbd_l - mbd_l.mean()
    mbd0_c = mbd_c - mbd_c.mean()
    
    pl.figure(fig2)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, mbd0_l , label='Lin_I, '+bl)
    pl.plot(tim, mbd0_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    #
    # Plot SNR
    #
    pl.figure(fig3)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, snr_l , label='Lin_I, '+bl)
    pl.plot(tim, snr_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    #
    # Subtract SNR means
    #
    snr0_l = snr_l - snr_l.mean()
    snr0_c = snr_c - snr_c.mean()
    
    pl.figure(fig4)
    pl.subplot(5, 3, ibl)
    pl.plot(tim, snr0_l , label='Lin_I, '+bl)
    pl.plot(tim, snr0_c , label='Cir_I, '+bl)
    pl.grid(True)
    pl.legend(loc='upper right')

    ibl = ibl + 1
    
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



pl.show()





    
