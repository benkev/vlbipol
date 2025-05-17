import matplotlib.pyplot as pl
import numpy as np
import pickle

with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)


sr = '1803+784'

thr_ghs = closl[sr]['GHS']['thour']
thr_ght = closl[sr]['GHT']['thour']

clp_ghs = closl[sr]['GHS']['cloph']
clp_ght = closl[sr]['GHT']['cloph']

clm_ghs = closl[sr]['GHS']['tau_mbd']
clm_ght = closl[sr]['GHT']['tau_mbd']


f1 = pl.figure(figsize=(6.4, 7))

ax1 = pl.subplot(2,1,1)
ax1.plot(thr_ghs, clp_ghs, 'r.', ms=8, label='GHS')
ax1.plot(thr_ght, clp_ght, 'g.', ms=8, label='GHT')
ax1.grid(1)
ax1.set_yticks(-120 + 30*np.arange(9))
ax1.set_title("VO2187, Source 1803+784: Phase Closures")
ax1.set_xlabel("t (hours)", labelpad=-35)
ax1.set_ylabel("deg")
ax1.legend()

ax2 = pl.subplot(2,1,2)
ax2.plot(thr_ghs, clm_ghs, 'r.', ms=8, label='GHS')
ax2.plot(thr_ght, clm_ght, 'g.', ms=8, label='GHT')
ax2.grid(1)
ax2.set_ylim(-55, 55)
ax2.set_title("VO2187, Source 1803+784: MBD Closures")
ax2.set_xlabel("t (hours)", labelpad=-35)
ax2.set_ylabel("ps")
ax2.legend()


