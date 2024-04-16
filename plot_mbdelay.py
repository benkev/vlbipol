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
    idx3819cI_1 = pickle.load(finp)


bls = list(idx3819l_1.keys())   # Baselines
bls.sort()                      # Lexigraphically sorted baselines

# fig, ax = pl.subplots(figsize=(8, 12))

fig = pl.figure(figsize=(8, 12))

ibl = 1
for bl in bls:
    t_l =   np.array(idx3819l_1[bl]['I']['time']) / 60
    t_cI = np.array(idx3819cI_1[bl]['I']['time']) / 60
    t_l = t_l - t_l[0]
    t_cI = t_cI - t_cI[0]
    mbd_l =   np.array(idx3819l_1[bl]['I']['mbdelay'])
    mbd_cI = np.array(idx3819cI_1[bl]['I']['mbdelay'])
   
    pl.subplot(5, 3, ibl)
    pl.plot(t_l,  mbd_l ,  label='Lin_I, '+bl)
    pl.plot(t_cI, mbd_cI , label='Cir_I, '+bl)
    pl.legend(loc='upper right')
    ibl = ibl + 1
    
# fig.tight_layout()
fig.tight_layout(rect=(0,0,1, 0.95))

pl.figtext(0.02, 0.97, "3819 Pseudo-Stokes I MultiBand Delays vs Time (min), " \
           "Linear & Circular Polarization after PolConvert", fontsize=11)
pl.show()





    
