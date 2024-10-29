help_text = '''
closure_delay.py - plot closure delay for  MBD or SBD.
    mbd:
    rmbd:
    tmbd:
    sbd:
    rsbd:
    tsbd:
'''

plotColorLegend =   False
plotAvailableTime = False #True


import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_time_var.py <par> [save], ")
    print("       where <par> is either MBD or SBD or SNR.")
    print("       save (optional): save  figures in pdf format.")
    sys.exit(0)

# arg_to_par = {'mbd':'mbdelay', 'sbd':'sbdelay', 'tmbd':'tot_mbd',
#               'tsbd':'tot_sbd', 'rmbd':'resid_mbd', 'rsbd':'resid_sbd'}
arg_to_par = {'mbd':'mbdelay', 'sbd':'sbdelay'}
    
pararg = (sys.argv[1]).lower()
if pararg not in arg_to_par.keys():
    print("Argument can be either %s or %s. Entered '%s'. "\
          "Exiting." % (*arg_to_par.keys(), sys.argv[1]))
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
from matplotlib.pyplot import cm
import matplotlib.patches as patches
from itertools import combinations
import copy

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

#
# Unpickle it:
#
with open('idx3819l.pkl', 'rb') as finp:
    idx3819l_1 = pickle.load(finp)

# with open('idx3819c.pkl', 'rb') as finp:
#     idx3819c_1 = pickle.load(finp)

with open('idx3819cI.pkl', 'rb') as finp:
    idx3819c_1 = pickle.load(finp)

#
# Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or 'snr' or
#     'tmbd' or 'tsbd' or 'rmbd' or 'rsbd'
#
parname = arg_to_par[pararg]

ps = "(ps)"

bls = list(idx3819l_1.keys())   # Baselines
bls.sort()                      # Lexigraphically sorted baselines
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
# Find all the baseline triplets with 3 stations (ie station triangles)
#
# trist: a string of three station letters, like 'EMS', 'MSY', 'TVY' etc.
#
#trian = {}
trians = [] # List of station riangles: 'EVY', 'EMV', 'ESV', 'ETV' etc
tribl = {}  # Dict to translate triangle to baselines, like MSV -> MS, SV, MV
ntri = 0   # Number of baseline triangles
for ab, bc, ca in combinations(bls, 3):
    stset = set(''.join(ab+bc+ca))
    trist = ''.join(sorted(stset))
    #print("trist = ", trist)
    #if len(stset) == 3:
    if len(trist) == 3:
        #print(stset)
        #print(trist)
        #print(ab, bc, ca)
        #trian[trist] = trist
        trians.append(trist)
        tribl[trist] = (ab, bc, ca)
        ntri = ntri + 1   # Number of baseline triangles

#
# Reorder tuples of baselines in tribl into pattern 'AB', 'BC', 'AC' 
#

tribl1 = {}
for trist in tribl.keys():
    l3bl = tribl[trist]  # List of the 3 baselines
    trian = l3bl
    
    st_end = l3bl[0][1]
    if l3bl[2][0] == st_end:
        trian = (l3bl[0], l3bl[2], l3bl[1])
        tribl1[trist] = trian
        print(tribl[trist], '->', trian)

    st_end = l3bl[1][1]
    if l3bl[0][0] == st_end:
        trian = (l3bl[1], l3bl[0], l3bl[2])
        tribl1[trist] = trian
        print(tribl[trist], '->', trian)
    elif l3bl[2][0] == st_end:
        trian = (l3bl[1], l3bl[2], l3bl[0])
        tribl1[trist] = trian
        print(tribl[trist], '->', trian)

tribl = tribl1
        
#sys.exit()


#
# To start processing from istart;  exclude bad data before istart.
#
istart = 2

tim = {}      # Time points. The gaps will be replaced with NaNs
tim1 = {}     # Original time points with some of them missing. 
# trul = 605.*np.arange(35) # Time ruler
par_l = {}
par_c = {}
snr_l = {}
snr_c = {}
snr_a = {}

#
# The arrays to contain the same data as the dictionaries par_l ans par_c etc.
#
atim   = np.zeros((ntri,35), dtype=float)
apar_l = np.zeros((ntri,35), dtype=float)
apar_c = np.zeros((ntri,35), dtype=float)
asnr_l = np.zeros((ntri,35), dtype=float)
asnr_c = np.zeros((ntri,35), dtype=float)

itri = 0
for bl in bls:
    tim1[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
    tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

    snr1_l = np.array(idx3819l_1[bl]['I']['snr'])[istart:]
    snr1_c = np.array(idx3819c_1[bl]['I']['snr'])[istart:]
    snr1_a = (abs(snr1_l.mean()) + abs(snr1_c.mean()))/2 # Avg Lin and Cir
    
    if parname == 'snr':
        par1_l = np.copy(snr1_l)
        par1_c = np.copy(snr1_c)
    else:
        par_l_us = np.array(idx3819l_1[bl]['I'][parname])[istart:] # In useconds
        par_c_us = np.array(idx3819c_1[bl]['I'][parname])[istart:] # In useconds
        par1_l = par_l_us*1e6           # Convert us to ps
        par1_c = par_c_us*1e6           # Convert us to ps
        
    #
    # Insert NaNs in the time gaps (ie over 605. seconds away)
    # Accordingly, insert NaNs in the parameter arrays
    #
    itim = np.int64(tim1[bl]/605) # Indices of non-NaN elements into tim and par
    tim[bl] = np.zeros(35)    
    tim[bl][:] = np.NaN
    tim[bl][itim] = tim1[bl]
    atim[itri,:] = tim[bl]
    
    snr_l[bl] = np.zeros(35)    
    snr_l[bl][:] = np.NaN
    snr_l[bl][itim] = snr1_l
    asnr_l[itri,:] = snr_l[bl]

    snr_c[bl] = np.zeros(35)    
    snr_c[bl][:] = np.NaN
    snr_c[bl][itim] = snr1_c
    asnr_c[itri,:] = snr_c[bl]

    par_l[bl] = np.zeros(35)    
    par_l[bl][:] = np.NaN
    par_l[bl][itim] = par1_l
    apar_l[itri,:] = par_l[bl]

    par_c[bl] = np.zeros(35)    
    par_c[bl][:] = np.NaN
    par_c[bl][itim] = par1_c
    apar_c[itri,:] = par_c[bl]

    itri = itri + 1

ntim = len(tim[bl]) # Any baseline bl, now tim sizes are the same for all

#
# Loop over the baseline triangles (bla, blb, blc) to find closure delays
#
tau_l = {} # Dict trist : array of closure delays (ab+bc+ca) for 35 times
tau_c = {} # Dict trist : array of closure delays (ab+bc+ca) for 35 times

#
# Arrays for tau
#
atau_l = np.zeros((ntri,35), dtype=float)
atau_c = np.zeros((ntri,35), dtype=float)

itri = 0
for trist in trians:
    print(trist, ': ', tribl[trist])
    ab, bc, ac = tribl[trist]
    tau_l[trist] = par_l[ab] + par_l[bc] - par_l[ac]
    tau_c[trist] = par_c[ab] + par_c[bc] - par_c[ac]

    atau_l[itri,:] = np.copy(tau_l[trist])
    atau_c[itri,:] = np.copy(tau_c[trist])
    
    itri = itri + 1

upar = parname.upper()

#cols = cm.rainbow(np.linspace(0, 1, ntri))
#cols = cm.gist_rainbow(np.linspace(0, 1, ntri))
#cols = cm.brg(1 - np.linspace(0, 1, ntri))
#cols = cm.jet(1 - np.linspace(0, 1, ntri))
cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))


# fig1, (ax11, ax12, ax13) = pl.subplots(3, 1, figsize=(8.4, 9),
#                               gridspec_kw={'height_ratios': [0.2, 0.6, 0.3]})
# GridSpecs
#gs = fig1.add_gridspec(ncols=2, nrows=3,
#                       width_ratios=[1, 1], height_ratios=[0.2, 0.6, 0.3])
#ax5 = pl.subplot2grid(shape=(3, 2), loc=(0, 0), colspan=3)
#fig1 = pl.figure(figsize=(8.4, 9))   #, constrained_layout=True)


ylms = -50000000   # +- ylimits, if ylms > 0


# fig1, axs = plt.subplots(ncols=2, nrows=3, figsize=(8.4, 9),
#                          constrained_layout=True,
#                          gridspec_kw={'height_ratios' : [0.2, 0.6, 0.3]})

# for ax in axs[0,:]:   # Remove all the 2 top axes
#     ax.remove()

# gs = axs[0, 0].get_gridspec()  # Get GridSpec from ANY of the axes

# ax0 = fig1.add_subplot(gs[0,:])  # Create one top axis

gs_kw1 = dict(width_ratios=[1, 1], height_ratios=[0.15, 0.6, 0.25])
fig1, axd1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                               ['distr_frfit', 'distr_pconv'],
                               ['hist_frfit', 'hist_pconv']],
                              gridspec_kw=gs_kw1, figsize=(8.4, 8),
                              layout="constrained")
#
# Plot color legend on top
#
ax_col = axd1['col_legend']

qntri = ntri//4         # Quarter of the number of triangles
ax_col.set_xlim(-0.1, 7.5)
ax_col.set_ylim(-0.2, qntri+0.5)

j = 0 
for i in range(qntri):
    y = qntri - i - 1      # Vertical patch position
    k = i                  # Index into trians[k] and cols[k,:]
    rect = patches.Rectangle((0, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(1.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

    k = i + qntri
    rect = patches.Rectangle((2, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(3.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

    k = i + 2*qntri
    rect = patches.Rectangle((4, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(5.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

    k = i + 3*qntri
    rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(7.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

ax_col.set_axis_off()
ax_col.set_title("%s Closure" % upar, fontsize=16)

timx = tim['VY']/3600

ax_ffd = axd1['distr_frfit']  # Plot distr of FourFit pseudo-I param vs Time

for ic in range(ntri):
    trist = trians[ic]
    ax_ffd.plot(timx, tau_c[trist], '.', color=cols[ic,:])
ax_ffd.grid(1)
ax_ffd.set_title("Fourfit Pseudo-I, %s vs Time" % upar)
ax_ffd.set_xlabel("hours", fontsize=14)
ax_ffd.set_ylabel("ps", fontsize=14)
ax_ffd.yaxis.set_label_coords(-0.05, 0.58)
ax_ffd.xaxis.set_label_coords(0.55, 0.07)
if ylms > 0: ax_ffd.set_ylim(-ylms, ylms)
if pararg == 'mbd':
    ax_ffd.set_ylim(-530, 620)
elif pararg == 'sbd':
    ax_ffd.set_ylim(-1990, 870)

    
ax_ffh = axd1['hist_frfit']  # Plot hist of FourFit I param

ax_ffh.hist(abs(atau_c.flatten()), 50)
ax_ffh.grid(1)
ax_ffh.set_xlabel("ps", fontsize=14)
ax_ffh.set_title("Fourfit Pseudo-I, abs(%s)" % upar)
ax_ffh.xaxis.set_label_coords(0.5, -0.12)


ax_pcd = axd1['distr_pconv'] # Plot distr of PolConvert  pseudo-I param vs Time

for ic in range(ntri):
    trist = trians[ic]
    ax_pcd.plot(timx, tau_l[trist], '.', color=cols[ic,:])
ax_pcd.grid(1)
ax_pcd.set_title("PolConvert I, %s vs Time" % upar)
ax_pcd.set_xlabel("hours", fontsize=14)
ax_pcd.set_ylabel("ps", fontsize=14)
ax_pcd.yaxis.set_label_coords(-0.05, 0.58)
ax_pcd.xaxis.set_label_coords(0.55, 0.07)
if pararg == 'mbd':
    ax_pcd.set_ylim(-530, 620)
elif pararg == 'sbd':
    ax_pcd.set_ylim(-1990, 870)


ax_pch = axd1['hist_pconv']  # Plot hist of PolConvert I param

ax_pch.hist(abs(atau_l.flatten()), 50)
ax_pch.grid(1)
ax_pch.set_xlabel("ps", fontsize=14)
ax_pch.set_title("PolConvert I, abs(%s)" % upar)
ax_pch.xaxis.set_label_coords(0.5, -0.12)


pl.savefig("%s_Closure_Delay.pdf" % upar, format='pdf')


#
# Plot available times for each baseline
#
if plotAvailableTime:
    #cols_bl = cm.rainbow(np.linspace(0, 1, nbls))
    #cols_bl = cm.nipy_spectral(np.linspace(0, 1, nbls))
    #cols_bl = cm.gist_rainbow(np.linspace(0, 1, nbls))
    cols_bl = cm.jet(np.linspace(0, 1, nbls))
    
    fig4, ax41 =  pl.subplots()
    
    fig4.text(0.22, 0.95, "Baseline Times with Missed Scans", fontsize=14)

    #sh = 1 # Just arbitrary shift to plot the lines 
    sh = np.ones(ntim) # Horizontal line with gaps
    for ib in range(nbls):
        bl = bls[ib]
        t = tim[bl]/3600
        y = nbls - ib
        yy = y*sh                    # Array of heights
        ax41.plot(t, yy, color=cols_bl[ib,:], lw=3)
        ax41.plot(t, yy, 'k.', markersize=5)
        ax41.text(-0.55, y-0.35, bl, fontsize=14)
        print("ib = %d, y = %d" % (ib, y))

    ax41.grid(True)
    ax41.set_xlabel("hours", fontsize=14)
    ax41.set_yticks([])
    ax41.set_ylim(0, nbls+1)
    ax41.set_xlim(-0.8, 6)
    fig4.tight_layout(rect=(0.00, 0.00, 0.98, 0.95))

    pl.savefig("Gaps_in_Time.pdf", format='pdf')

    

#
# Plot triangle color legend
#
if plotColorLegend:
    fig3, ax6 = pl.subplots()

    qntri = ntri//4         # Quarter of the number of triangles
    j = 0 
    for i in range(qntri):
        y = qntri - i - 1      # Vertical patch position
        k = i                  # Index into trians[k] and cols[k,:]
        rect = patches.Rectangle((0, y), .95, .95, facecolor=cols[k,:])
        ax6.add_patch(rect)
        ax6.text(1.1, y+0.3, trians[k], fontsize=14)
        print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + qntri
        rect = patches.Rectangle((2, y), .95, .95, facecolor=cols[k,:])
        ax6.add_patch(rect)
        ax6.text(3.1, y+0.3, trians[k], fontsize=14)
        print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 2*qntri
        rect = patches.Rectangle((4, y), .95, .95, facecolor=cols[k,:])
        ax6.add_patch(rect)
        ax6.text(5.1, y+0.3, trians[k], fontsize=14)
        print("i = %d, y = %d, k = %2d" % (i, y, k))
        
        k = i + 3*qntri
        rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
        ax6.add_patch(rect)
        ax6.text(7.1, y+0.3, trians[k], fontsize=14)
        print("i = %d, y = %d, k = %2d" % (i, y, k))
        

    pl.ylim(0, qntri)
    pl.xlim(0, 7.5)

    ax6.set_axis_off()

    pl.savefig("Triangle_color_legend.pdf", format='pdf')




#============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#-------------------   With Y and without Y  --------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#============================================================================

#
# Create triangle lists with and without station 'Y'
#
trians_with_y = []
trians_sans_y = []

for ic in range(ntri):
    trist = trians[ic]
    print('trist = ', trist)
    if 'Y' in trist:
        trians_sans_y.append(trist)
    else:
        trians_with_y.append(trist)

print("trians_sans 'Y': ", trians_sans_y)
print("trians_with 'Y':    ", trians_with_y)






!!!!!!!!!!!! ????????????????? 







#
# Plot closures with and without station 'Y'
#
gs_kw1 = dict(width_ratios=[1, 1, 1, 1, 1],
              height_ratios=[0.10, 0.3, 0.15, 0.3, 0.15])
fig1, axd1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                               ['distr_frfit_sansY', 'distr_pconv_sansY'],
                               ['hist_frfit_sansY', 'hist_pconv_sansY'],
                               ['distr_frfit_withY', 'distr_pconv_withY'],
                               ['hist_frfit_withY', 'hist_pconv_withY']],
                              gridspec_kw=gs_kw1, figsize=(8.4, 10),
                              layout="constrained")
#
# Plot color legend on top
#
ax_col = axd1['col_legend']

qntri = ntri//4         # Quarter of the number of triangles
ax_col.set_xlim(-0.1, 7.5)
ax_col.set_ylim(-0.2, qntri+0.5)

j = 0 
for i in range(qntri):
    y = qntri - i - 1      # Vertical patch position
    k = i                  # Index into trians[k] and cols[k,:]
    rect = patches.Rectangle((0, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(1.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

    k = i + qntri
    rect = patches.Rectangle((2, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(3.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

    k = i + 2*qntri
    rect = patches.Rectangle((4, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(5.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

    k = i + 3*qntri
    rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(7.1, y+0.2, trians[k], fontsize=12)
    print("i = %d, y = %d, k = %2d" % (i, y, k))

ax_col.set_axis_off()
ax_col.set_title("%s Closure" % upar, fontsize=16)

timx = tim['VY']/3600

ax_ffd = axd1['distr_frfit_sansY']  # Plot distr of FourFit pseudo-I param no Y

for ic in range(ntri):
    trist = trians[ic]
    ax_ffd.plot(timx, tau_c[trist], '.', color=cols[ic,:])
ax_ffd.grid(1)
ax_ffd.set_title("Fourfit Pseudo-I, %s vs Time" % upar)
ax_ffd.set_xlabel("hours", fontsize=14)
ax_ffd.set_ylabel("ps", fontsize=14)
ax_ffd.yaxis.set_label_coords(-0.05, 0.58)
ax_ffd.xaxis.set_label_coords(0.55, 0.07)
if ylms > 0: ax_ffd.set_ylim(-ylms, ylms)
if pararg == 'mbd':
    ax_ffd.set_ylim(-530, 620)
elif pararg == 'sbd':
    ax_ffd.set_ylim(-1990, 870)

    
ax_ffhy = axd1['hist_frfit_sansY']  # Plot hist of FourFit I param

ax_ffh.hist(abs(atau_c.flatten()), 50)
ax_ffh.grid(1)
ax_ffh.set_xlabel("ps", fontsize=14)
ax_ffh.set_title("Fourfit Pseudo-I, abs(%s)" % upar)
ax_ffh.xaxis.set_label_coords(0.5, -0.12)




ax_ffd = axd1['distr_frfit_withY']  # Plot distr of FourFit pseudo-I param w/Y

for ic in range(ntri):
    trist = trians[ic]
    ax_ffd.plot(timx, tau_c[trist], '.', color=cols[ic,:])
ax_ffd.grid(1)
ax_ffd.set_title("Fourfit Pseudo-I, %s vs Time" % upar)
ax_ffd.set_xlabel("hours", fontsize=14)
ax_ffd.set_ylabel("ps", fontsize=14)
ax_ffd.yaxis.set_label_coords(-0.05, 0.58)
ax_ffd.xaxis.set_label_coords(0.55, 0.07)
if ylms > 0: ax_ffd.set_ylim(-ylms, ylms)
if pararg == 'mbd':
    ax_ffd.set_ylim(-530, 620)
elif pararg == 'sbd':
    ax_ffd.set_ylim(-1990, 870)

    
ax_ffh = axd1['hist_frfit_withY']  # Plot hist of FourFit I param

ax_ffh.hist(abs(atau_c.flatten()), 50)
ax_ffh.grid(1)
ax_ffh.set_xlabel("ps", fontsize=14)
ax_ffh.set_title("Fourfit Pseudo-I, abs(%s)" % upar)
ax_ffh.xaxis.set_label_coords(0.5, -0.12)






ax_pcd = axd1['distr_pconv_sansY'] # Plot distr of PolConv ps-I param no Y

for ic in range(ntri):
    trist = trians[ic]
    ax_pcd.plot(timx, tau_l[trist], '.', color=cols[ic,:])
ax_pcd.grid(1)
ax_pcd.set_title("PolConvert I, %s vs Time" % upar)
ax_pcd.set_xlabel("hours", fontsize=14)
ax_pcd.set_ylabel("ps", fontsize=14)
ax_pcd.yaxis.set_label_coords(-0.05, 0.58)
ax_pcd.xaxis.set_label_coords(0.55, 0.07)
if pararg == 'mbd':
    ax_pcd.set_ylim(-530, 620)
elif pararg == 'sbd':
    ax_pcd.set_ylim(-1990, 870)


ax_pch = axd1['hist_pconv_sansY']  # Plot hist of PolConvert I param

ax_pch.hist(abs(atau_l.flatten()), 50)
ax_pch.grid(1)
ax_pch.set_xlabel("ps", fontsize=14)
ax_pch.set_title("PolConvert I, abs(%s)" % upar)
ax_pch.xaxis.set_label_coords(0.5, -0.12)


pl.savefig("%s_Closure_Delay_y_no_y.pdf" % upar, format='pdf')





    
pl.show()




