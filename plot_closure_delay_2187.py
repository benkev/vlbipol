help_text = '''
plot_closure_delay_2187.py - plot closure delay for residual and
                             total MBD or SBD, session 2187 (VO2187).
Usage:
%run plot_closure_delay_2187.py <par> [save]

where <par> is one of the following parameters:
    mbd:  
    sbd:
    tmbd:
    tsbd:

    save is optional: save figures in pdf format.
'''

plotColorLegend =   False
plotAvailableTime = False #True


import sys

if len(sys.argv) < 2  or sys.argv[1] == '--help':
    print(help_text)
    print("Usage:")
    print("python plot_closure_delay.py <par> [save], ")
    print("       where <par> is either MBD or SBD or TMBD or TSBD.")
    print("       save (optional): save  figures in pdf format.")
    sys.exit(0)

arg_to_par = {'mbd':'mbdelay', 'sbd':'sbdelay', 'tmbd':'tot_mbd',
              'tsbd':'tot_sbd'}
    
pararg = (sys.argv[1]).lower()
if pararg not in arg_to_par.keys():
    print("Argument can be either %s or %s or %s or %s. Entered '%s'. "\
          "Exiting." % (*arg_to_par.keys(), sys.argv[1]))
    sys.exit(0)

sf = False  # Save figure request
if len(sys.argv) == 3:
    if sys.argv[2] == 'save':
        sf = True
    else:
        print("Argument can only be 'save'. Entered '%s'. Exiting." %
              sys.argv[2])
        sys.exit(0)
    
    
import pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib.patches as patches
from itertools import combinations
import copy


def find_baseline_triangles(bls):
    '''
    Find all the baseline triplets with 3 stations (ie station triangles)
    such that their order makes correct closure, like 'MT', 'TE', 'ME'.
    The last baseline is always in inverted order (here not 'EM' but 'ME'),
    so in the closure delay it must be present with the minus sign. For example,

        tau_MTE = tau_MT + tau_TE - tau_ME
    
    Inpit:
        bls: list of baseline names (of a particular session).
             Example: ['EV', 'EY', 'ME', 'MS', 'MT', 'MV', 'MY', 'SE', 'SV', \
                       'SY', 'TE', 'TV', 'TY', 'VY']
    Returns:
        
        tribl: dictionary of baseline triplets, the keys being triangle names
               as 3-letter strings, like 'MTE', and values being the 3-element
               lists of the baselines comprising the closure triangle.
    Example:
        tribl['MTE'] -> ['MT', 'TE', 'ME'].
    '''

    bls.sort()                       # Lexigraphically sort baselines

    # Set of station letters, stset
    ststr = ''
    for bl in bls: ststr = ststr + bl  # Concatenate baseline strings in ststr
    stset = set(ststr)  # Leave only unique station letters in the sts set

    # String of station letters ststr
    nsts = len(stset)
    ststr = ''.join(sorted(stset))

    #
    # trist: a string of three station letters, like 'EMS', 'MSY', 'TVY' etc.
    #
    #trian = {}
    trians = [] # List of station riangles: 'EVY', 'EMV', 'ESV', 'ETV' etc
    tribl = {}  # Dict triangle -> baselines, like MSV -> MS, SV, MV
    ntri = 0    # Number of baseline triangles
    for ab, bc, ca in combinations(bls, 3):
        stset = set(''.join(ab+bc+ca))
        trist = ''.join(sorted(stset))
        if len(trist) == 3:
            trians.append(trist)
            tribl[trist] = (ab, bc, ca)
            ntri = ntri + 1   # Number of baseline triangles

    #
    # Reorder tuples of baselines in tribl into pattern 'AB', 'BC', 'AC' 
    #
    print("The baseline triplets are reordered:")

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

    return tribl




def plot_closure_legend(ax_col, trians, cols, par, fs=12):
    '''
    Plot color legend for the closure triangles in a separate axis ax_col.
    For 16 closure triangles, it plots 4 colimns by 4 in height, each showing
    the color and the triangle name on the right.
    Inputs:
        ax_col: axis; better to have shorter in height and wider.
                It is assumed to be on top of other plotsand its title
                describes the whole figure.
        trians: list of 3-letter closure triangle names.
        cols:   Array of colors, ntri by 4, from pyplot.cm(numbers 0 to 1)
                For example:
                    cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))
        par:    The parameter name to print in the title in uppercase.
        fz:     Fontsise of the printed triangle names
    '''
    ntri = len(trians)
    qntri = ntri//4         # Quarter of the number of triangles
    #cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))           # Colors
    upar = par.upper() # The uppercase parameter name to print in the title.
    
    ax_col.set_xlim(-0.1, 7.5)
    ax_col.set_ylim(-0.2, qntri+0.5)

    j = 0 
    for i in range(qntri):
        y = qntri - i - 1      # Vertical patch position
        k = i                  # Index into trians[k] and cols[k,:]
        rect = patches.Rectangle((0, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(1.1, y+0.2, trians[k], fontsize=fs)
        #print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + qntri
        rect = patches.Rectangle((2, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(3.1, y+0.2, trians[k], fontsize=fs)
        #print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 2*qntri
        rect = patches.Rectangle((4, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(5.1, y+0.2, trians[k], fontsize=fs)
        #print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 3*qntri
        rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(7.1, y+0.2, trians[k], fontsize=fs)
        #print("i = %d, y = %d, k = %2d" % (i, y, k))

    ax_col.set_axis_off()
    ax_col.set_title("%s Closure" % upar, fontsize=16)




def plot_closures_distr(ax_distr, timh, atau, seltri, cols, ylim, pararg, ttl):
    '''
    Plot distribution of a delay closures.
    '''

    if isinstance(seltri, (list, tuple, np.ndarray)):
        sel = np.array(seltri, dtype='bool')  # Any sequence into bool array sel
        for ic in range(ntri):
            if sel[ic]:
                ax_distr.plot(timh, atau[ic,:], '.', color=cols[ic,:])
    else: # if seltri is elemental, e.g. any number, plot all
        for ic in range(ntri):
            ax_distr.plot(timh, atau[ic,:], '.', color=cols[ic,:])
        
    ax_distr.grid(1)
    ax_distr.set_title(ttl)
    ax_distr.set_xlabel("hours", fontsize=14)
    ax_distr.set_ylabel("ps", fontsize=14)
    ax_distr.yaxis.set_label_coords(-0.05, 0.58)
    ax_distr.xaxis.set_label_coords(0.55, 0.07)

    ax_distr.set_ylim(ylim)
    


    

# def plot_closures_hist(ax_hist, timh, atau, seltri, pararg, ttl):
#     '''
#     Plot distribution and histogram of a delay closures.
#     '''

#     if isinstance(seltri, (list, tuple, np.ndarray)):
#         sel = np.array(seltri, dtype='bool')  # Any seq. into bool array sel
#         # print("atau[sel,:].flatten().shape = ", atau[sel,:].flatten().shape)
#         ax_hist.hist(abs(atau[sel,:].flatten()), 50)
#     else: # if seltri is elemental, e.g. any number:
#         ax_hist.hist(abs(atau.flatten()), 50)

#     ax_hist.grid(1)
#     ax_hist.set_xlabel("ps", fontsize=14)
#     ax_hist.set_title(ttl)
#     ax_hist.xaxis.set_label_coords(0.5, -0.12)




    
def plot_closures_hist_horiz(ax_hist, timh, atau, seltri, pararg, xlims, ylims,
                             colr, ttl):
    '''
    Plot distribution and histogram of a delay closures.

    ax_hist: axis
    atau[(ntri,ntim]: closure delay array
    seltri[ntri]: bool array with True at the selected triangles to plot.
                  If seltri is not a sequence (say, any number), all tau
                  are plotted.
    
    '''
 
    # print("atau.flatten()" % )
    
    if isinstance(seltri, (list, tuple, np.ndarray)):
        sel = np.array(seltri, dtype='bool')  # Any sequence into bool array sel
        # print("atau[sel,:].flatten().shape = ", atau[sel,:].flatten().shape)
        ax_hist.hist(atau[sel,:].flatten(), 50, color=colr,
                     orientation='horizontal')
    else: # if seltri is elemental, e.g. any number:
        ax_hist.hist(atau.flatten(), 50, color=colr,
                     orientation='horizontal')
    ax_hist.set_xlim(xlims)
    ax_hist.set_ylim(ylims)
    ax_hist.grid(1)
    #ax_hist.set_xlabel("ps", fontsize=14)
    ax_hist.set_title(ttl)
    #ax_hist.xaxis.set_label_coords(0.5, -0.12)
    #ax_hist.set_xticks([]);
    #                ax_hist.tick_params(axis='x', rotation=90)
    #xl = ax_hist.get_xticklabels()
    #ax_hist.set_xticklabels(xl, rotation=-90); 
    #ax_hist.set_yticks([]); 
    ax_hist.set_yticklabels('')


    
# def plot_closures_distr_and_hist_sel(ax_distr, ax_hist, timh, tau, tau_sel,
#                                 trians, trisel, cols, pararg, ttl):
#     '''
#     Plot distribution and histogram of a delay closures for only selected
#     triangles listed in trisel.
#     '''
#     for ic in range(ntri):
#         trist = trians[ic]
#         if trist in trisel:
#             ax_distr.plot(timh, tau_l[trist], '.', color=cols[ic,:])
#     ax_distr.grid(1)
#     ax_distr.set_title(ttl)
#     ax_distr.set_xlabel("hours", fontsize=14)
#     ax_distr.set_ylabel("ps", fontsize=14)
#     ax_distr.yaxis.set_label_coords(-0.05, 0.58)
#     ax_distr.xaxis.set_label_coords(0.55, 0.07)
#     if ylms > 0: ax_distr.set_ylim(-ylms, ylms)
#     if pararg == 'mbd':
#         ax_distr.set_ylim(-530, 620)
#     elif pararg == 'sbd':
#         ax_distr.set_ylim(-1990, 870)


#     ax_hist.hist(abs(tau_sel), 50)
#     ax_hist.grid(1)
#     ax_hist.set_xlabel("ps", fontsize=14)
#     ax_hist.set_title("Fourfit Pseudo-I, abs(%s) sans Y" % upar)
#     ax_hist.xaxis.set_label_coords(0.5, -0.12)

































    
# ========================== Body of the Code ================================

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

#
# Unpickle it:
#
with open('idx2187lI.pkl', 'rb') as finp:
    idxl = pickle.load(finp)

with open('idx2187cI.pkl', 'rb') as finp:
    idxc = pickle.load(finp)

#
# Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or
#     'tmbd' or 'tsbd' or 'rmbd' or 'rsbd'
#
parname = arg_to_par[pararg]

ps = "(ps)"

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

#
# Find tribl, the dictionary with all the baseline triplets with 3 stations
# (ie station triangles) such that their order makes a correct closure,
# like 'MT', 'TE', 'ME'.
# The last baseline is always in inverted order (here not 'EM' but 'ME'),
# so in the closure delay it must be present with the minus sign. For example,
#
#        tau_MTE = tau_MT + tau_TE - tau_ME
#

tribl = find_baseline_triangles(bls)   # ============== CALL ============== >>

trians = list(tribl.keys())
ntri = len(tribl)

print()
print('bls = ', bls, '\n')
print('tribl:')
for trist in tribl.keys():
    print(trist, " -> ", tribl[trist])
    
#sys.exit(0)

tim1 = {}     # Original time points with some of them missing. 
tim = {}      # Time points. The gaps will be replaced with NaNs

for bl in bls:
    tim1[bl] = np.array(idxl[bl]['I']['time'])  #/ 3600 # Sec -> hours
    tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

#
# Find the minimum time between the scans over all the baselines
#
min_t_scan = 1000000000

for bl in bls:
    t_scans = np.diff(tim1[bl])
    bl_min_t_scan = np.min(t_scans)
    if bl_min_t_scan < min_t_scan:
        min_t_scan = bl_min_t_scan

#
# Search for the maximum number of time counts, ntim,
# among all the baselines. At the end of the loop below, ntim will hold
# the length for the parameter storage in arrays.
# The parameter sets shorter than ntim contain time gaps to be filled with NaNs.
#
ntim = -1       # Will contain the maximum time count
for bl in bls:
    bl_ntim = np.max(tim1[bl]/min_t_scan)    # Max t counts for the baseline
    #print("bl_ntim = %f" % bl_ntim)
    if bl_ntim > ntim:
        ntim = np.int64(np.ceil(bl_ntim))

    print("len(tim1['%s']) = %d; Max t counts = %f" %
                                   (bl, len(tim1[bl]), bl_ntim))

# ?? ntim = ntim + 1 # !!!!!!! I DO NOT KNOW WHY 3819 NEEDS IT ???????????????
    
print("Max time counts: %d;  min scan time: %d s." % (ntim, min_t_scan))


#
# Create array of all the sources
#
srcl = []
srcc = []
for bl in bls:
    srcl.extend(idxl[bl]['I']['source'])
    srcc.extend(idxl[bl]['I']['source'])
lsrc = srcl + srcc
asrc = np.unique(lsrc)     # Turn into np.array leaving only unique source names
asrc.sort()                # Sort the source names lexicographically

nsrc = len(asrc)    # Number of sources

#
# Create hash-table (dictionary) of <source> : <its index into asrc>
#
idxs = {}
for i in range(nsrc):
    s = asrc[i]
    idxs[s] = i


#
# FOR A CERTAIN BASELINE bl !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
bl = 'GE'
#
# In dict src_tim each source points at the list os scan times
#
src_tim = {}
for s in asrc:         # Give empty list values to each source key 
    src_tim[s] = []

bl_ntim = len(idxl[bl]['I']['time']) # Number of time counts for baseline bl

for i in range(bl_ntim):  # Gather the scan times in lists for each source
    src_tim[idxl[bl]['I']['source'][i]].append(tim1[bl][i])











par_l = {}
par_c = {}
snr_l = {}
snr_c = {}
snr_a = {}

#
# The arrays to contain the same data as the dictionaries par_l ans par_c etc.
#
atim   = np.zeros((ntri,ntim), dtype=float)
apar_l = np.zeros((ntri,ntim), dtype=float)
apar_c = np.zeros((ntri,ntim), dtype=float)
asnr_l = np.zeros((ntri,ntim), dtype=float)
asnr_c = np.zeros((ntri,ntim), dtype=float)

itri = 0       # Array index of a baseline triangle
for bl in bls:
    snr1_l = np.array(idxl[bl]['I']['snr'])
    snr1_c = np.array(idxc[bl]['I']['snr'])
    snr1_a = (abs(snr1_l.mean()) + abs(snr1_c.mean()))/2 # Avg Lin and Cir
    
    if parname == 'snr':
        par1_l = np.copy(snr1_l)
        par1_c = np.copy(snr1_c)
    else:
        par1_l = np.array(idxl[bl]['I'][parname])*1e6    # In picoseconds
        par1_c = np.array(idxc[bl]['I'][parname])*1e6    # In picoseconds
        
    #
    # Insert NaNs in the time gaps (ie over 605. seconds without a time count)
    # Accordingly, insert NaNs in the parameter arrays
    #
    itim = np.int64(tim1[bl]/min_t_scan) # Indices of non-NaN elements
    tim[bl] = np.zeros(ntim)    
    tim[bl][:] = np.nan
    tim[bl][itim] = tim1[bl]
    atim[itri,:] = tim[bl]
    
    snr_l[bl] = np.zeros(ntim)    
    snr_l[bl][:] = np.nan
    snr_l[bl][itim] = snr1_l
    asnr_l[itri,:] = snr_l[bl]

    snr_c[bl] = np.zeros(ntim)    
    snr_c[bl][:] = np.nan
    snr_c[bl][itim] = snr1_c
    asnr_c[itri,:] = snr_c[bl]

    par_l[bl] = np.zeros(ntim)    
    par_l[bl][:] = np.nan
    par_l[bl][itim] = par1_l
    apar_l[itri,:] = par_l[bl]

    par_c[bl] = np.zeros(ntim)    
    par_c[bl][:] = np.nan
    par_c[bl][itim] = par1_c
    apar_c[itri,:] = par_c[bl]

    itri = itri + 1


#
# Loop over the baseline triangles (bla, blb, blc) to find closure delays
#

tau_l = {} # Dict trist : array of closure delays (ab+bc+ca) for ntim times
tau_c = {} # Dict trist : array of closure delays (ab+bc+ca) for ntim times

#
# Arrays for tau
#
atau_l = np.zeros((ntri,ntim), dtype=float)
atau_c = np.zeros((ntri,ntim), dtype=float)

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

#
# Bool arrays to mask trians with and without station Y
#

sel_Y = np.zeros(ntri, dtype='bool')    # With Y
sel_noY = np.zeros(ntri, dtype='bool')  # Without Y

for itri in range(ntri):
    if 'Y' in trians[itri]:
        sel_Y[itri] = True
        print("Y in trians[%d] = %s" % (itri, trians[itri]))
    else:
        sel_noY[itri] = True
        print("Y not in trians[%d] = %s" % (itri, trians[itri]))



#
# ========================  PLOTTING  =============================
# 


#cols = cm.rainbow(np.linspace(0, 1, ntri))
#cols = cm.gist_rainbow(np.linspace(0, 1, ntri))
#cols = cm.brg(1 - np.linspace(0, 1, ntri))
#cols = cm.jet(1 - np.linspace(0, 1, ntri))

cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))


#
# Do not plot total_mbd and total_sbd closures with Y station 
#
if pararg == 'mbd' or pararg == 'sbd':

    # gs_kw1 = dict(width_ratios=[0.3, 0.7], height_ratios=[0.15, 0.425, 0.425])
    # fig1, axd1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
    #                                ['hist_frfit', 'distr_frfit'],
    #                                ['hist_pconv', 'distr_pconv']],
    #                                gridspec_kw=gs_kw1, figsize=(8.4, 8),
    #                                layout="constrained")

    gs_kw1 = dict(width_ratios=[0.75, 0.25], height_ratios=[0.15, 0.425, 0.425])
    fig1, axd1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                                   ['distr_frfit', 'hist_frfit'],
                                   ['distr_pconv', 'hist_pconv']],
                                   gridspec_kw=gs_kw1, figsize=(8.4, 8),
                                   layout="constrained")
    #
    # Plot color legend on top
    #
    ax_col = axd1['col_legend']    # Get the axis for color legend

    plot_closure_legend(ax_col, trians, cols, upar)  # =======  CALL ======= >>



    if pararg == 'mbd':
        # ylim_distr = (-530, 620)  # It includes the outliers beyond +-400 ps
        ylim_distr = (-220, 220)
        xlim_hist = (-5, 170)  # Now it is 'xlim' after 90-deg rotation
    elif pararg == 'sbd':
        ylim_distr = (-1200, 1000)
        xlim_hist = (-5, 100)
    # elif pararg == 'tmbd':
    #     ylim_distr = (-220, 220)
    #     xlim_hist = (-5, 170)  # Now it is 'xlim' after 90-deg rotation
    # elif pararg == 'tsbd':
    #     ylim_distr = (-1200, 1000)
    #     xlim_hist = (-5, 100)



    timh = tim['TV']/3600   # Time counts (in hours) from baseline with no NaNs

    nfinite = np.count_nonzero(np.isfinite(atau_l)) # Non-NaNs in atau_l & _c

    hist_colr = 'red'

    #fig1b = pl.figure(); ax_ffd = pl.gca()
    ax_ffd = axd1['distr_frfit'] # Plot distr of FourFit ps.-I param vs Time  
    ttl_ffd = "Fourfit Pseudo-I, %s vs Time (%d triangles)" % (upar, ntri)
    plot_closures_distr(ax_ffd, timh, atau_l,  # ========== CALL ============ >>
                       1, cols, ylim_distr, pararg, ttl_ffd)

    #fig1a = pl.figure(); ax_ffh = pl.gca()
    ax_ffh = axd1['hist_frfit']  # Plot hist of FourFit I param
    ttl_ffh = "%s (%d points)" % (upar, nfinite)
    plot_closures_hist_horiz(ax_ffh, timh, atau_l, # ======= CALL ========== >>
                       1, pararg, xlim_hist, ylim_distr, hist_colr, ttl_ffh)

    #fig1d = pl.figure(); ax_pcd = pl.gca()
    ax_pcd = axd1['distr_pconv'] # Plot distr of PolConvert I param vs Time
    ttl_pcd = "PolConvert I, %s vs Time (%d triangles)" % (upar, ntri)
    plot_closures_distr(ax_pcd, timh, atau_c,  # ========== CALL ============ >>
                       1, cols, ylim_distr, pararg, ttl_pcd)

    #fig1c = pl.figure(); ax_pch = pl.gca()
    ax_pch = axd1['hist_pconv']  # Plot hist of PolConvert I param
    ttl_pch = "%s (%d points)" % (upar, nfinite)
    plot_closures_hist_horiz(ax_pch, timh, atau_c, # ======= CALL ========== >>
                        1, pararg, xlim_hist, ylim_distr, hist_colr, ttl_pch)
    print("atau_c.flatten().shape = ", atau_c.flatten().shape)








    
#============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#-------------------   With Y and without Y  --------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#============================================================================

#
# Separate tau's with and without station 'Y'
# and create triangle lists with and without station 'Y'
#
tau_l_sans_y = np.empty(0, dtype=float)
tau_l_with_y = np.empty(0, dtype=float)
tau_c_sans_y = np.empty(0, dtype=float)
tau_c_with_y = np.empty(0, dtype=float)
tau_l_flat = np.empty(0, dtype=float)
tau_c_flat = np.empty(0, dtype=float)

trians_with_y = []
trians_sans_y = []

for ic in range(ntri):
    trist = trians[ic]
    print('trist = ', trist)
    if 'Y' in trist:
        tau_l_with_y = np.append(tau_l_with_y, tau_l[trist])
        tau_c_with_y = np.append(tau_c_with_y, tau_c[trist])
        trians_with_y.append(trist)
    else:
        tau_l_sans_y = np.append(tau_l_sans_y, tau_l[trist])
        tau_c_sans_y = np.append(tau_c_sans_y, tau_c[trist])
        trians_sans_y.append(trist)

ntri_noY = len(trians_sans_y)
ntri_Y = len(trians_with_y)

print("trians_sans 'Y': ", trians_sans_y)
print("trians_with 'Y':    ", trians_with_y)




#
# Plot closures with and without station 'Y'
#

gs_kw2 = dict(width_ratios=[0.75, 0.25],
              height_ratios=[0.1, 0.225, 0.225, 0.225, 0.225])
              # height_ratios=[0.08, 0.23, 0.23, 0.23, 0.23])
fig2, axd2 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                               ['distr_frfit_sansY', 'hist_frfit_sansY'],
                               ['distr_pconv_sansY', 'hist_pconv_sansY'],
                               ['distr_frfit_withY', 'hist_frfit_withY'],
                               ['distr_pconv_withY', 'hist_pconv_withY']],
                               gridspec_kw=gs_kw2, figsize=(8.4, 10),
                               layout="constrained")
#
# Plot color legend on top
#
ax_col = axd2['col_legend']
plot_closure_legend(ax_col, trians, cols, upar, fs=10)  # ======  CALL ==== >>

timh = tim['TV']/3600  # Time counts (in hours) from baseline with no NaNs

# Non-NaNs with and without Y in atau_l and atau_c
nfinite_noY = np.count_nonzero(np.isfinite(atau_l[sel_noY,:]))
nfinite_Y = np.count_nonzero(np.isfinite(atau_l[sel_Y,:]))

hist_colr = 'red'

#
# Closure distributions and histograms without the Y station
#

if pararg == 'mbd':
    # ylim_distr = (-530, 620)  # It includes the outliers beyond +-400 ps
    ylim_distr = (-180, 180)
    xlim_hist = (-3, 100)  # Now it is 'xlim' after 90-deg rotation
elif pararg == 'sbd':
    ylim_distr = (-500, 500)
    xlim_hist = (-1, 40)
elif pararg == 'tmbd':
    ylim_distr = (-180, 180)
    xlim_hist = (-3, 100)  # Now it is 'xlim' after 90-deg rotation
elif pararg == 'tsbd':
    ylim_distr = (-500, 500)
    xlim_hist = (-1, 40)

ax_ffd = axd2['distr_frfit_sansY'] # Plot distr of FourFit pseudo-I no Y  
ttl_ffd = "Fourfit Pseudo-I, %s vs Time, no Y (%d triangles)" % \
                                                           (upar, ntri_noY) 
plot_closures_distr(ax_ffd, timh, atau_l,  # ============ CALL ============= >>
                   sel_noY, cols, ylim_distr, pararg, ttl_ffd)

ax_ffh = axd2['hist_frfit_sansY']  # Plot hist of FourFit I param
ttl_ffh = "%s (%d points)" % (upar, nfinite_noY)
plot_closures_hist_horiz(ax_ffh, timh, atau_l, # ========= CALL =========== >>
                sel_noY, pararg, xlim_hist, ylim_distr, hist_colr, ttl_ffh)

ax_pcd = axd2['distr_pconv_sansY'] # Plot distr of PolConvert pseudo-I no Y
ttl_pcd = "PolConvert I, %s vs Time, no Y (%d triangles)" % \
                                                        (upar, ntri_noY)
plot_closures_distr(ax_pcd, timh, atau_c,  # ============ CALL ============= >>
                   sel_noY, cols, ylim_distr, pararg, ttl_pcd)

ax_pch = axd2['hist_pconv_sansY']  # Plot hist of PolConvert I param
ttl_pch = "%s (%d points)" % (upar, nfinite_noY)
plot_closures_hist_horiz(ax_pch, timh, atau_c, # ========= CALL =========== >>
                sel_noY, pararg, xlim_hist, ylim_distr, hist_colr, ttl_pch)


#
# Closure distributions and histograms with the Y station only
#

if pararg == 'mbd':
    # ylim_distr = (-530, 620)  # It includes the outliers beyond +-400 ps
    ylim_distr = (-180, 180)
    xlim_hist = (-3, 100)  # Now it is 'xlim' after 90-deg rotation
elif pararg == 'sbd':
    ylim_distr = (-1200, 1000)
    xlim_hist = (-1, 40)
elif pararg == 'tmbd':
    ylim_distr = (-1e9, 1e9)
    xlim_hist = (-1, 60)  # Now it is 'xlim' after 90-deg rotation
elif pararg == 'tsbd':
    ylim_distr = (-1e9, 1e9)
    xlim_hist = (-1, 60)


ax_ffd = axd2['distr_frfit_withY'] # Plot distr of FourFit pseudo-I with Y only
ttl_ffd = "Fourfit Pseudo-I, %s vs Time, with Y only (%d triangles)" % \
                                                                (upar, ntri_Y)
plot_closures_distr(ax_ffd, timh, atau_l,  # ============ CALL ============= >>
                   sel_Y, cols, ylim_distr, pararg, ttl_ffd)

ax_ffh = axd2['hist_frfit_withY']  # Plot hist of FourFit I param
ttl_ffh = "%s (%d points)" % (upar, nfinite_Y)
plot_closures_hist_horiz(ax_ffh, timh, atau_l, # ========= CALL =========== >>
                sel_Y, pararg, xlim_hist, ylim_distr, hist_colr, ttl_ffh)

ax_pcd = axd2['distr_pconv_withY'] # Plot distr of PolConvert ps.-I with Y only
ttl_pcd = "PolConvert I, %s vs Time, with Y only (%d triangles)" % \
                                                                (upar, ntri_Y)
plot_closures_distr(ax_pcd, timh, atau_c,  # ============ CALL ============= >>
                   sel_Y, cols, ylim_distr, pararg, ttl_pcd)

ax_pch = axd2['hist_pconv_withY']  # Plot hist of PolConvert I param with Y only
ttl_pch = "%s (%d points)" %  (upar, nfinite_Y)
plot_closures_hist_horiz(ax_pch, timh, atau_c, # ========= CALL =========== >>
                sel_Y, pararg, xlim_hist, ylim_distr, hist_colr, ttl_pch)



pl.show()

#
# Save figures on request
#
if sf:
    # fig1 is not created for total_mbd or total_sbd; omit saving its pdf
    if 'fig1' in locals(): # In the current scope's local variables
        pl.figure(fig1)
        pl.savefig("%s_Closure_Delay.pdf" % upar, format='pdf')
        
    pl.figure(fig2)
    pl.savefig("%s_Closure_Delay_Y_no_Y.pdf" % upar, format='pdf')



