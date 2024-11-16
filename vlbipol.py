help_text = \
'''
vlbipol.py - a module containing the class vlbipol to compare
             fourfit Pseudo-I and PolConvert products.
'''

import pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib.patches as patches
from itertools import combinations
import copy

pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.


class vlbipol:

    #
    # Unpickle the dictionaries with fourfit and PolConvert priducts
    #
    with open('idx3819l.pkl', 'rb') as finp:
        idx3819l_1 = pickle.load(finp)

    with open('idx3819cI.pkl', 'rb') as finp:
        idx3819c_1 = pickle.load(finp)
    del finp                     # Keep the namespace as clean as possible

    # arg_to_par = {'mbd':'mbdelay', 'sbd':'sbdelay', 'tmbd':'tot_mbd',
    #               'tsbd':'tot_sbd', 'rmbd':'resid_mbd', 'rsbd':'resid_sbd'}
    arg_to_par = {'mbd':'mbdelay', 'sbd':'sbdelay'}
    istart = 2 # Start processing from istart;  exclude bad data before istart.
    ps = "(ps)"

    bls = list(idx3819l_1.keys())   # Baselines
    bls.sort()                      # Lexigraphically sorted baselines
    nbls = len(bls)

    #
    # Find tribl, the dictionary with all the baseline triplets with 3 stations
    # (ie station triangles) such that their order makes a correct closure,
    # like 'MT', 'TE', 'ME'.
    # The last baseline is always in inverted order (here not 'EM' but 'ME'),
    # so in the closure delay it must be present with the minus sign.
    # For example,
    #
    #        tau_MTE = tau_MT + tau_TE - tau_ME
    #

    tribl = find_baseline_triangles(bls)   # =========== CALL ============== >>

    trians = list(tribl.keys())
    ntri = len(tribl)

    print()
    print('bls = ', bls, '\n')
    print('tribl:')
    for trist in tribl.keys():
        print(trist, " -> ", tribl[trist])


    @classmethod
    def closure_residual_delay(cls)
    #
    # Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or 'snr' or
    #     'tmbd' or 'tsbd' or 'rmbd' or 'rsbd'
    #
    parname = cls.arg_to_par[prm]


    #sys.exit(0)

    #
    # To start processing from istart;  exclude bad data before istart.
    #
    istart = cls.istart

    tim1 = {}     # Original time points with some of them missing. 
    tim = {}      # Time points. The gaps will be replaced with NaNs

    #
    # Search for the maximum number of time counts, ntim,
    # among all the baselines. At the end of the loop below, ntim will hold
    # the length for the parameter storage in arrays.
    # The parameter sets shorter than ntim contain time gaps to be filled with NaNs.
    #
    ntim = -1 
    for bl in bls:
        tim1[bl] = np.array(idx3819l_1[bl]['I']['time'])[istart:] #/ 60 # Sec -> min
        tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

        if len(tim1[bl]) > ntim:
            ntim = len(tim1[bl])

        print("len(tim1['%s']) = %d" % (bl, len(tim1[bl])))



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
        # Insert NaNs in the time gaps (ie over 605. seconds without a time count)
        # Accordingly, insert NaNs in the parameter arrays
        #
        itim = np.int64(tim1[bl]/605) # Indices of non-NaN elements into tim and par
        tim[bl] = np.zeros(ntim)    
        tim[bl][:] = np.NaN
        tim[bl][itim] = tim1[bl]
        atim[itri,:] = tim[bl]

        snr_l[bl] = np.zeros(ntim)    
        snr_l[bl][:] = np.NaN
        snr_l[bl][itim] = snr1_l
        asnr_l[itri,:] = snr_l[bl]

        snr_c[bl] = np.zeros(ntim)    
        snr_c[bl][:] = np.NaN
        snr_c[bl][itim] = snr1_c
        asnr_c[itri,:] = snr_c[bl]

        par_l[bl] = np.zeros(ntim)    
        par_l[bl][:] = np.NaN
        par_l[bl][itim] = par1_l
        apar_l[itri,:] = par_l[bl]

        par_c[bl] = np.zeros(ntim)    
        par_c[bl][:] = np.NaN
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




    ylms = -50000000   # +- ylimits, if ylms > 0

    # gs_kw1 = dict(width_ratios=[1, 1], height_ratios=[0.15, 0.6, 0.25])
    # fig1, axd1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
    #                                ['distr_frfit', 'distr_pconv'],
    #                                ['hist_frfit', 'hist_pconv']],
    #                                gridspec_kw=gs_kw1, figsize=(8.4, 8),
    #                                layout="constrained")

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

    plot_closure_legend(ax_col, trians, cols, upar)  # =========  CALL ========= >>



    if prm == 'mbd':
        # ylim_distr = (-530, 620)  # It includes the outliers beyond +-400 ps
        ylim_distr = (-220, 220)
        xlim_hist = (-5, 170)  # Now it is 'xlim' after 90-deg rotation
    elif prm == 'sbd':
        ylim_distr = (-1990, 870)
        xlim_hist = (-10, 1000)



    timh = tim['TV']/3600   # Time counts (in hours) from baseline with no NaNs

    hist_colr = 'red'

    #fig1b = pl.figure(); ax_ffd = pl.gca()
    ax_ffd = axd1['distr_frfit'] # Plot distr of FourFit pseudo-I param vs Time  
    ttl_ffd = "Fourfit Pseudo-I, %s vs Time" % upar
    plot_closures_dist(ax_ffd, timh, atau_l,  # ============ CALL ============= >>
                       1, cols, ylim_distr, prm, ttl_ffd)

    #fig1a = pl.figure(); ax_ffh = pl.gca()
    ax_ffh = axd1['hist_frfit']  # Plot hist of FourFit I param
    ttl_ffh = "%s" % upar
    plot_closures_hist_horiz(ax_ffh, timh, atau_l, # ========= CALL =========== >>
                             1, prm, xlim_hist, ylim_distr, hist_colr, ttl_ffh)

    #fig1d = pl.figure(); ax_pcd = pl.gca()
    ax_pcd = axd1['distr_pconv'] # Plot distr of PolConvert  pseudo-I param vs Time
    ttl_pcd = "PolConvert I, %s vs Time" % upar
    plot_closures_dist(ax_pcd, timh, atau_c,  # ============ CALL ============= >>
                       1, cols, ylim_distr, prm, ttl_pcd)

    #fig1c = pl.figure(); ax_pch = pl.gca()
    ax_pch = axd1['hist_pconv']  # Plot hist of PolConvert I param
    ttl_pch = "%s" % upar
    plot_closures_hist_horiz(ax_pch, timh, atau_c, # ========= CALL =========== >>
                             1, prm, xlim_hist, ylim_distr, hist_colr, ttl_pch)





    sys.exit(0)





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

    print("trians_sans 'Y': ", trians_sans_y)
    print("trians_with 'Y':    ", trians_with_y)

    # #
    # # 
    # #

    # sel_Y = np.zeros(ntri, dtype='bool')
    # sel_noY = np.zeros(ntri, dtype='bool')

    # for itri in range(ntri):
    #     if 'Y' in trians[itri]:
    #         sel_Y[itri] = True
    #         print("Y in trians[%d] = %s" % (itri, trians[itri]))
    #     else:
    #         sel_noY[itri] = True
    #         print("Y not in trians[%d] = %s" % (itri, trians[itri]))





    #!!!!!!!!!!!! ????????????????? 






    #
    # Plot closures with and without station 'Y'
    #
    gs_kw2 = dict(width_ratios=[1, 1],
                  height_ratios=[0.12, 0.3, 0.14, 0.3, 0.14])
    fig2, axd2 = pl.subplot_mosaic([['col_legend',       'col_legend'],
                                   ['distr_frfit_sansY', 'distr_pconv_sansY'],
                                   ['hist_frfit_sansY',  'hist_pconv_sansY'],
                                   ['distr_frfit_withY', 'distr_pconv_withY'],
                                   ['hist_frfit_withY',  'hist_pconv_withY']],
                                   gridspec_kw=gs_kw2, figsize=(8.4, 10),
                                   layout="constrained")
    #
    # Plot color legend on top
    #
    ax_col = axd2['col_legend']

    plot_closure_legend(ax_col, trians, cols, upar)  # =========  CALL ========= >>

    timh = tim['TV']/3600  # Time counts (in hours) from baseline with no NaNs

    #
    # Plot distribution and histogram of mbdelay closures without station Y
    #

    ax_ffd = axd2['distr_frfit_sansY']  # Plot distr of FourFit pseudo-I param no Y

    for ic in range(ntri):
        trist = trians[ic]
        if trist in trians_sans_y:
            ax_ffd.plot(timh, tau_l[trist], '.', color=cols[ic,:])
    ax_ffd.grid(1)
    ax_ffd.set_title("Fourfit Pseudo-I, %s vs Time" % upar)
    ax_ffd.set_xlabel("hours", fontsize=14)
    ax_ffd.set_ylabel("ps", fontsize=14)
    ax_ffd.yaxis.set_label_coords(-0.05, 0.58)
    ax_ffd.xaxis.set_label_coords(0.55, 0.07)
    if ylms > 0: ax_ffd.set_ylim(-ylms, ylms)
    if prm == 'mbd':
        ax_ffd.set_ylim(-530, 620)
    elif prm == 'sbd':
        ax_ffd.set_ylim(-1990, 870)


    ax_ffh = axd2['hist_frfit_sansY']  # Plot hist of FourFit I param

    ax_ffh.hist(abs(tau_l_sans_y), 50)
    ax_ffh.grid(1)
    ax_ffh.set_xlabel("ps", fontsize=14)
    ax_ffh.set_title("Fourfit Pseudo-I, abs(%s) sans Y" % upar)
    ax_ffh.xaxis.set_label_coords(0.5, -0.12)




    ax_ffd = axd2['distr_frfit_withY']  # Plot distr of FourFit pseudo-I param w/Y

    for ic in range(ntri):
        trist = trians[ic]
        ax_ffd.plot(timh, tau_c[trist], '.', color=cols[ic,:])
    ax_ffd.grid(1)
    ax_ffd.set_title("Fourfit Pseudo-I, %s vs Time" % upar)
    ax_ffd.set_xlabel("hours", fontsize=14)
    ax_ffd.set_ylabel("ps", fontsize=14)
    ax_ffd.yaxis.set_label_coords(-0.05, 0.58)
    ax_ffd.xaxis.set_label_coords(0.55, 0.07)
    if ylms > 0: ax_ffd.set_ylim(-ylms, ylms)
    if prm == 'mbd':
        ax_ffd.set_ylim(-530, 620)
    elif prm == 'sbd':
        ax_ffd.set_ylim(-1990, 870)


    ax_ffh = axd2['hist_frfit_withY']  # Plot hist of FourFit I param

    ax_ffh.hist(abs(atau_c.flatten()), 50)
    ax_ffh.grid(1)
    ax_ffh.set_xlabel("ps", fontsize=14)
    ax_ffh.set_title("Fourfit Pseudo-I, abs(%s)" % upar)
    ax_ffh.xaxis.set_label_coords(0.5, -0.12)






    ax_pcd = axd2['distr_pconv_sansY'] # Plot distr of PolConv ps-I param no Y

    for ic in range(ntri):
        trist = trians[ic]
        ax_pcd.plot(timh, tau_l[trist], '.', color=cols[ic,:])
    ax_pcd.grid(1)
    ax_pcd.set_title("PolConvert I, %s vs Time" % upar)
    ax_pcd.set_xlabel("hours", fontsize=14)
    ax_pcd.set_ylabel("ps", fontsize=14)
    ax_pcd.yaxis.set_label_coords(-0.05, 0.58)
    ax_pcd.xaxis.set_label_coords(0.55, 0.07)
    if prm == 'mbd':
        ax_pcd.set_ylim(-530, 620)
    elif prm == 'sbd':
        ax_pcd.set_ylim(-1990, 870)


    ax_pch = axd2['hist_pconv_sansY']  # Plot hist of PolConvert I param

    ax_pch.hist(abs(atau_l.flatten()), 50)
    ax_pch.grid(1)
    ax_pch.set_xlabel("ps", fontsize=14)
    ax_pch.set_title("PolConvert I, abs(%s)" % upar)
    ax_pch.xaxis.set_label_coords(0.5, -0.12)





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


    pl.show()

    #
    # Save figures on request
    #
    if sf:
        pl.figure(fig1)
        pl.savefig("%s_Closure_Delay.pdf" % upar, format='pdf')
        pl.figure(fig2)
        pl.savefig("%s_Closure_Delay_y_no_y.pdf" % upar, format='pdf')


