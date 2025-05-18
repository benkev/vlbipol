'''
libvp.py - a set of functions commomly used.

find_baseline_triangles(bls):
    Find all the baseline triplets from baseline list bls, such that the
    baseline order in them makes a correct closure
session_time_start(idx):
    Find session time start stim0 (time of the first scan) from the dictionary
    idx[bl][polprod][data_item]
make_closure_dic(idxs, bls=None):
    Create dictionary of all possible closures
        clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
    from the dictionary
        idxs[src][time][bl][data_name]
clos_to_clot(clos, tribl=None, bls=None):
    clos[sr][tr][di] --> clot[tr][sr][di], di - "data item"
    Rearrange a clos dict into clot by permuting the first two indices sr and tr
    Create dictionary of all possible closures with a triangle as first key
        clot[tri][src]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
    from the dictionary
        clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
    The clot dictionary contains the same data as clos, but the first two
    indices, sr and tr, are permuted. See make_closure_dic().  
find_tail_bounds(ni, thr=5):
    Find histogram ni indices, l_idx and r_idx, at which the sparse tails
    are to be cut off.
group_tails(ni, lr_inds):
    Group the tail data in histogram array ni: cut the tails and place
    what was in the tails to bins 0 and -1 (end).

'''

import os, sys, re
import copy, pickle
import numpy as np
from itertools import combinations
from bisect import bisect_right  # Bisection algorithm for efficient search



def find_baseline_triangles(bls):
    '''
    Find all the baseline triplets from baseline list bls, such that the
    baseline order in them makes a correct closure, like 'MT', 'TE', 'ME'.
    The last baseline is always in inverted order (here not 'EM' but 'ME'),
    so in the closure it must be present with the minus sign. For example,

        tau_MTE = tau_MT + tau_TE - tau_ME           - for closure delay
    or
        cloph__MTE = cloph_MT + cloph_TE - cloph_ME  - for closure phase
    
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
    ###  print("The baseline triplets are reordered:")

    tribl1 = {}
    for trist in tribl.keys():
        l3bl = tribl[trist]  # List of the 3 baselines
        trian = l3bl

        st_end = l3bl[0][1]
        if l3bl[2][0] == st_end:
            trian = (l3bl[0], l3bl[2], l3bl[1])
            tribl1[trist] = trian
            ### print(tribl[trist], '->', trian)

        st_end = l3bl[1][1]
        if l3bl[0][0] == st_end:
            trian = (l3bl[1], l3bl[0], l3bl[2])
            tribl1[trist] = trian
            ### print(tribl[trist], '->', trian)
        elif l3bl[2][0] == st_end:
            trian = (l3bl[1], l3bl[2], l3bl[0])
            tribl1[trist] = trian
            ### print(tribl[trist], '->', trian)

    tribl = tribl1

    return tribl




def session_time_start(idx):
    '''
    Find session time start stim0 (time of the first scan)
    '''
    bls = list(idx.keys())    # All the baselines
    tim_total = []     # List of all the time counts for all the baselines

    for bl in bls:
        timbl = idx[bl]['I']['time_tag']
        tim_total.extend(timbl)

    ttim = np.unique(tim_total)   # Unite all the time counts in one array
    ttim.sort()
    stim0 = ttim[0]    # Session time start (the fitst scan)

    return stim0





def make_closure_dic(idxs, bls=None):
    '''
    Create dictionary of all possible closures
        clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
    from the dictionary
        idxs[src][time][bl][data_name]

    The bls parameter is a list of allowed baselines. If not given,
    all the baselines are involved. In the VO2187 experiment, for example,
    the baseline ST and all the baselines with station Y are excluded.

    The clos dictionary returned contains not only the closures, but also
    the data triplets used to compute the closures. All the numeric data are
    in arrays sorted in time ascending order. For example,

    closl = make_closure_dic(idxsl, bls)     # Use linear polarization data 

    closl['0955+476']['EGM']['time'] -->    (s, time from the experiment start)
        array([ 2248.,  3619.,  6700.,  9360., 10685., 13321., 21292., 22588.,
                23881., 25343., 79412., 80982., 83646.])
    
    closl['0955+476']['EGM']['thour'] --> (hours, time from the experim. start)
        array([  0.6244,   1.0053,   1.8611,   2.6,   2.9681, 3.7003,   5.9144,
                 6.2744,   6.6336,   7.0397,  22.0589,  22.495,  23.235])
    
    closl['0955+476']['EGM']['cloph'] -->                  (deg, closure phase)
        array([-15.1 ,   6.59,   1.01,   0.,   9.22,   2.56,  19.77,  11.09,
                10.36,   5.06,  17.94, -18.66,  11.66])

    closl['0955+476']['EGM']['phase'] -->       (deg, phases used for closures)
        array([[146.61, 176.44, 338.15],
               [154.64, 201.7 , 349.76],
               [ 12.88, 208.22, 220.09],
               [ 25.11,  62.5 ,  87.61],
               [298.5 ,  68.46, 357.74],
               [113.03,  33.49, 143.97],
               [290.49, 255.25, 165.97],
               [318.21, 215.77, 162.88],
               [307.54,  87.32,  24.49],
               [330.73, 330.7 , 296.37],
               [ 29.32, 304.32, 315.71],
               [102.33, 166.78, 287.77],
               [305.16, 293.44, 226.94]])

    closl['0955+476']['EGM']['bl'] -->             (baselines of EGM triangle )
        [('GM', 'ME', 'GE')]
    
    closl['0955+476']['EGM']['file'] -->                (files data taken from)
        [['GM.X.2.3HJQDH', 'ME.X.1.3HJQDH', 'GE.X.6.3HJQDH'],
        ['GM.X.1.3HJQFI', 'ME.X.2.3HJQFI', 'GE.X.3.3HJQFI'],
        ['GM.X.1.3HJQKN', 'ME.X.4.3HJQKN', 'GE.X.6.3HJQKN'],
        ['GM.X.1.3HJQP5', 'ME.X.4.3HJQP5', 'GE.X.5.3HJQP5'],
        ['GM.X.2.3HJQRH', 'ME.X.1.3HJQRH', 'GE.X.5.3HJQRH'],
        ['GM.X.2.3HJQW6', 'ME.X.7.3HJQW6', 'GE.X.28.3HJQW6'],
        ['GM.X.2.3HJRA1', 'ME.X.5.3HJRA1', 'GE.X.27.3HJRA1'],
        ['GM.X.2.3HJRC6', 'ME.X.5.3HJRC6', 'GE.X.28.3HJRC6'],
        ['GM.X.1.3HJRE6', 'ME.X.4.3HJRE6', 'GE.X.15.3HJRE6'],
        ['GM.X.2.3HJRGM', 'ME.X.4.3HJRGM', 'GE.X.10.3HJRGM'],
        ['GM.X.2.3HJU22', 'ME.X.1.3HJU22', 'GE.X.3.3HJU22'],
        ['GM.X.1.3HJU4K', 'ME.X.2.3HJU4K', 'GE.X.3.3HJU4K'],
        ['GM.X.1.3HJU95', 'ME.X.3.3HJU95', 'GE.X.2.3HJU95']]

    All the data items available:
    
    closl['0955+476']['EGM'].keys()
        dict_keys(['bl', 'time', 'thour', 'time_tag', 'cloph', 'tau_mbd',
                   'tau_sbd', 'tau_tmbd', 'tau_tsbd', 'phase', 'mbd', 'sbd',
                   'tmbd', 'tsbd', 'snr', 'pol_prod', 'file', 'dir'])
    
    In the code below, the following variables are introduced for brevity:
        xst = idxs[sr][tm]
        cst = clos[sr][tr]
    '''
    
    clos = {}

    for sr in idxs.keys():
        for tm in idxs[sr].keys():

            # Find baselines for the current source and time
            xbls_all = list(idxs[sr][tm].keys()) 

            if bls:
                xbls = []
                for bl in xbls_all:
                    if bl in bls: xbls.append(bl)
            else:
                xbls = xbls_all        

            if len(xbls) < 3: continue  # === Ignore less than 3 bls === >>

            # Find baseline triangles for 3 or more baselines in xbls
            xtris = find_baseline_triangles(xbls)

            xst = idxs[sr][tm]

            for tr in xtris.keys():
                #
                # For this source and this time, save triads of parameters
                # involved in calculation of closures`
                #
                bl3 = xtris[tr]    # list ot 3 bls making up a triangle tr
                ph3 =   [xst[bl]['phase'] for bl in bl3]
                dtec3 = [xst[bl]['dtec'] for bl in bl3]
                mbd3 =  [xst[bl]['mbdelay'] for bl in bl3]
                sbd3 =  [xst[bl]['sbdelay'] for bl in bl3]
                tmbd3 = [xst[bl]['tot_mbd'] for bl in bl3]
                tsbd3 = [xst[bl]['tot_sbd'] for bl in bl3]
                snr3 =  [xst[bl]['snr'] for bl in bl3]
                fl3 =   [xst[bl]['file'] for bl in bl3]
                dir3 =  [xst[bl]['dir'] for bl in bl3]
                pp3 =   [xst[bl]['pol_prod'] for bl in bl3]
                pp = pp3[0]
                ttag = xst[bl3[0]]['time_tag']

                #
                # Compute all the possible closures
                #
                ab, bc, ac = xtris[tr]   # 3 bls making up a triangle tr
                cloph = xst[ab]['phase'] + xst[bc]['phase'] - \
                        xst[ac]['phase']
                cloph = ((cloph + 180) % 360) - 180  # Reduce to [-180 .. +180]
                tau_mbd = xst[ab]['mbdelay'] + \
                          xst[bc]['mbdelay'] - \
                          xst[ac]['mbdelay']
                tau_sbd = xst[ab]['sbdelay'] + \
                          xst[bc]['sbdelay'] - \
                          xst[ac]['sbdelay']
                tau_tmbd = xst[ab]['tot_mbd'] + \
                          xst[bc]['tot_mbd'] - \
                          xst[ac]['tot_mbd']
                tau_tsbd = xst[ab]['tot_sbd'] + \
                          xst[bc]['tot_sbd'] - \
                          xst[ac]['tot_sbd']
                if sr in clos.keys():
                    if tr in clos[sr].keys():
                        cst = clos[sr][tr]
                        #
                        # Find index insr into the time array using fast 
                        # dichotomy (or bisection) algorithm.
                        # The insr index points at the location to insert the
                        # time tag keeping the time ascending order.
                        #
                        insr = bisect_right(cst['time'], tm)

                        # cst['bl'].insert(insr, bl3) # Insert bls only once!
                        cst['time'] = np.insert(cst['time'], insr, tm)
                        cst['thour'] = np.insert(cst['thour'], insr, tm/3600)
                        cst['time_tag'] = np.insert(cst['time_tag'],
                                                    insr, ttag, 0)
                        cst['cloph'] = np.insert(cst['cloph'], insr, cloph)
                        cst['tau_mbd'] = np.insert(cst['tau_mbd'],insr, tau_mbd)
                        cst['tau_sbd'] = np.insert(cst['tau_sbd'],insr, tau_sbd)
                        cst['tau_tmbd'] = np.insert(cst['tau_tmbd'],
                                                    insr, tau_tmbd)
                        cst['tau_tsbd'] = np.insert(cst['tau_tsbd'],
                                                    insr, tau_tsbd)
                        cst['phase'] = np.insert(cst['phase'], insr, ph3, 0)
                        cst['dtec'] = np.insert(cst['dtec'], insr, dtec3, 0)
                        cst['mbd'] = np.insert(cst['mbd'], insr, mbd3, 0)
                        cst['sbd'] = np.insert(cst['sbd'], insr, sbd3, 0)
                        cst['tmbd'] = np.insert(cst['tmbd'], insr, tmbd3, 0)
                        cst['tsbd'] = np.insert(cst['tsbd'], insr, tsbd3, 0)
                        cst['snr'] = np.insert(cst['snr'], insr, snr3, 0)
                        #cst['pol_prod'].insert(insr, pp) # Insert pp only once!
                        cst['file'].insert(insr, fl3)
                        cst['dir'].insert(insr, dir3)
                    else:
                        clos[sr][tr] =  {
                            'bl':[bl3], # Insert bls only once!
                            'time': np.array([tm], dtype=float),
                            'thour': np.array([tm], dtype=float)/3600,
                            'time_tag': np.array([ttag], dtype=float),
                            'cloph': np.array([cloph], dtype=float),
                            'tau_mbd': np.array([tau_mbd], dtype=float),
                            'tau_sbd': np.array([tau_sbd], dtype=float),
                            'tau_tmbd': np.array([tau_tmbd], dtype=float),
                            'tau_tsbd': np.array([tau_tsbd], dtype=float),
                            'phase': np.array([ph3], dtype=float),
                            'dtec': np.array([dtec3], dtype=float),
                            'mbd': np.array([mbd3], dtype=float),
                            'sbd': np.array([sbd3], dtype=float),
                            'tmbd': np.array([tmbd3], dtype=float),
                            'tsbd': np.array([tsbd3], dtype=float),
                            'snr': np.array([snr3], dtype=float),
                            'pol_prod':[pp], # Insert polprod only once!
                            'file':[fl3],
                            'dir':[dir3]
                        }
                else:
                    clos[sr] = {}
                    clos[sr][tr] = {
                        'bl':[bl3], # Insert bls only once!
                        'time': np.array([tm], dtype=float),
                        'thour': np.array([tm], dtype=float)/3600,
                        'time_tag': np.array([ttag], dtype=float),
                        'cloph': np.array([cloph], dtype=float),
                        'tau_mbd': np.array([tau_mbd], dtype=float),
                        'tau_sbd': np.array([tau_sbd], dtype=float),
                        'tau_tmbd': np.array([tau_tmbd], dtype=float),
                        'tau_tsbd': np.array([tau_tsbd], dtype=float),
                        'phase': np.array([ph3], dtype=float),
                        'dtec': np.array([dtec3], dtype=float),
                        'mbd': np.array([mbd3], dtype=float),
                        'sbd': np.array([sbd3], dtype=float),
                        'tmbd': np.array([tmbd3], dtype=float),
                        'tsbd': np.array([tsbd3], dtype=float),
                        'snr': np.array([snr3], dtype=float),
                        'pol_prod':[pp], # Insert polprod only once!
                        'file':[fl3],
                        'dir':[dir3]
                    }
    return clos



def clos_to_clot(clos, tribl=None, bls=None):
    '''
    clos[sr][tr][di] --> clot[tr][sr][di], di - "data item"

    Rearrange a clos dict into clot by permuting the first two indices sr and tr
    
    Create dictionary of all possible closures with a triangle as first key
        clot[tri][src]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
    from the dictionary
        clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]

    The tribl parameter is a dictionary tribl[tr] --> (bl1, bl2, bl3)
    If not provided, it is created from the bls list.
    
    The bls parameter is a list of allowed baselines.

    If neither tribl nor bls are provided, tribl is loaded from disk,
    from file tribl_2107.pkl.

    The returned clot dictionary contains not only the closures, but also
    the data triplets used to compute the closures including Mark4 fringe-fit
    file directories and file names.
    All the numeric data are in arrays sorted in time ascensing order.

    The clot dictionary contains the same data as clos, but the first two
    indices, sr and tr, are permuted. See make_closure_dic().  
    '''
      
    if tribl and bls:   # Both tribl and bls provided: use tribl only
        pass
    elif bls:  # Only bls provided: compute tribl from bls:
        tribl = find_baseline_triangles(bls)
    else:      # Neither tribl nor bls provided:
        with open('tribl_2187.pkl', 'rb') as finp:  # Load tribl from disk
            tribl = pickle.load(finp)
        
    clot = {tr: {} for tr in tribl.keys()} 

    for sr in clos.keys():
        for tr in clos[sr].keys():
            for di in clos[sr][tr].keys():
                if sr in clot[tr].keys():
                    clot[tr][sr][di] = copy.deepcopy(clos[sr][tr][di])
                else:
                    clot[tr][sr] = {}
                    clot[tr][sr][di] = copy.deepcopy(clos[sr][tr][di])

    return clot




def find_tail_bounds(ni, thr=5):
    '''
    Find histogram ni indices, l_idx and r_idx, at which the sparse tails of ni
    are to be cut off.
      ni: histogram array with frequencies
      thr: frequency threshold, below which the tail is considered "sparse".
    Returns l_idx and r_idx, where
      l_idx points AT the first ni element greater than thr;
      r_idx points AFTER the last ni element greater than thr,
   so ni[l_idx : r_idx] is the subarray whose elements are > thr.
    
    '''
    nbin = len(ni)

    l_idx = 0
    r_idx = nbin

    for i in range(nbin):
        if ni[i] > thr:
            l_idx = i
            break

    for i in range(nbin-1,-1,-1): # nbin-1 downto 0 
        if ni[i] > thr:
            r_idx = i + 1
            break

    return l_idx, r_idx




def group_tails(ni, lr_inds):
    '''
    Group the tail data in histogram array ni:
    cut the tails and place what was in the tails to
    bins 0 and -1 (end).
      ni: histogram array with frequencies
      lr_inds: a sequence of two indices, (l_idx, r_idx), left and right.
      l_idx points AT the first ni element greater than thr;
      r_idx points AFTER the last ni element greater than thr,
    so ni[l_idx : r_idx] is the subarray whose elements are > thr.
    
    Returns ni with the sparse tails grouped.
    '''

    nbin = len(ni)
    l_idx, r_idx = lr_inds
    
    ni_grp = np.copy(ni[l_idx : r_idx])   

    if l_idx > 0:
        l_ni = ni[:l_idx+1]
        ni_grp[0] = np.sum(l_ni)

    if r_idx < nbin:
        r_ni = ni[r_idx-1:]
        ni_grp[-1] = np.sum(r_ni)

    return ni_grp














