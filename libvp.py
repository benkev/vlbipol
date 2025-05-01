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
make_idx(base_dir, pol='lin', max_depth=2):
    Returns an index dictionary for all the fringe-fit files found in 
    any directory under the base_directory at max_depth.
    pol: 'lin' - linear polarization, 'cir' - circular polarization

'''

import os, sys, re
import copy, pickle
from itertools import combinations
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib.patches as patches
from itertools import combinations
from bisect import bisect_right  # Bisection algorithm for efficient search
#from group_tails import find_tail_bounds, group_tails

import hopstestb as ht
import ffcontrol
from vpal import fringe_file_manipulation as ffm
from vpal.utility import int_to_time, time_to_int
import mk4b



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

    The returned the clos dictionary contains not only the closures, but also
    the data triplets used to compute the closures. All the numeric data are
    in arrays sorted in time ascensing order. For example,

    closl = make_closure_dic(idxsl, bls)     # Use linear polarization data 

    closl['0955+476']['EGM']['time'] -->    (s, time from the experiment start)
        array([ 2248.,  3619.,  6700.,  9360., 10685., 13321., 21292., 22588.,
                23881., 25343., 79412., 80982., 83646.])
    
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
        dict_keys(['bl', 'time', 'time_tag', 'cloph', 'tau_mbd', 'tau_sbd',
                   'tau_tmbd', 'tau_tsbd', 'phase', 'mbd', 'sbd', 'tmbd',
                   'tsbd', 'snr', 'pol_prod', 'file', 'dir'])
    
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
                            'time_tag': np.array([ttag], dtype=float),
                            'cloph': np.array([cloph], dtype=float),
                            'tau_mbd': np.array([tau_mbd], dtype=float),
                            'tau_sbd': np.array([tau_sbd], dtype=float),
                            'tau_tmbd': np.array([tau_tmbd], dtype=float),
                            'tau_tsbd': np.array([tau_tsbd], dtype=float),
                            'phase': np.array([ph3], dtype=float),
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
                        'time_tag': np.array([ttag], dtype=float),
                        'cloph': np.array([cloph], dtype=float),
                        'tau_mbd': np.array([tau_mbd], dtype=float),
                        'tau_sbd': np.array([tau_sbd], dtype=float),
                        'tau_tmbd': np.array([tau_tmbd], dtype=float),
                        'tau_tsbd': np.array([tau_tsbd], dtype=float),
                        'phase': np.array([ph3], dtype=float),
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




def make_idx(base_dir, pol='lin', max_depth=2):
    '''
    Returns an index dictionary for all the fringe-fit files found in 
    any directory under the base_directory at max_depth.

    pol: 'lin' - linear polarization
         'cir' - circular polarization

     Creates and returns dictionaries to keep Mark4 data in convenient
     format. Returns the following dictionaries:

             1. idx[baseline][polproduct][data_item]
             
             is a dictionary of baselines, each baseline being a
             subdictionary of polproducts. The polproduct is in turn a
             subdictionary with data_item keys for 'time', 'file',
             'mbdelay' etc.
             The value for each lowest-level key is a list or an array. 
             The 'time' key points at the array of times in ascensing order.
                        Values in 'time' are counted from the experiment start.
             The 'file' key points at the list of file names corresponding to
                        the times in the 'time' list.
             Other keys point to their respective lists and arrays also
             strictly in ascending time order.

             Accessing the lists
             For example, the files, times, and more data for baseline 'EV'
                 are in subdictionary idx['EV'].
             The files, times, and more data for baseline 'EV' and polarization
                 'XY' are in subdictionary idx['EV']['XY'].
             idx['EV']['XY']['time'] is the list of time tags for 'EV' and 'XY'
             idx['EV']['XY']['file'] is the list of file names for 'EV' and 'XY'
             idx['EV']['XY']['mbdelay'] is the list of multiband delays
                                        for 'EV' and 'XY'.
             idx['EV']['XY']['sbdelay'] is the list of single-band delays
                                        for 'EV' and 'XY'.
             idx['EV']['XY']['snr'] is the list of signal to noise ratios
                                        for 'EV' and 'XY'.

             
             2. idxs[source][time][baseline][data_item]

             is a dictionary of celestial sources, each being a dictionary
             of their observation times, each being a dictionary of 
             the baselines involved. Each baseline is a subdictionary of
             data items. For example:

             idxsl['1803+784'][4065.0]['IT'] ->
                 {'mbdelay': -291.9238177128136,
                  'sbdelay': 273.9999908953905,
                  'tot_mbd': -6427.2611338087445,
                  'tot_sbd': -6427.260567884917,
                  'phase': 193.77931213378906,
                  'snr': 112.79531860351562,
                  'pol_prod': 'I',
                  'dir': '187-1907a',
                  'file': 'IT.X.2.3HJQG3',
                  'full_fname': '/home/benkev/Work/2187/scratch/Lin_I/2187/' \
                                '187-1907a/IT.X.2.3HJQG3',
                  'time_tag': 1341601665.5}


             3. idxf[directory][file]

             is a dictionary of Mark4 data directories, each being a dictionary
             of files therein, each being a dictionary of data items.
             For example:

             idxfl['188-0435a']['HT.X.4.3HJS31'] -->
                 {'source': '1803+784',
                  'time': 38118.0,
                  'bl': 'HT',
                  'mbdelay': -2092.3358388245106,
                  'sbdelay': -3383.500035852194,
                  'phase': 100.13414764404297,
                  'snr': 141.68231201171875,
                  'pol_prod': 'I',
                  'tot_mbd': -8538.048102473467,
                  'tot_sbd': -8538.04939363753,
                  'full_fname': '/home/benkev/Work/2187/scratch/Lin_I/2187/' \
                                '188-0435a/HT.X.4.3HJS31',
                  'time_tag': 1341635718.5}
    '''

    #
    # NOTE:
    #       mf.t208.contents.resid_mbd is the same as f_obj.mbdelay
    #       mf.t208.contents.resid_sbd is the same as f_obj.sbdelay
    # Therefore, adding resid_mbd and resid_sbd to the index makes no sense.
    # The lines with them have been commented out.
    #
    
    assert os.path.isdir(base_dir)

    # Remove trailing "/", if any (os.path.sep is usually "/")
    base_dir = base_dir.rstrip(os.path.sep)
    # base_dir = os.path.abspath(base_dir)
    num_sep = base_dir.count(os.path.sep)
    # ff_list = []

    #
    # Polarization notation table
    #
    lin2cir = {'XX':'LL', 'XY':'LR', 'YX':'RL', 'YY':'RR'}

    idx = {}
    idxs1 = {}      # To be a dict [src][time][bl][dn] with unsorted time
    idxf= {}
    
    # print("base_dir = ", base_dir)

    dir_line = []   # List to accumulate 10 dirnames for progress print
    n_linepad = 10  # Number of dirnames to print in one line
    i_linepad = 0    

    for root_dir, subdirs, files in os.walk(base_dir):

        # Apply the max depth filter
        if root_dir.count(os.path.sep) > num_sep + max_depth:
            continue

        f_obj = ffm.FringeFileHandle()
        
        for file in files:
            # print("root_dir = ", root_dir)
            # print("file = ", file)
            
            #abs_filename = os.path.abspath(filename)
            #filename_base = os.path.split(abs_filename)[1]
            
            #
            # Only if the file matches the pattern, its name is 
            # assigned to filename. The pattern is: two baseline capital 
            # letters, dot, letter X, dot, one ore more digits, dot,
            # six alphanumeric characters.
            #
            filename = re.findall(r"^[A-Z]{2}\.X\.[0-9]+\.\w{6}$", file)
            if filename == []:
                continue

            filename = filename[0]
            full_name = os.path.join(root_dir, filename)
            dir_name = root_dir.split('/')[-1]

            bl = filename[:2]  # Baseline is first two letters of filename

            #
            # Exclude (occasional) autocorrelations
            #
            if bl[0] == bl[1]:
                continue

            # print("full_name = ", full_name)
            # print("dir_name = ", dir_name, ", filename = ", filename)
            
            pp_list = ht.get_file_polarization_product_provisional(full_name)

            if len(pp_list) == 1:
                pp = pp_list[0]
                if pol == 'cir': # For circular polarization change X->L, Y->R
                    pp = lin2cir[pp]
            elif pp_list == ['XX', 'YY']: # For circular polarization only
                pp = 'I'
            else:
                continue

            f_obj.load(full_name)
            src = f_obj.source
            phase = f_obj.resid_phas

            #
            # Check the source correctness using regex
            #
            # mobj = re.fullmatch(r"^[0-9A-Z+-]{5,8}$", src) # Match object
            # if mobj is None:
            #     print("+++ Source '%s' does not match parretn! +++" % src)
            
            ttag = f_obj.time_tag         # Float, time or measurement 
            mbdelay = f_obj.mbdelay*1e6   # Float, multiband delay, us -> ps 
            sbdelay = f_obj.sbdelay*1e6   # Float, single-band delay, us -> ps 
            snr = f_obj.snr               # Float, signal to noise ratio

            mf = mk4b.mk4fringe(full_name)
            tot_mbd = mf.t208.contents.tot_mbd  # Total multiband delay, us? 
            tot_sbd = mf.t208.contents.tot_sbd  # Total single-band delay, us?

            if bl in idx.keys():
                if pp in idx[bl].keys():

                    #
                    # Insert time tag, full_name, mbdelay, sbdelay, and snr
                    # into the lists so that to keep ascending time order
                    #

                    if 'time' in idx[bl][pp].keys(): # Just one of the keys
                        #
                        # Find index insr into the time list using fast 
                        # dichotomy (or bisection) algorithm.
                        # The insr index points at the location to insert the
                        # time tag keeping time ascending order.
                        #
                        insr = bisect_right(idx[bl][pp]['time'], ttag)

                        idx[bl][pp]['time'].insert(insr, ttag)
                        idx[bl][pp]['source'].insert(insr, src)
                        idx[bl][pp]['dir'].insert(insr, dir_name)
                        idx[bl][pp]['file'].insert(insr, filename)
                        idx[bl][pp]['full_fname'].insert(insr, full_name)
                        idx[bl][pp]['mbdelay'].insert(insr, mbdelay)
                        idx[bl][pp]['sbdelay'].insert(insr, sbdelay)
                        idx[bl][pp]['snr'].insert(insr, snr)
                        idx[bl][pp]['tot_mbd'].insert(insr, tot_mbd)
                        idx[bl][pp]['tot_sbd'].insert(insr, tot_sbd)
                        idx[bl][pp]['phase'].insert(insr, phase)
                        idx[bl][pp]['time_tag'].insert(insr, ttag)

                    else:
                        idx[bl][pp] = {'time':[ttag], 'source':[src],
                                       'dir': [dir_name],
                                       'file': [filename],
                                       'full_fname':[full_name], \
                                       'mbdelay': [mbdelay], \
                                       'sbdelay': [sbdelay], \
                                       'snr': [snr], \
                                       'tot_mbd': [tot_mbd], \
                                       'tot_sbd': [tot_sbd], \
                                       'phase': [phase], 'time_tag': [ttag]}

                else: # Polproduct subdictionary does not exist in the baseline
                      # subdictionary yet. Create it.
                      # New dict {time,name,mbdelay,sbdelay,snr} 
                      # for polproduct pp
                    idx[bl][pp] = {'time':[ttag], 'source':[src],
                                   'dir': [dir_name], 'file': [filename],
                                   'full_fname':[full_name], \
                                   'mbdelay': [mbdelay], 'sbdelay': [sbdelay], 
                                   'snr': [snr], \
                                   'tot_mbd': [tot_mbd], \
                                   'tot_sbd': [tot_sbd], \
                                   'phase': [phase],
                                   'time_tag': [ttag]}
            else: # Baseline subdictionary does not exist in the idx
                  # dictionary yet. Create new baseline subdictionary with 
                  # a new polproduct subdictionary inside.
                idx[bl] = {}                      # New dict for baseline
                # New dict {time,name,mbdelay,sbdelay,snr} for polproduct pp
                idx[bl][pp] = {'time':[ttag], 'source':[src],
                               'dir': [dir_name], 'file': [filename],
                               'full_fname':[full_name], \
                               'mbdelay': [mbdelay], 'sbdelay': [sbdelay], 
                               'snr': [snr], \
                               'tot_mbd': [tot_mbd], 'tot_sbd': [tot_sbd], \
                               'phase': [phase], 'time_tag': [ttag]}

# ========================= idxs[sr][tm][bl][data_name]========================


            if src in idxs1.keys():
                if ttag in idxs1[src].keys():
                    idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                           'sbdelay': sbdelay,
                                           'tot_mbd': tot_mbd,
                                           'tot_sbd': tot_sbd, 
                                           'phase': phase,
                                           'snr': snr,
                                           'pol_prod': pp,
                                           'dir': dir_name,
                                           'file': filename,
                                           'full_fname': full_name,
                                           'time_tag': ttag}
                else: # ttag subdictionary does not exist in the idxs1[src]
                      # subdictionary yet. Create new time subdictionary with 
                      # a new baseline subdictionary inside.
                    idxs1[src][ttag] = {}           # New dict for time
                    idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                           'sbdelay': sbdelay,
                                           'tot_mbd': tot_mbd,
                                           'tot_sbd': tot_sbd,
                                           'phase': phase,
                                           'snr': snr,
                                           'pol_prod': pp,
                                           'dir': dir_name,
                                           'file': filename,
                                           'full_fname': full_name,
                                           'time_tag': ttag}
                    
            else: # Source subdictionary does not exist in the idxs1 dictionary
                  # yet. Create.
                idxs1[src] = {}           # New subdict for source
                idxs1[src][ttag] = {}     # New subdict for time
                idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                       'sbdelay': sbdelay,
                                       'tot_mbd': tot_mbd,
                                       'tot_sbd': tot_sbd,
                                       'phase': phase,
                                       'snr': snr,
                                       'pol_prod': pp,
                                       'dir': dir_name,
                                       'file': filename,
                                       'full_fname': full_name,
                                       'time_tag': ttag}

# ========================= idxf[dir][file][data_name]=======================

            if dir_name in idxf.keys():
                idxf[dir_name][filename] = {'source':src,
                                'time':ttag,
                                'bl': bl,
                                'mbdelay': mbdelay,
                                'sbdelay': sbdelay,
                                'phase': phase,
                                'snr': snr,
                                'pol_prod': pp,
                                'tot_mbd': tot_mbd,
                                'tot_sbd': tot_sbd,
                                'full_fname': full_name,
                                'time_tag': ttag}
            else:
                idxf[dir_name] = {}
                idxf[dir_name][filename] = {'source':src,
                                'time':ttag,
                                'bl': bl,
                                'mbdelay': mbdelay,
                                'sbdelay': sbdelay,
                                'phase': phase,
                                'snr': snr,
                                'pol_prod': pp,
                                'tot_mbd': tot_mbd,
                                'tot_sbd': tot_sbd,
                                'full_fname': full_name,
                                'time_tag': ttag}
        #
        # Show progress printing the data directory name just processed
        #
        data_dir = os.path.basename(os.path.normpath(root_dir))
        # print("%s/ done ..." % data_dir)
        if i_linepad < n_linepad: # Accumulate dirs in dir_line
            dir_line.append(data_dir)
            i_linepad = i_linepad + 1
        else: # Dump dir_line with dirs left-aligned in 11-char field
            fline = ''.join(f'{ds:<11}' for ds in dir_line)
            print(fline, '...')
            dir_line = []
            i_linepad = 0

        # End of loop 'for file in files:'
    # End of loop 'for root_dir, subdirs, files in os.walk(base_dir):'

    if dir_line: # Dump the rest of dirnames (if any)
        fline = ''.join(f'{ds:<11}' for ds in dir_line)
        print(fline, ' ...')


    #
    # Find session time start ttim0 (time of the first scan)
    #
    bls = list(idx.keys())    # All the baselines
    tim_total = []     # List of all the time counts for all the baselines

    for bl in bls:
        timbl = idx[bl]['I']['time_tag']
        tim_total.extend(timbl)

    ttim = np.unique(tim_total)   # Unite all the time counts in one array
    ttim.sort()
    ttim0 = ttim[0]    # Session time start (the fitst scan)

    print("Session time start ttim0 = ", ttim0)

        
    #
    # The dict idxs1[src][time][bl][data_name] has been created with
    # generally unsorted times. The dict idxs is a copy of idxs1, but
    # with times sorted in ascending order.
    #
    # Init idxs with the source keys only
    #
    idxs = {sr : None for sr in idxs1.keys()}

    #
    # Rearrange each source subdictionary in idxs1[sr][tm][bl][data_name]
    # into time ascending order in idxs[sr][tm][bl][data_name],
    # changing the time keys to the times counted from the experiment start
    #
    for sr in idxs1.keys():
        idxs_tm = {tm-ttim0: idxs1[sr][tm] for tm in sorted(idxs1[sr].keys())}
        idxs[sr] = idxs_tm

    #
    # In idx, change 'time' to the session time
    #
    for bl in idx.keys():
        for pp in idx[bl].keys():
            # idx[bl][pp]['time_tag'] = idx[bl][pp]['time']
            idx[bl][pp]['time'] -= ttim0
            
    #
    # In idx, replace all the numeric lists with numpy arrays
    #
    for bl in idx.keys():
        for pp in idx[bl].keys():
            dat = idx[bl][pp]
            dat['mbdelay'] = np.array(dat['mbdelay'])
            dat['sbdelay'] = np.array(dat['sbdelay'])
            dat['tot_mbd'] = np.array(dat['tot_mbd'])
            dat['tot_sbd'] = np.array(dat['tot_sbd'])
            dat['snr'] =     np.array(dat['snr'])
            dat['time'] =    np.array(dat['time'])
            dat['time_tag'] = np.array(dat['time_tag'])
            dat['phase'] =   np.array(dat['phase'])

    #
    # In idxf, add data name 'time_tag', and change 'time' to the session time
    #
    for dr in idxf.keys():
        for fl in idxf[dr].keys():
            # idxf[dr][fl]['time_tag'] = idxf[dr][fl]['time']
            idxf[dr][fl]['time'] -= ttim0

            
    return idx, idxs, idxf















