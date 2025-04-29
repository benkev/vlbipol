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
plot_1803_784 = False

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
import copy
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib.patches as patches
from itertools import combinations
from bisect import bisect_right  # Bisection algorithm for efficient search
from group_tails import find_tail_bounds, group_tails

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



            
def make_idxst_bl(idx, bls, ttim0):
    '''
    Create dictionary idxst_bl with the time keys in ascending order:

        idxst_bl[source][time] --> list of baselines

    For each celestial source and each time it has a list of the baselines
    pointing at the celestial source at the time. The time is in seconds
    counted from the session start time ttim0. For absolute time (as in the
    fringe-fit files) simply use 0 for ttim0.

    Parameters:
        idx: the index dictionary created by one of the scripts
             make_sorted_idx_<...>.py from the fringe-fit files of a session
        bls: list of the baselines selected
        ttim0: session start time. For absolute time (as in the fringe-fit
               files) simply use 0 for ttim0.

    Returns: idxst_bl
    
    Examples of using idxst_bl:
    
        idxst_bl['2113+293'][23137.0] --> ['GE', 'GS', 'GT', 'SE', 'TE']
        idxst_bl['0529+483'][57456.0] --> ['GI', 'GM', 'GS', 'GT', 'IM',
                                       'IS', 'IT', 'MS', 'MT']
    '''

    idxst_bl1 = {}  # Time unsorted dictionary

    for bl in bls:
        srcs = idx[bl]['I']['source']
        atms =  np.array(idx[bl]['I']['time']) - ttim0
        nt = len(atms)

        for i in range(nt):
            sr = srcs[i]
            tm = atms[i]

            if sr in idxst_bl1.keys():
                if tm in idxst_bl1[sr].keys():
                    idxst_bl1[sr][tm].append(bl)
                else:
                    idxst_bl1[sr][tm] = [bl]
            else:
                idxst_bl1[sr] = {}
                idxst_bl1[sr][tm] = [bl]

    #
    # Init idxst_bl with source keys only
    #
    idxst_bl = {sr : None for sr in idxst_bl1.keys()}
                
    #
    # Rearrange each source subdictionary in idxst_bl1 into time ascending
    # order in idxst_bl
    #
    for sr in idxst_bl1.keys():
        idxst_bl_tm = {tm: idxst_bl1[sr][tm] \
                       for tm in sorted(idxst_bl1[sr].keys())}
        idxst_bl[sr] = idxst_bl_tm

    return idxst_bl




            
def make_idxst_file(idx, bls, ttim0):
    '''
    Create dictionary idxs_bl with the time keys in ascending order:

        idxst_file[source][time] --> list of fringe-fit files 

    For each celestial source and each time it has a list of the fringe-fit
    files with the data for the celestial source at the time. The time is
    in seconds counted from the session start time ttim0. For absolute time
    (as in the fringe-fit files) simply use 0 for ttim0. 

    Parameters:
        idx: the index dictionary created by one of the scripts
             make_sorted_idx_<...>.py from the fringe-fit files of a session
        bls: list of the baselines selected
        ttim0: session start time. For absolute time (as in the fringe-fit
               files) simply use 0 for ttim0.

    Returns: idxst_file
    
    Examples of using idxst_file:

        idxst_ff['1519-273'][32420.0] -->
                             [('HE', './2187/188-0300a/HE.X.1.3HJRT3'),
                              ('HM', './2187/188-0300a/HM.X.3.3HJRT3'),
                              ('ME', './2187/188-0300a/ME.X.2.3HJRT3')]
    
        idxstt_file['2113+293'][23137.0] -->
                             [('GE', './2187/188-0025b/GE.X.10.3HJRD4'),
                              ('GS', './2187/188-0025b/GS.X.7.3HJRD4'),
                              ('GT', './2187/188-0025b/GT.X.6.3HJRD4'),
                              ('SE', './2187/188-0025b/SE.X.2.3HJRD4'),
                              ('TE', './2187/188-0025b/TE.X.8.3HJRD4')]

        idxst_ile['1849+670'][6040.0] -->
                             [('GE', './2187/187-1940/GE.X.6.3HJQJL'),
                              ('GI', './2187/187-1940/GI.X.4.3HJQJL'),
                              ('GM', './2187/187-1940/GM.X.1.3HJQJL'),
                              ('IE', './2187/187-1940/IE.X.5.3HJQJL'),
                              ('IM', './2187/187-1940/IM.X.3.3HJQJL'),
                              ('ME', './2187/187-1940/ME.X.2.3HJQJL')]
   
    '''
    idxst_file1 = {}  # Time unsorted dictionary

    for bl in bls:
        srcs = idx[bl]['I']['source']
        files = idx[bl]['I']['file']
        atms =  np.array(idx[bl]['I']['time']) - ttim0
        nt = len(atms)

        for i in range(nt):
            sr = srcs[i]
            tm = atms[i]
            file = files[i]

            if sr in idxst_file1.keys():
                if tm in idxst_file1[sr].keys():
                    idxst_file1[sr][tm].append((bl, file))
                else:
                    idxst_file1[sr][tm] = [(bl, file)]
            else:
                idxst_file1[sr] = {}
                idxst_file1[sr][tm] = [(bl, file)]

    #
    # Init idxst_bl with source keys only
    #
    idxst_file = {sr : None for sr in idxst_file1.keys()}

    #
    # Rearrange each source subdictionary in idxst_file1 into time ascending
    # order in idxst_file
    #
    for sr in idxst_file1.keys():
        idxst_file_tm = {tm: idxst_file1[sr][tm] \
                       for tm in sorted(idxst_file1[sr].keys())}
        idxst_file[sr] = idxst_file_tm

    

    return idxst_file



            
def make_idxst_tri(idxst_bl):
    '''
    Create dictionary of baseline triplets idxst_tri
    with the time keys in ascending order:

        idxst_tri[source][time][triangle] --> list of baselines in triangle

    For each celestial source, available time and available triangle it has
    a list of the baseline triplets making up the triangle.

    For each celestial source and each time it has a list of the baselines
    pointing at the celestial source at the time. The time is in seconds
    counted from the session start time.

    Parameters:
        idxs: the index dictionary created by one of the scripts
             make_sorted_idx_<...>.py from the fringe-fit files of a session

    Returns: idxst_tri
    
    Examples of using idxst_bl:
    
        idxst_tri['2113+293'][23137.0]['EGS'] --> ('GS', 'SE', 'GE')
        idxst_tri['2113+293'][23137.0]['EGT'] --> ('GT', 'TE', 'GE')
    or
        idxst_tri['0529+483'][57456.0] -->  {'GIM': ('GI', 'IM', 'GM'),
                                           'GIS': ('GI', 'IS', 'GS'),
                                           'GIT': ('GI', 'IT', 'GT'),
                                           'GMS': ('GM', 'MS', 'GS'),
               
                                           'IMS': ('IM', 'MS', 'IS'),
                                           'IMT': ('IM', 'MT', 'IT')}
    '''

    idxst_tri = copy.deepcopy(idxst_bl)

    for sr in idxst_bl.keys():
        for tm in idxst_bl[sr].keys():
            sr_tm_bls = idxst_bl[sr][tm]

            # Find baseline triangles if only 3 or more bls present
            if len(sr_tm_bls) >= 3:
                sr_tm_tris = find_baseline_triangles(sr_tm_bls)
                idxst_tri[sr][tm] = sr_tm_tris
            else:
                del idxst_tri[sr][tm]
                if not idxst_tri[sr]:  # If idxst_tri[sr] == {} (i.e. empty):
                    del idxst_tri[sr]

    return idxst_tri


            
def make_param_dict(idx, parname,  bls, ttim0):
    '''
    Create dictionary par:

        par[source][time][baseline] --> parameter value

    For each celestial source and each time and each available baseline
    it has the value of parname parameter ('mbdelay', 'sbdelay', 'tot_mbd',
    or 'tot_sbd'). The time is in seconds
    counted from the session start time ttim0. For absolute time (as in the
    fringe-fit files) simply use 0 for ttim0.

    Note: the band delay data are converted from microseconds to picoseconds.
          SNR is stored as is. 

    Parameters:
        idx: the index dictionary created by one of the scripts
             make_sorted_idx_<...>.py from the fringe-fit files of a session
        parname: 'mbdelay', 'sbdelay', 'tot_mbd', or 'tot_sbd'
        bls: list of the baselines selected
        ttim0:s ession start time. For absolute time (as in the fringe-fit
               files) simply use 0 for ttim0.
        us2ps: 

    Returns: 
        par: dictionary par[source][time][baseline] with the parname values.
    
    Examples of getting and using par_l:
    
        idxl: dict with linear polprod data
        parname = 'mbdelay'

        par_l = make_param_dict(idxl, parname,  bls, ttim0)
    
        par_l['2113+293'][23137.0] -->      {'GE': -1689.4368454813957,
                     'GS': -2856.121165677905, 'GT': -2838.0430303514004,
                     'SE': 1164.5995546132326, 'TE': 1150.638679973781}

        par_l['2113+293'][23137.0]['GE'] --> -1689.4368454813957 

        par_l['0529+483'][57456.0] -->
                     {'GI': -14689.7342056036,   'GM': 3380.787093192339,
                      'GS': 1807.482447475195,   'GT': 1858.755829744041,
                      'IM': -3876.4160126447678, 'IS': -5924.326367676258,
                      'IT': -5494.215991348028,  'MS': -1583.591802045703,
                      'MT': -1535.640680231154}

        par_l['0529+483'][57456.0]['MT'] --> -1535.640680231154
    '''
    par = {}

    for bl in bls:
        srcs = idx[bl]['I']['source']
        prms = np.array(idx[bl]['I'][parname])
        if parname in ['tot_mbd', 'tot_sbd']:
            prms = prms * 1e6   # Micro- to picoseconds
        atms =  np.array(idx[bl]['I']['time']) - ttim0
        nt = len(atms)

        for i in range(nt):
            sr = srcs[i]    # Source
            tm = atms[i]    # Time
            pr = prms[i]    # Parameter: mbd, sbd, tot_mbd, tot_sbd, or snr

            if sr in par.keys():
                if tm in par[sr].keys():
                    par[sr][tm][bl] = pr
                else:
                    par[sr][tm] = {}
                    par[sr][tm][bl] = pr
            else:
                par[sr] = {}
                par[sr][tm] = {}
                par[sr][tm][bl] = pr

    return par




def make_closure_delay_stt_dict(idxst_tri, par):
    '''
    Create closure delay dictionary tau_stt (source-time-triangle):

        tau_stt[source][time][triangle] --> closure delay value

    It has the structure similar to that of idxst_tri[source][time][triangle]
    dictionary, but the baseline triplets are substituted with the
    closure delays computed for those triplets.

    Returns:
        tau_stt: closure delay dictionary 

    Examples of getting and using tau_stt_l, closure delay for linear polprods:

        par_l: dictionary with linear polprod parameter values

        tau_stt_l = make_closure_delay_stt_dict(idxst_tri, par_l)

        tau_stt_l['2113+293'][23137.0] --> {'EGS': -2.0847655832767487,
                                        'EGT': 2.032495103776455}

        tau_stt_l['2113+293'][23137.0]['EGS'] --> -2.0847655832767487
        tau_stt_l['2113+293'][23137.0]['EGS'] --> 2.032495103776455
    
        tau_stt_l['0529+483'][57456.0] -->      {'GIM': -21946.937311440706,
                 'GIS': -22421.543020755053, 'GIT': -22042.70602669567,
                 'GMS': -10.287156328558922, 'GMT': -13.609416782855988,
                 'IMS': 464.3185529857874,   'IMT': 82.15929847210646}

        tau_stt_l['0529+483'][57456.0]['IMT'] --> 82.15929847210646
        tau_stt_l['0529+483'][57456.0]['GIS'] --> -22421.543020755053
        tau_stt_l['0529+483'][57456.0]['GMT'] --> -13.609416782855988
    
    '''
    
    tau_stt = copy.deepcopy(idxst_tri)

    for sr in idxst_tri.keys():
        for tm in idxst_tri[sr].keys():
            for tri in idxst_tri[sr][tm].keys():
                # del tau_stt[sr][tm][tri]
                ab, bc, ac = idxst_tri[sr][tm][tri]
                pst = par[sr][tm]  # Dictionary {baseline : parameter}
                tau_stt[sr][tm][tri] = pst[ab] + pst[bc] - pst[ac]

                # print(pst[ab], pst[bc], pst[ac])

    return tau_stt




def make_closure_par_stt_dict(idxst_tri, par):
    '''
    Create closure par dictionary par_stt (source-time-triangle):

        par_stt[source][time][triangle] --> closure par value

    It has the structure similar to that of idxst_tri[source][time][triangle]
    dictionary, but the baseline triplets are substituted with the
    closure pars computed for those triplets.

    Returns:
        par_stt: closure par dictionary 

    Examples of getting and using par_stt_l, closure par
    for linear polprods:

        par_l: dictionary with linear polprod parameter values

        par_stt_l = make_closure_par_stt_dict(idxst_tri, par_l)

        par_stt_l['2113+293'][23137.0] --> {'EGS': -2.0847655832767487,
                                        'EGT': 2.032495103776455}

        par_stt_l['2113+293'][23137.0]['EGS'] --> -2.0847655832767487
        par_stt_l['2113+293'][23137.0]['EGS'] --> 2.032495103776455
    
        par_stt_l['0529+483'][57456.0] -->      {'GIM': -21946.937311440706,
                 'GIS': -22421.543020755053, 'GIT': -22042.70602669567,
                 'GMS': -10.287156328558922, 'GMT': -13.609416782855988,
                 'IMS': 464.3185529857874,   'IMT': 82.15929847210646}

        par_stt_l['0529+483'][57456.0]['IMT'] --> 82.15929847210646
        par_stt_l['0529+483'][57456.0]['GIS'] --> -22421.543020755053
        par_stt_l['0529+483'][57456.0]['GMT'] --> -13.609416782855988
    
    '''
    
    par_stt = copy.deepcopy(idxst_tri)

    for sr in idxst_tri.keys():
        for tm in idxst_tri[sr].keys():
            for tri in idxst_tri[sr][tm].keys():
                # del par_stt[sr][tm][tri]
                ab, bc, ac = idxst_tri[sr][tm][tri]
                pst = par[sr][tm]  # Dictionary {baseline : parameter}
                par_stt[sr][tm][tri] = pst[ab] + pst[bc] - pst[ac]

    return par_stt




def make_closure_delay_tri_dict(tau_stt, trians):
    '''
    Create dictionary tau[triangle][] from tau_stt[source][time][triangle].
    The triangle keys have the order if thiangles in the trians list and
    the same as tribl.keys(). 
    Each triangle key points at a subdictionary
    with three keys, 'time', 'source', and 'tau', pointing at lists
    in ascending time order.

    '''
    
    tau = {}
    for tr in trians:  # Fill in the keys in trians order
        tau[tr] = ()

    for sr in tau_stt.keys():
        for tm in tau_stt[sr].keys():
            for tr in tau_stt[sr][tm].keys():

                clod = tau_stt[sr][tm][tr]  # Closure delay value

                # if tr in tau.keys():

                if 'time' in tau[tr]: # Just one of the keys
                    #
                    # Find index insr into the time list using fast 
                    # dichotomy (or bisection) algorithm.
                    # The insr index points at the location to insert the
                    # time value keeping time ascending order.
                    #
                    insr = bisect_right(tau[tr]['time'], tm)

                    tau[tr]['time'].insert(insr, tm)
                    tau[tr]['source'].insert(insr, sr)
                    tau[tr]['tau'].insert(insr, clod)
                else:
                    tau[tr] = {'time':[tm], 'source':[sr], 'tau':[clod]}

                # else:
                #    tau[tr] = {'time':[tm], 'source':[sr], 'tau':[clod]}

    return tau




def make_closure_phase_tri_dict(phase_stt, trians):
    '''
    Create dictionary cloph[triangle][] from phase_stt[source][time][triangle].
    The triangle keys have the order if thiangles in the trians list and
    the same as tribl.keys(). 
    Each triangle key points at a subdictionary
    with three keys, 'time', 'source', and 'cloph', pointing at lists
    in ascending time order.

    '''
    
    cloph = {}
    for tr in trians:  # Fill in the keys in trians order
        cloph[tr] = ()

    for sr in phase_stt.keys():
        for tm in phase_stt[sr].keys():
            for tr in phase_stt[sr][tm].keys():

                clop = phase_stt[sr][tm][tr]  # Closure phase value

                # Reduce to [-180 .. +180]
                clop = ((clop + 180) % 360) - 180

                if 'time' in cloph[tr]: # Just one of the keys
                    #
                    # Find index insr into the time list using fast 
                    # dichotomy (or bisection) algorithm.
                    # The insr index points at the location to insert the
                    # time value keeping time ascending order.
                    #
                    insr = bisect_right(cloph[tr]['time'], tm)

                    cloph[tr]['time'].insert(insr, tm)
                    cloph[tr]['source'].insert(insr, sr)
                    cloph[tr]['cloph'].insert(insr, clop)
                else:
                    cloph[tr] = {'time':[tm], 'source':[sr], 'cloph':[clop]}

    return cloph



#
# IDEA: based on make_closure_delay_tri_dict(tau_stt, trians), write
#       a "universal" function make_closure_tri_dict(tau_stt, trians, clopar)
#       creating either closure delay or closure phase dictionary depending on
#       the clopar value, either 'tau' or 'phase'
#


def make_closure_tri_dict(tau_stt, trians, clopar):
    '''
    Create dictionary clos[triangle][] from tau_stt[source][time][triangle].
    The triangle keys have the order if thiangles in the trians list and
    the same as tribl.keys(). 
    Each triangle key points at a subdictionary
    with three keys, 'time', 'source', and 'tau', pointing at lists
    in ascending time order.

    '''
    
    clos = {}
    
    for tr in trians:  # Fill in the keys in trians order
        clos[tr] = ()

    for sr in tau_stt.keys():
        for tm in tau_stt[sr].keys():
            for tr in tau_stt[sr][tm].keys():

                clod = tau_stt[sr][tm][tr]  # Closure delay value

                # if tr in clos.keys():

                if 'time' in clos[tr]: # Just one of the keys
                    #
                    # Find index insr into the time list using fast 
                    # dichotomy (or bisection) algorithm.
                    # The insr index points at the location to insert the
                    # time value keeping time ascending order.
                    #
                    insr = bisect_right(clos[tr]['time'], tm)

                    clos[tr]['time'].insert(insr, tm)
                    clos[tr]['source'].insert(insr, sr)
                    clos[tr]['tau'].insert(insr, clod)
                    clos[tr]['phase'].insert(insr, phase)
                else:
                    clos[tr] = {'time':[tm], 'source':[sr], 'tau':[clod],
                                'phase':[phase]}

                # else:
                #    clos[tr] = {'time':[tm], 'source':[sr], 'tau':[clod]}

    return clos




def make_dic_3var(idxst_tri, par):
    '''
    Create dictionary

      dic_3v[source][time][triangle] -- > array([p1, p2, p3]),

    where p1, p2, and p3 are the values of par parameter at the 3 baselines
    that make up the triangle. For example,

      phase_c = make_param_dict(idxc, 'phase', bls, ttim0)
      dic_3ph_c = make_dic_3var(idxst_tri, phase_c)

    or

      mbd_c = make_param_dict(idxc, 'mbdelay', bls, ttim0)
      dic_3mbd_c = make_dic_3var(idxst_tri, mbd_c)

    Example:

        dic_3ph_l['0529+483'][57456.0] -->
            {'GIM': array([ 314.698517,  118.204346,  185.23468 ]),
             'GIS': array([ 314.698517,  238.002098,  359.85682 ]),
             'GIT': array([ 314.698517,  335.924393,  104.039719]),
             'GMS': array([ 185.23468 ,  191.607224,  359.85682 ]),
             'GMT': array([ 185.23468 ,  302.702412,  104.039719]),
             'IMS': array([ 118.204346,  191.607224,  238.002098]),
             'IMT': array([ 118.204346,  302.702412,  335.924393])}

    
    
    '''
    
    dic_3v = copy.deepcopy(idxst_tri)

    for sr in idxst_tri.keys():
        for tm in idxst_tri[sr].keys():
            for tri in idxst_tri[sr][tm].keys():
                ab, bc, ac = idxst_tri[sr][tm][tri]
                pst = par[sr][tm]  # Dictionary {baseline : parameter}
                dic_3v[sr][tm][tri] = np.array((pst[ab], pst[bc], pst[ac]))

    return dic_3v




def plot_closure_legend(ax_col, trians, cols, par, fs=12):
    '''
    Plot color legend for the closure triangles in a separate axis ax_col.
    For 16 closure triangles, it plots 4 colimns by 4 in height, each showing
    the color and the triangle name on the right.
    Inputs:
        ax_col: axis; better to have shorter in height and wider.
                It is assumed to be on top of other plots and its title
                describes the whole figure.
        trians: list of 3-letter closure triangle names.
        cols:   Array of colors, ntri by 4, from pyplot.cm(numbers 0 to 1)
                For example:
                    cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))
        par:    The parameter name to print in the title in uppercase.
        fs:     Fontsise of the printed triangle names
    '''
    ntri = len(trians)
    qntri = ntri // 4        # Quarter of the number of triangles
    #rntri = ntri % 4         # Residue of triangles
    #cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))           # Colors
    upar = par.upper() # The uppercase parameter name to print in the title.
    
    ax_col.set_xlim(-0.1, 7.5)
    # ax_col.set_ylim(-0.2, qntri+0.5)
    ax_col.set_ylim(-1.2, qntri+0.5)

    #qntri = qntri + 1      # Leave room for the last triangle
    for i in range(qntri):
        y = qntri - i - 1      # Vertical patch position
        k = i                  # Index into trians[k] and cols[k,:]
        rect = patches.Rectangle((0, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(1.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + qntri
        rect = patches.Rectangle((2, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(3.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 2*qntri
        rect = patches.Rectangle((4, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(5.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 3*qntri
        rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(7.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

    # Remaining triangle
    k = 28
    y = -1
    rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(7.1, y+0.2, trians[k], fontsize=fs)

    ax_col.set_axis_off()
    ax_col.set_title("%s Closure" % upar, fontsize=16)


    

def plot_closures_distr(ax_distr, tau, cols, parname, ttl,
                        seltri=None, ms=3, yl=None):
    '''
    Plot distributions of the delay closures for the station triangles
    in tau.keys().

    ms: marker size
    yl: ylim
    '''
    trians = list(tau.keys())
    ntri = len(trians)

    if isinstance(seltri, (list, tuple, np.ndarray)): # Plot selected trians
        for ic in range(ntri):
            tr = trians[ic]
            if tr in seltri:
                timh = np.array(tau[tr]['time'])/3600  # Time in hours
                clod = tau[tr]['tau']  # Closure delays in time ascending order
                ax_distr.plot(timh, clod, '.', color=cols[ic,:], ms=ms)
    else: # Plot all triangles if seltri is not a sequence
        for ic in range(ntri):
            tr = trians[ic]
            timh = np.array(tau[tr]['time'])/3600  # Time in hours
            clod = tau[tr]['tau']  # Closure delays in time ascending order
            ax_distr.plot(timh, clod, '.', color=cols[ic,:], ms=ms)
                    
    ax_distr.grid(1)
    ax_distr.set_title(ttl % (parname.upper(), ntri))
    ax_distr.set_xlabel("hours", fontsize=14)
    ax_distr.set_ylabel("ps", fontsize=14)
    ax_distr.yaxis.set_label_coords(-0.05, 0.58)
    ax_distr.xaxis.set_label_coords(0.55, 0.07)
    if yl:
        ax_distr.set_ylim(yl)


    
def plot_closures_hist_horiz(ax_hist, tau, parname, colr, ttl, nbins=101,
                             seltri=None, xl=None, yl=None):
    '''
    Plot distribution and histogram of a delay closures.
    Parameters:
        ax_hist: axis for plotting
        tau:     closure delay dictionary
        parname:  'mbd', 'sbd', 'tot_mbd', or 'tot_sbd'
        colr:    color of the histogram bars
        ttl:     axis title
        nbins:   number of the histogram bins
        seltri:  list (or tuple, or array) with the selected triangles to plot.
                  If seltri is not a sequence (say, any number), the closure
                  delays for all the station triangles in tau are plotted.
        xl: x limits
        yl: y limits
    
    Returns:
        clod: the sample used to plot the histogram
    
    '''
    #
    # Create zero-element array clod to append closure delays
    #
    clod = np.zeros((0,), dtype=np.float64) 

    if isinstance(seltri, (list, tuple, np.ndarray)): # Gather selected only
        for tr in seltri:
            clod = np.append(clod, tau[tr]['tau'])
    else: # Gather all closures in clod if seltri is not a sequence
        for tr in tau.keys():
            clod = np.append(clod, tau[tr]['tau'])
            
    nclod = len(clod)

    ax_hist.hist(clod, nbins, color=colr,
                     orientation='horizontal')

    
    if xl:
        ax_hist.set_xlim(xl)
    if yl:
        ax_hist.set_ylim(yl)

    ax_hist.grid(1)
    ax_hist.set_title(ttl % (parname.upper(), nclod))
    #ax_hist.xaxis.set_label_coords(0.5, -0.12)
    #ax_hist.set_xticks([]);
    #                ax_hist.tick_params(axis='x', rotation=90)
    #xl = ax_hist.get_xticklabels()
    #ax_hist.set_xticklabels(xl, rotation=-90); 
    #ax_hist.set_yticks([]); 
    ax_hist.set_yticklabels('')

    return nclod
    






























    
# ========================== Body of the Code ================================

# pl.rcParams['text.usetex'] = True # Use LaTeX in Matplotlib text
pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.
#print("pl.isinteractive() -> ", pl.isinteractive())

#
# Unpickle it:
#
with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
with open('idx2187cI.pkl', 'rb') as finp: idxc = pickle.load(finp)

with open('idxs2187lI.pkl', 'rb') as finp: idxsl = pickle.load(finp)
with open('idxs2187cI.pkl', 'rb') as finp: idxsc = pickle.load(finp)

with open('idxf2187lI.pkl', 'rb') as finp: idxfl = pickle.load(finp)
with open('idxf2187cI.pkl', 'rb') as finp: idxfc = pickle.load(finp)

#
# Determine the parameter name 'parname': 'mbdelay', 'sbdelay', or
#     'tmbd' or 'tsbd' or 'rmbd' or 'rsbd'
#
parname = arg_to_par[pararg]

ps = "(ps)"

#
# Circular polarization data uses FEWER baselines than linear pol.
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

#
# Find session time start ttim0 (time of the first scan)
#
stim0 = session_time_start(idxl)
print('Session Time Start = stim0', stim0, '\n')

ttim0 = 0.0

#
# Create array of all the sources
#
srcl = []
srcc = []
for bl in bls:
    srcl.extend(idxl[bl]['I']['source'])
    srcc.extend(idxl[bl]['I']['source'])
rcls = srcl + srcc
asrc = np.unique(rcls)     # Turn into np.array leaving only unique source names
asrc.sort()                # Sort the source names lexicographically

nsrc = len(asrc)    # Number of sources

#
# Create dictionary idxs_ff:
#    idxst_file[source][time] --> list of tuples (baseline, fringe-fit file)
# For each source and each time it has a list of tuples
#     (baseline, fringe-fit file) 
# for the source at the time.
#

# idxst_ff = make_idxst_file(idxl, bls, ttim0)


#
# Create dictionary idxst_bl:
#    idxst_bl[source][time] --> list of baselines
# For each source and each time it  have a list of the baselines pointing
# at the source at the time.
#

idxst_bl = make_idxst_bl(idxl, bls, ttim0)

#
# Create dictionary of baseline triplets idxst_tri:
#    idxst_tri[source][time][triangle]
# For each source, available time and available triangle it has a list of
# the baselines triplets making up the triangle.
#

# idxst_tri replace with idxstt_3bl

idxst_tri = make_idxst_tri(idxst_bl)

#
# Create dictionaries par_l and par_c with parameter values:
#
# Linear polprods:
#     par_l[source][time][baseline] --> mbd, sbd, tot_mbd, tot_sbd, or snr value
#
# Circular polprods:
#     par_c[source][time][baseline] --> mbd, sbd, tot_mbd, tot_sbd, or snr value
#

par_l = make_param_dict(idxl, parname,  bls, ttim0)
par_c = make_param_dict(idxc, parname,  bls, ttim0)

#
# Create dictionaries tau_stt_l and tau_stt_c with closure delay values for the
#    parameter in parname, mbd, sbd, tot_mbd, or tot_sbd.
#
# Linear polprods:
#     tau_stt_l[source][time][triangle] --> closure delay value
#
# Circular polprods:
#     tau_stt_l[source][time][triangle] --> closure delay value
#

tau_stt_l = make_closure_delay_stt_dict(idxst_tri, par_l)
tau_stt_c = make_closure_delay_stt_dict(idxst_tri, par_c)


#
# Create dictionaries tau_l and  tau_c
#
tau_l = make_closure_delay_tri_dict(tau_stt_l, trians)
tau_c = make_closure_delay_tri_dict(tau_stt_c, trians)


phase_l = make_param_dict(idxl, 'phase', bls, ttim0)
phase_c = make_param_dict(idxc, 'phase', bls, ttim0)

cloph_stt_l = make_closure_par_stt_dict(idxst_tri, phase_l)
cloph_stt_c = make_closure_par_stt_dict(idxst_tri, phase_c)

cloph_l = make_closure_phase_tri_dict(cloph_stt_l, trians)
cloph_c = make_closure_phase_tri_dict(cloph_stt_c, trians)



mbd_l = make_param_dict(idxl, 'mbdelay', bls, ttim0)
d3mbd_l = make_dic_3var(idxst_tri, mbd_l)

d3ph_l = make_dic_3var(idxst_tri, phase_l)


def reduce_angle_stt_to_180(cloph_stt):
    for sr in cloph_stt.keys():
        for tm in cloph_stt[sr].keys():
            for tr in cloph_stt[sr][tm].keys():
                clop = cloph_stt[sr][tm][tr]  # Closure phase value
                clop_180 = ((clop + 180) % 360) - 180 # Reduce to [-180 .. +180]
                cloph_stt[sr][tm][tr] = clop_180


reduce_angle_stt_to_180(cloph_stt_l)
reduce_angle_stt_to_180(cloph_stt_c)


#
# ============================= PLOTTING =================================
#

def plot_cloph_stt(cloph_stt, src, tri, col):         # , ttl):
    plt1 = True
    for tm in cloph_stt[src].keys():
        for tr in cloph_stt[src][tm].keys():
            thr = tm/3600          # Time, hours
            if plt1:               # Create label for the first plot only
                if tr == tri:
                    pl.plot(thr, cloph_stt[src][tm][tr], '.', label=src,
                            color=col)
                plt1 = False    
            else:
                if tr == tri:
                    pl.plot(thr, cloph_stt[src][tm][tr], '.', color=col)
                
            # print(thr, cloph_stt[src][tm][tr])
#    pl.title(ttl)



#
# Create dict [triangle]['time', 'cloph] of lists for source sr = '1803+784'
#
src1 = '1803+784'

tric_l = {}
tric_c = {}

for tr in cloph_l.keys():            # Triangles
    nt = len(cloph_l[tr]['time'])    # Number of points for all sources
    for i in range(nt):
        if cloph_l[tr]['source'][i] == src1:
            tm = cloph_l[tr]['time'][i]
            cp_l = cloph_l[tr]['cloph'][i]
            cp_c = cloph_c[tr]['cloph'][i]
            if tr in tric_l.keys():
                tric_l[tr]['time'].append(tm)
                tric_l[tr]['cloph'].append(cp_l)
                tric_c[tr]['time'].append(tm)
                tric_c[tr]['cloph'].append(cp_c)
            else:
                tric_l[tr] = {'time': [tm], 'cloph': [cp_l]}
                tric_c[tr] = {'time': [tm], 'cloph': [cp_c]}

#
# Print numbers of points for each of the triangles (for source src1)
#
# for tr in tric_l.keys():            # Triangles
#     print("%s: %3d points" % (tr, len(tric_l[tr]['time'])))


if plot_1803_784:
#
# Make HUGE number of figures with closure phases of all the triangles
# looking at 1803+784
#

    src_1 = '1803+784'

    for tri_1 in tric_l.keys():
        plt1 = True

        th = np.array(tric_l[tri_1]['time'])/3600   # Time (hours)
        cp_l = np.array(tric_l[tri_1]['cloph'])  # Closure phase, linpol
        cp_c = np.array(tric_c[tri_1]['cloph'])  # Closure phase, cirpol

        npt = len(th)

        pl.figure()

        pl.plot(th, cp_l, 'r-', lw=0.2)
        pl.plot(th, cp_c, 'b-', lw=0.2)

        if plt1:               # Create label for the first plot only
            pl.plot(th, cp_l, 'r.', label="Lin")
            pl.plot(th, cp_c, 'b.', label="Cir")
            plt1 = False
        else:
            pl.plot(th, cp_l, 'r.')
            pl.plot(th, cp_c, 'b.')

        pl.plot([0, 24], [0, 0], 'k', lw=0.5)

        ttl = "VO2187 %s Closure Phase of %s (%d pt)" % \
            (src_1, tri_1, npt)

        pl.title(ttl)
        pl.xlabel("hours")
        pl.ylabel("degrees")
        pl.ylim(-170, 170)
        pl.legend()

        pl.savefig("VO2187_1803+784_Closure_Phase_of_%s_%d_pt.pdf" % \
                    (tri_1, npt), format='pdf')


# sys.exit(0)
        
# ???????????????????
#
# Combine tribl and idxs[src][time][bl][data_name]
# into clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
#
        
clos = {}

for sr in idxsl.keys():
    for tm in idxsl[sr].keys():

        srtm_bls = list(idxsl[sr][tm].keys())

        if len(srtm_bls) < 3: continue  # ===== Ignore less than 3 bls ==== >>

        # Find baseline triangles for 3 or more baselines present
        srtm_tris = find_baseline_triangles(srtm_bls)
        
        ix_srtm = idxsl[sr][tm]

        for tr in srtm_tris.keys():
            #
            # For this source and this time, save triads of parameters involved
            # in calculation of closures`
            #
            bl3 = srtm_tris[tr]     # list ot 3 bls making up a triangle tr
            ph3 = [ix_srtm[bl]['phase'] for bl in bl3]
            mbd3 = [ix_srtm[bl]['mbdelay'] for bl in bl3]
            sbd3 = [ix_srtm[bl]['sbdelay'] for bl in bl3]
            tmbd3 = [ix_srtm[bl]['tot_mbd'] for bl in bl3]
            tsbd3 = [ix_srtm[bl]['tot_sbd'] for bl in bl3]
            snr3 = [ix_srtm[bl]['snr'] for bl in bl3]
            fl3 = [ix_srtm[bl]['file'] for bl in bl3]
            dir3 = [ix_srtm[bl]['dir'] for bl in bl3]
            pp3 = [ix_srtm[bl]['pol_prod'] for bl in bl3]
            pp = pp3[0]

            #
            # Compute all the possible closures
            #
            ab, bc, ac = srtm_tris[tr]   # 3 bls making up a triangle tr
            cloph = ix_srtm[ab]['phase'] + ix_srtm[bc]['phase'] - \
                    ix_srtm[ac]['phase']
            cloph = ((cloph + 180) % 360) - 180   # Reduce to [-180 .. +180]
            tau_mbd = ix_srtm[ab]['mbdelay'] + \
                      ix_srtm[bc]['mbdelay'] - \
                      ix_srtm[ac]['mbdelay']
            tau_sbd = ix_srtm[ab]['sbdelay'] + \
                      ix_srtm[bc]['sbdelay'] - \
                      ix_srtm[ac]['sbdelay']
            tau_tmbd = ix_srtm[ab]['tot_mbd'] + \
                      ix_srtm[bc]['tot_mbd'] - \
                      ix_srtm[ac]['tot_mbd']
            tau_tsbd = ix_srtm[ab]['tot_sbd'] + \
                      ix_srtm[bc]['tot_sbd'] - \
                      ix_srtm[ac]['tot_sbd']
            if sr in clos.keys():
                if tr in clos[sr].keys():
                   #
                    # Find index insr into the time list using fast 
                    # dichotomy (or bisection) algorithm.
                    # The insr index points at the location to insert the
                    # time tag keeping time ascending order.
                    #
                    insr = bisect_right(clos[sr][tr]['time'], tm)

                    clos[sr][tr]['time'].insert(insr, tm)
                    clos[sr][tr]['cloph'].insert(insr, cloph)
                    clos[sr][tr]['tau_mbd'].insert(insr, tau_mbd)
                    clos[sr][tr]['tau_sbd'].insert(insr, tau_sbd)
                    clos[sr][tr]['tau_tmbd'].insert(insr, tau_tmbd)
                    clos[sr][tr]['tau_tsbd'].insert(insr, tau_tsbd)
                    clos[sr][tr]['bl'].insert(insr, bl3)
                    clos[sr][tr]['phase'].insert(insr, ph3)
                    clos[sr][tr]['mbd'].insert(insr, mbd3)
                    clos[sr][tr]['sbd'].insert(insr, sbd3)
                    clos[sr][tr]['tmbd'].insert(insr, tmbd3)
                    clos[sr][tr]['tsbd'].insert(insr, tsbd3)
                    clos[sr][tr]['snr'].insert(insr, snr3)
                    clos[sr][tr]['pol_prod'].insert(insr, pp)
                    clos[sr][tr]['file'].insert(insr, fl3)
                    clos[sr][tr]['dir'].insert(insr, dir3)
                else:
                    clos[sr][tr] = {'time':[tm],
                                    'cloph':[cloph],
                                    'tau_mbd':[tau_mbd],
                                    'tau_sbd':[tau_sbd],
                                    'tau_tmbd':[tau_tmbd],
                                    'tau_tsbd':[tau_tsbd],
                                    'bl':[bl3],
                                    'phase':[ph3],
                                    'mbd':[mbd3],
                                    'sbd':[sbd3],
                                    'tmbd':[tmbd3],
                                    'tsbd':[tsbd3],
                                    'snr':[snr3],
                                    'pol_prod':[pp],
                                    'file':[fl3],
                                    'dir':[dir3]
                                    }
            else:
                clos[sr] = {}
                clos[sr][tr] = {'time':[tm],
                                'cloph':[cloph],
                                'tau_mbd':[tau_mbd],
                                'tau_sbd':[tau_sbd],
                                'tau_tmbd':[tau_tmbd],
                                'tau_tsbd':[tau_tsbd],
                                'bl':[bl3],
                                'phase':[ph3],
                                'mbd':[mbd3],
                                'sbd':[sbd3],
                                'tmbd':[tmbd3],
                                'tsbd':[tsbd3],
                                'snr':[snr3],
                                'pol_prod':[pp],
                                'file':[fl3],
                                'dir':[dir3]
                                }



def make_closure_dic(idxs, bls=None):
    '''
    Create dictionary of all possible closures
        clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
    from the dictionary
        idxs[src][time][bl][data_name]

    The bls parameter is a list of allowed baselines. If not given,
    all the baselines are involved. In the VO2187 experiment, for example,
    the baseline ST and all the baselines with station Y are excluded.

    In the code, the following variables are for brevity:
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









#src_1 = '1803+784'
src_1 = '0059+581'

pl.figure()

plt1 = True
for tm in cloph_stt_l[src_1].keys():
    for tr in cloph_stt_l[src_1][tm].keys():
        thr = tm/3600
        if plt1:
            pl.plot(thr, cloph_stt_l[src_1][tm][tr], 'r.', label=src_1)
        else:
            pl.plot(thr, cloph_stt_l[src_1][tm][tr], 'r.')
        plt1 = False
#         print(thr, cloph_stt_l[src_1][tm][tr])
pl.legend()


# plt1 = True
# for tm in cloph_stt_l[src_2].keys():
#     for tr in cloph_stt_l[src_2][tm].keys():
#         thr = tm/3600
#         if plt1:
#             pl.plot(thr, cloph_stt_l[src_2][tm][tr], 'b.', label='1803+784')
#         else:
#             pl.plot(thr, cloph_stt_l[src_2][tm][tr], 'b.')
#         plt1 = False
#         print(thr, cloph_stt_l[src_2][tm][tr])
# pl.legend()

        
# pl.show()


closl = make_closure_dic(idxsl, bls)
closc = make_closure_dic(idxsc, bls)

dat_1 = 'tau_mbd'

pl.figure()

plt1 = True
for tr in closl[src_1].keys():
    th = closl[src_1][tr]['time']/3600
    if plt1:
        pl.plot(th, closl[src_1][tr][dat_1], 'r.', label=src_1)
        plt1 = False
    else:
        pl.plot(th, closl[src_1][tr][dat_1], 'r.')
    
pl.legend()


pl.figure()

plt1 = True
for tr in closc[src_1].keys():
    th = closc[src_1][tr]['time']/3600
    if plt1:
        pl.plot(th, closc[src_1][tr][dat_1], 'r.', label=src_1)
        plt1 = False
    else:
        pl.plot(th, closc[src_1][tr][dat_1], 'r.')
    
pl.legend()


pl.show()


# for i in range(len(tau_l['EGH']['time'])):
#     print("%7.1f '%8s' %f" % (tau_l[tr]['time'][i],
#                               tau_l[tr]['source'][i],
#                               tau_l[tr]['tau'][i]))

# t = np.array(tau_l['EGH']['time'])
# pl.plot(t)

#
# Colors for the station triangles in the trians list and in
#     tau_l.keys() and tau_c.keys()
#
cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))


upar = parname.upper()

if pararg == "mbd":
    yl = 200
elif pararg == "sbd":  
    yl = 2000
    
dist_ylim = (-yl,yl)
#hist_xlim = (-50,6000)   # SBD


gs_kw1 = dict(width_ratios=[0.75, 0.25], height_ratios=[0.15, 0.425, 0.425])
fig1, axd1 = pl.subplot_mosaic([['col_legend', 'col_legend'],
                               ['distr_frfit', 'hist_frfit'],
                               ['distr_pconv', 'hist_pconv']],
                               gridspec_kw=gs_kw1, figsize=(8.4, 10),
                               layout="constrained")
#
# Plot color legend on top
#
ax_col = axd1['col_legend']    # Get the axis for color legend

plot_closure_legend(ax_col, trians, cols, upar, fs=9)  

hist_colr = 'red'

ax_ffd = axd1['distr_frfit'] # Plot distr of FourFit ps.-I param vs Time  
ttl_ffd = "Fourfit Pseudo-I, %s vs Time (%d triangles)" # % (upar, ntri)

plot_closures_distr(ax_ffd, tau_l, cols, parname, ttl_ffd, yl=dist_ylim)


ax_ffh = axd1['hist_frfit']  # Plot hist of FourFit I param
ttl_ffh = "%s (%d points)" # % (upar, nfinite)

# nclod_l = plot_closures_hist_horiz(ax_ffh, tau_l, parname, hist_colr, ttl_ffh,
#                                  yl=dist_ylim, xl=hist_xlim)
# nclod_l = plot_closures_hist_horiz(ax_ffh, tau_l, parname, hist_colr, ttl_ffh)


clod0 = np.zeros((0,), dtype=np.float64)
for tr in tau_l.keys():
    clod0 = np.append(clod0, tau_l[tr]['tau'])
    
ixc = np.where(abs(clod0) < yl)
clod = np.copy(clod0[ixc])
ax_ffh.hist(clod, 101, color=hist_colr, orientation='horizontal')
ax_ffh.set_ylim(dist_ylim)
# ax_ffh.set_ylim(hist_xlim)


ax_pcd = axd1['distr_pconv'] # Plot distr of PolConvert I param vs Time
ttl_pcd = "PolConvert I, %s vs Time (%d triangles)" # % (upar, ntri)

plot_closures_distr(ax_pcd, tau_c, cols, parname, ttl_pcd, yl=dist_ylim)

ax_pch = axd1['hist_pconv']  # Plot hist of PolConvert I param
ttl_pch = "%s (%d points)" # % (upar, nfinite)

# nclod_c = plot_closures_hist_horiz(ax_pch, tau_c, parname, hist_colr, ttl_pch,
#                          yl=dist_ylim, xl=hist_xlim)
#nclod_c = plot_closures_hist_horiz(ax_pch, tau_c, parname, hist_colr, ttl_pch)

clod0 = np.zeros((0,), dtype=np.float64)
for tr in tau_c.keys():
    clod0 = np.append(clod0, tau_c[tr]['tau'])
    
ixc = np.where(abs(clod0) < yl)
clod = np.copy(clod0[ixc])
ax_pch.hist(clod, 101, color=hist_colr, orientation='horizontal')
ax_pch.set_ylim(dist_ylim)
# ax_pch.set_ylim(hist_xlim)



























pl.show()





sys.exit(0)





# #
# # Find sources with maximal time counts
# #
# tc = []     # Time counts
# stc = []    # Sources for the time counts
# for sr in tau_stt_l.keys():
#     stc.append(sr)
#     tc.append(len(tau_stt_l[sr].keys()))

# tc = np.array(tc)
# itc = np.argsort(-tc)  # Find indices of the time counts in descending order
# tc = tc[itc]           # Sort time counts in descending order
# stc = np.array(stc)
# stc = stc[itc]         # Sort sources in descending order of their time counts

# ns = 58  # Number of sources to be plotted
# cols = cm.nipy_spectral(1 - np.linspace(0, 1, ns)) # Colors for each source
# upar = parname.upper()

# #
# # Plot closure delay
# #
# # fig_cols = pl.figure(figsize=(10,4)); ax = pl.subplot()
# # fig_cols.tight_layout(rect=(0,0,1, 0.95))
# # plot_closure_legend(ax, stc[:ns], cols, parname, fs=10)

# f1 = pl.figure()
# pl.plot([0, 24], [0, 0], 'k')

# atau_stt_l = []
# # trs_l = set()         # Triangles involved
# ic = 0                  # Color index
# for sr in stc[:ns]:
#     for tm in tau_stt_l[sr].keys():
#         # trs_l.update(tau_stt_l[sr][tm].keys())
#         for tr in tau_stt_l[sr][tm].keys():
#             atau_stt_l.append(tau_stt_l[sr][tm][tr])
#             pl.plot(tm/3600, tau_stt_l[sr][tm][tr], '.', color=cols[ic,:], ms=3)
#     ic = ic + 1        

# atau_stt_l = np.array(atau_stt_l)

# pl.ylim(-1200, 1200)
# pl.title("VO2187_Closure Delay (Linear PolProds)");
# pl.xlabel("hours")

# pl.savefig("VO2187_%s_Closure_Delay_Lin.pdf" % upar, format='pdf')

# nbin_ini = 101     # Initial number of histogram bins (before tail grouping)
# ni_l, bedges = np.histogram(atau_stt_l, nbin_ini)
# # Compute bin centers and bin width
# xi_l = (bedges[:-1] + bedges[1:]) / 2
# bw_l = bedges[1] - bedges[0]

# l_idx_l, r_idx_l = find_tail_bounds(ni_l, thr=10)
# lr_l = l_idx_l, r_idx_l
# ni_grp_l =  group_tails(ni_l, lr_l)
# xi_grp_l = xi_l[l_idx_l : r_idx_l]

# f2 = pl.figure()
# pl.bar(xi_grp_l, ni_grp_l, width=bw_l, color='blue', ec='w', align='center')
# #pl.bar(xi_l, ni_l, width=bw_l, color='brown', ec='w', align='center')
# pl.grid(1)
# # pl.hist(atau_stt_l, 101); pl.grid(1)
# pl.xlim(-1200, 1200)
# pl.title("VO2187 Distribution of Closure Delay (Linear PolProds)");
# pl.xlabel("ps")

# pl.savefig("VO2187_%s_Distr_of_Closure_Delay_Lin.pdf" % upar, format='pdf')



# f3 = pl.figure()
# pl.plot([0, 24], [0, 0], 'k', lw=0.4)

# atau_stt_c = []
# # trs_c = set()         # Triangles involved
# ic = 0                  # Color index
# for sr in stc[:ns]:
#     for tm in tau_stt_c[sr].keys():
#         # trs_c.update(tau_stt_c[sr][tm].keys())
#         for tr in tau_stt_c[sr][tm].keys():
#             atau_stt_c.append(tau_stt_c[sr][tm][tr])
#             pl.plot(tm/3600, tau_stt_c[sr][tm][tr], '.', color=cols[ic,:], ms=3)
#     ic = ic + 1        

# atau_stt_c = np.array(atau_stt_c)

    
# pl.ylim(-1200, 1200)
# pl.title("VO2187_Closure Delay (Circular PolProds)");
# pl.xlabel("hours")

# pl.savefig("VO2187_%s_Closure_Delay_Cir.pdf" % upar, format='pdf')

# ni_c, bedges = np.histogram(atau_stt_c, nbin_ini)
# # Compute bin centers and bin width
# xi_c = (bedges[:-1] + bedges[1:]) / 2
# bw_c = bedges[1] - bedges[0]

# l_idx_c, r_idx_c = find_tail_bounds(ni_c, thr=10)
# lr_c = l_idx_c, r_idx_c
# ni_grp_c =  group_tails(ni_c, lr_c)
# xi_grp_c = xi_c[l_idx_c : r_idx_c]


# f4 = pl.figure()
# pl.bar(xi_grp_c, ni_grp_c, width=bw_c, color='blue', ec='w', align='center')
# #pl.bar(xi_c, ni_c, width=bw_c, color='brown', ec='w', align='center')
# pl.grid(1)
# # pl.hist(atau_stt_c, 101); pl.grid(1)
# pl.xlim(-1200, 1200)
# pl.title("VO2187 Distribution of Closure Delay (Circular PolProds)");
# pl.xlabel("ps")

# pl.savefig("VO2187_%s_Distr_of_Closure_Delay_Cir.pdf" % upar, format='pdf')







# # def plot_closures_hist(ax_hist, timh, atau, seltri, pararg, ttl):
# #     '''
# #     Plot distribution and histogram of a delay closures.
# #     '''

# #     if isinstance(seltri, (list, tuple, np.ndarray)):
# #         sel = np.array(seltri, dtype='bool')  # Any seq. into bool array sel
# #         # print("atau[sel,:].flatten().shape = ", atau[sel,:].flatten().shape)
# #         ax_hist.hist(abs(atau[sel,:].flatten()), 50)
# #     else: # if seltri is elemental, e.g. any number:
# #         ax_hist.hist(abs(atau.flatten()), 50)

# #     ax_hist.grid(1)
# #     ax_hist.set_xlabel("ps", fontsize=14)
# #     ax_hist.set_title(ttl)
# #     ax_hist.xaxis.set_label_coords(0.5, -0.12)


# # #
# # # Create dictionary idx_asrc with the source indices
# # # into the source array asrc
# # #
# # idx_asrc = {}  # Index into the source array asrc

# # for i in range(nsrc):
# #     s = asrc[i]
# #     idx_asrc[s] = i

# # #
# # # Test printouts 
# # #

# # b1 = idxst_bl['2113+293'][23137.0]
# # print(b1)
# # #             ['GE', 'GS', 'GT', 'SE', 'TE']
# # b1tri = find_baseline_triangles(b1)
# # print(b1tri)
# # # {'EGS': ('GS', 'SE', 'GE'), 'EGT': ('GT', 'TE', 'GE')}



# # b2 = idxst_bl['0529+483'][25087.0]
# # print(b2)
# # #             ['IM', 'IS', 'IT', 'MS', 'MT']
# # b2tri = find_baseline_triangles(b2)
# # print(b2tri)
# # # {'IMS': ('IM', 'MS', 'IS'), 'IMT': ('IM', 'MT', 'IT')}

# # b3 = idxst_bl['0529+483'][57456.0]
# # print(b3)
# # #            ['GI', 'GM', 'GS', 'GT', 'IM', 'IS', 'IT', 'MS', 'MT']
# # b3tri = find_baseline_triangles(b3)
# # print(b3tri)
# # # {'GIM': ('GI', 'IM', 'GM'),
# # #  'GIS': ('GI', 'IS', 'GS'),
# # #  'GIT': ('GI', 'IT', 'GT'),
# # #  'GMS': ('GM', 'MS', 'GS'),
# # #  'GMT': ('GM', 'MT', 'GT'),
# # #  'IMS': ('IM', 'MS', 'IS'),
# # #  'IMT': ('IM', 'MT', 'IT')}

# # print("\nidxst_bl dictionary:")
# # for sr in idxst_bl.keys():
# #     print("\nSource '%s':" % sr)
# #     for tm in idxst_bl[sr].keys():
# #         sr_tm_bls = idxst_bl[sr][tm]
# #         if len(sr_tm_bls) >= 3:
# #             print("    t = %.2f (%d bls): " % (tm, len(sr_tm_bls)),
# #                   idxst_bl[sr][tm])



# # print("\nidxst_bl dictionary: times with <3 baselines")
# # for sr in idxst_bl.keys():
# #     print("\nSource '%s':" % sr)
# #     for tm in idxst_bl[sr].keys():
# #         sr_tm_bls = idxst_bl[sr][tm]
# #         if len(sr_tm_bls) < 3:
# #             print("    t = %.2f: " % tm, idxst_bl[sr][tm])

# # print("\nidxst_tri dictionary:")
# # for sr in idxst_tri.keys():
# #     print("\nSource '%s':" % sr)
# #     for tm in idxst_tri[sr].keys():
# #         print("    t = %.1f: " % tm)
# #         for tri in idxst_tri[sr][tm].keys():
# #             sr_tm_tri = idxst_tri[sr][tm][tri]
# #             print("        '%s': " % tri, sr_tm_tri)



# #
# #
# # Plot closure phase for sources
# #
# src_1 = '1803+784'
# src_2 = '0059+581'


# pl.figure()
# ttl = "VO2187 Closure Phase, Lin PolProd"

# # tri_1 = 'EGS'
# tri_1 = 'HIT'
# plot_cloph_stt(cloph_stt_l, src_1, tri_1, 'r')
# # plot_cloph_stt(cloph_stt_l, src_2, tri_1, 'b')

# pl.title(ttl)
# pl.ylim(-170, 170)
# pl.legend()

# pl.savefig("VO2187_1803+784_Closure_Phase_Lin_Pol.pdf", format='pdf')

# pl.figure()
# ttl = "VO2187 Closure Phase, Cir PolProd"

# plot_cloph_stt(cloph_stt_c, src_1, 'r')
# plot_cloph_stt(cloph_stt_c, src_2, 'b')

# pl.title(ttl)
# pl.ylim(-200, 200)
# pl.legend()

# pl.savefig("VO2187_Closure_Phase_Cir_Pol.pdf", format='pdf')

# #
# # Count occurrencies of each triangle to find which one has 30 points
# #
# tric = {}
# for tm in cloph_stt_l['1803+784'].keys():
#     print(tm)
#     trs = list(cloph_stt_l['1803+784'][tm].keys())
#     print(trs)




# tim1 = {}     # Original time points with some of them missing. 
# tim = {}      # Time points. The gaps will be replaced with NaNs

# tim_total = []     # List of all the time counts for all the baselines

# for bl in bls:
#     timbl = idxl[bl]['I']['time']
#     tim_total.extend(timbl)
#     tim1[bl] = np.array(timbl)  #/ 3600 # Sec -> hours

# ttim = np.unique(tim_total)   # Unite all the time counts in one array
# ttim0 = ttim[0]    # 'global' time start

# # Set time for each baseline start at 'global' zero
# for bl in bls:
#     tim1[bl] = tim1[bl] - ttim0  # Set time start at 'global' zero
#     print("tim1[%s] = %.2f" % (bl, tim1[bl][0]))



# # for bl in bls:
# #     tim1[bl] = np.array(idxl[bl]['I']['time'])  #/ 3600 # Sec -> hours
# #     ??????????????????/
# #     tim1[bl] = tim1[bl] - tim1[bl][0]  # Set time start at zero

# #
# # Find the minimum time between the scans over all the baselines
# #
# min_t_scan = 1000000000

# for bl in bls:
#     t_scans = np.diff(tim1[bl])
#     bl_min_t_scan = np.min(t_scans)
#     if bl_min_t_scan < min_t_scan:
#         min_t_scan = bl_min_t_scan

# #
# # Search for the maximum number of time counts, ntim,
# # among all the baselines. At the end of the loop below, ntim will hold
# # the length for the parameter storage in arrays.
# # The parameter sets shorter than ntim contain time gaps to be filled with NaNs.
# #
# ntim = -1       # Will contain the maximum time count
# for bl in bls:
#     bl_ntim = np.max(tim1[bl]/min_t_scan)    # Max t counts for the baseline
#     #print("bl_ntim = %f" % bl_ntim)
#     if bl_ntim > ntim:
#         ntim = np.int64(np.ceil(bl_ntim))

#     print("len(tim1['%s']) = %d; Max t counts = %f" %
#                                    (bl, len(tim1[bl]), bl_ntim))

# # ?? ntim = ntim + 1 # !!!!!!! I DO NOT KNOW WHY 3819 NEEDS IT ???????????????
    
# print("Max time counts: %d;  min scan time: %d s." % (ntim, min_t_scan))




# # ???????????????????
# #
# # Combine tribl and idxs[src][time][bl][data_name]
# # into clos[src][tri]['time', 'cloph', 'tau_mbd', 'tau_sbd' etc.]
# #
        
# clos = {}

# for sr in idxsl.keys():
#     for tm in idxsl[sr].keys():

#         srtm_bls = list(idxsl[sr][tm].keys())

#         if len(srtm_bls) < 3: continue  # ===== Ignore less than 3 bls ==== >>

#         # Find baseline triangles for 3 or more baselines present

#         srtm_tris = find_baseline_triangles(srtm_bls)
#         ix_srtm = idxsl[sr][tm]

#         #
#         # For this source and this time, create dict
#         #   tri_prm[tr][triads of phase, mbd, etc.]
#         # to save triads of parameters involved in calculation of closures
#         # for each of the station triangles in srtm_tris
#         #
#         tri_prm = {}      # tri_prm[tr][triads of phase, mbd, etc.]

#         for tr in srtm_tris.keys():
#             bl3 = srtm_tris[tr]     # list ot 3 bls making up a triangle tr
#             tri_prm[tr] = {'phase': [ix_srtm[bl]['phase'] for bl in bl3],
#                            'mbd': [ix_srtm[bl]['mbdelay'] for bl in bl3],
#                            'sbd': [ix_srtm[bl]['sbdelay'] for bl in bl3],
#                            'tmbd': [ix_srtm[bl]['tot_mbd'] for bl in bl3],
#                            'tsbd': [ix_srtm[bl]['tot_sbd'] for bl in bl3],
#                            'snr': [ix_srtm[bl]['snr'] for bl in bl3],
#                            'file': [ix_srtm[bl]['file'] for bl in bl3],
#                            'dir': [ix_srtm[bl]['dir'] for bl in bl3],
#                            'bl': srtm_tris[tr]}

#         if sr in clos.keys():
#             if tr in clos[sr].keys():

#                 #
#                 # Compute all the possible closures
#                 #
                
#                 #
#                 # Find index insr into the time list using fast 
#                 # dichotomy (or bisection) algorithm.
#                 # The insr index points at the location to insert the
#                 # time tag keeping time ascending order.
#                 #
#                 insr = bisect_right(clos[sr][tr]['time'], tm)

#                 clos[sr][tr]['time'].insert(insr, tm)
#                 clos[sr][tr]['cloph'].insert(insr, cloph)
#                 clos[sr][tr]['tau_mbd'].insert(insr, tau_mbd)
#                 clos[sr][tr]['tau_sbd'].insert(insr, tau_sbd)
#                 clos[sr][tr]['tau_tmbd'].insert(insr, tau_tmbd)
#                 clos[sr][tr]['tau_tsbd'].insert(insr, tau_tsbd)
#                 clos[sr][tr]['snr'].insert(insr, tri_prm[tr]['snr'])
#                 clos[sr][tr]['bl'].insert(insr, srtm_tris[tr])
#                 clos[sr][tr]['mbd'].insert(insr, )
#                 clos[sr][tr]['sbd'].insert(insr, )
#                 clos[sr][tr]['tmbd'].insert(insr, )
#                 clos[sr][tr]['tsbd'].insert(insr, )
#                 clos[sr][tr]['phase'].insert(insr, )
#             else:
#                 clos[sr][tr] = {'time':[tm],
#                                 'cloph':[cloph],
#                                 'phase':[]
#                                 '':[]
#                                 '':[]
#                                 '':[]
#                                 '':[]
#                                 '':[]
#                                 '':[]
#                                 '':[]
#                                 '':[]
#                                 }

















