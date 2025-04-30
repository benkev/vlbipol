help_text = '''
make_sorted_idx.py: define function make_idx().

make_idx():  Creates and returns index dictionary to select data files
             by baselines and polarizations. It is a dictionary of baselines,
             each baseline being a subdictionary of polproducts. The polproduct
             is in turn a subdictionary with keys for 'time', 'file',
             'ngdelay', and possibly, in future, some more keys extracted from
             the data files.
             The value for each lowest-level key is a list. 
             The 'time' key points at the list of time tags in ascensing order.
             The 'file' key points at the list of file names corresponding 
                        the time tags in the 'time' list.
             Other keys point to their respective lists also strictly in 
                        ascending time order.

             Access to the lists
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
'''

import os, sys, re
import copy, pickle
from bisect import bisect_right  # Bisection algorithm to efficiently search
import numpy as np
# import matplotlib.pyplot as plt

import hopstestb as ht
import ffcontrol
from vpal import fringe_file_manipulation as ffm
from vpal.utility import int_to_time, time_to_int
import mk4b


def make_idx(base_dir, pol='lin', max_depth=2):
    """
    Returns an index dictionary for all the fringe (type2) files found in 
    any directory under the base_directory.
    pol: 'lin' - linear polarization
         'cir' - circular polarization

    """

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
        print("%s/ done ..." % data_dir)
        


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
    # In idx, rplace all the numeric lists with numpy arrays
    #
    for bl in isx.keys():
        for pp in isx[bl].keys():
            dat = isx[bl][pp]
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






if __name__ == '__main__':

    #
    # Linear polarization
    #
    
    # linI_2187 = "/media/benkev/Seagate_Backup_Plus_5TB_2/Work/" \
    #             "2187/scratch/Lin_I/2187"

    linI_2187 = "/home/benkev/Work/2187/scratch/Lin_I/2187"
       
    idxl, idxsl, idxfl = make_idx(linI_2187)
    print("Created dictionaries idxl, idxsl, and idxfl, linear polarization")
    
    with open('idx2187lI.pkl', 'wb') as fout: pickle.dump(idxl, fout)
    with open('idxs2187lI.pkl', 'wb') as fout: pickle.dump(idxsl, fout)
    with open('idxf2187lI.pkl', 'wb') as fout: pickle.dump(idxfl, fout)

    print("Linear polarization data source:")
    print("   ", linI_2187)
    print()
    print("Linear polarization data dictionaries pickled and saved on disk\n")
    
    # sys.exit(0)

    
    #
    # Circular polarization
    #

    cirI_2187 = "/home/benkev/Work/vo2187_exprm/DiFX_pconv/2187"

    idxc, idxsc, idxfc = make_idx(cirI_2187, 'cir')
    print("Created dictionaries idxc, idxsc, and idxfc, circular polarization")
    
    with open('idx2187cI.pkl', 'wb') as fout: pickle.dump(idxc, fout)
    with open('idxs2187cI.pkl', 'wb') as fout: pickle.dump(idxsc, fout)
    with open('idxf2187cI.pkl', 'wb') as fout: pickle.dump(idxfc, fout)
    
    print("Circular polarization data source:")
    print("   ", cirI_2187)
    print()
    print("Circular polarization data dictionaries pickled and saved on disk\n")

    
    # idx2187cI = make_idx(cirI_2187, 'cir')
    # print("Created idx2187cI, circular polarization")

    # with open('idx2187cI.pkl', 'wb') as fout:
    #     pickle.dump(idx2187cI, fout)            # Pickle the index dict

    # sys.exit(0)

    print("To load from disk:")
    print('''
    import pickle
    
    with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
    with open('idx2187cI.pkl', 'rb') as finp: idxc = pickle.load(finp)

    with open('idxs2187lI.pkl', 'rb') as finp: idxsl = pickle.load(finp)
    with open('idxs2187cI.pkl', 'rb') as finp: idxsc = pickle.load(finp)

    with open('idxf2187lI.pkl', 'rb') as finp: idxfl = pickle.load(finp)
    with open('idxf2187cI.pkl', 'rb') as finp: idxfc = pickle.load(finp)
    ''')
    print()

    # =======================================================================
    # ======= Creation of dictionaries with all the closures possible =======
    # ======= as clos[src][tri][data_item] (closl and closc),         =======
    # ======= where data_items are 'cloph', 'tau_mbd', 'tau_sbd' etc. =======
    # =======================================================================
    #
    #
    # Circular polarization data use FEWER baselines than linear pol.
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
    # Exclude the 'ST' baseline: the S and T stations are too close
    # to each other
    #
    if 'ST' in bls:
        iST = bls.index('ST')
        bls.pop(iST)

    nbls = len(bls)

    print("Found baselines common for linear and circular polarization")
    print("ST baseline excluded\n")
    
    print()
    print('nbls = ', nbls)
    print('bls = ', bls, '\n')
    print()

    closl = make_closure_dic(idxsl, bls)
    closc = make_closure_dic(idxsc, bls)

    print("Created dictionary of data closures, closl, linear polarization")
    print("Created dictionary of data closures, closc, circular polarization\n")
    
    print("Linear pol. data closure dictionary pickled and saved on disk")
    print("Circular pol. data closure dictionary pickled and saved on disk\n")

    print("To load from disk:")
    print('''
    import pickle
    
    with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)
    with open('clos2187cI.pkl', 'rb') as finp: closc = pickle.load(finp)
    ''')
    print("")





    # ======================================================================
    
    #
    # Unpickle it:
    #
    # import pickle
    #
    # with open('idx2187cI.pkl', 'rb') as finp: idx2187cI_1 = pickle.load(finp)

