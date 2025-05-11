'''
librd.py - functions used to read Mark4 fringe-fit data into
           Python dictionaries.

make_idx(base_dir, pol='lin', max_depth=2):
    Returns index dictionaries for all the fringe-fit files found in 
    any directory under the base_directory at max_depth.


NOTE: The vpal module imported here changes the matplotlib backend
      to "Agg", which is non-interactive and can’t show GUI windows, it is
      meant for saving images in files only. You can check which backend is
      in use:

      import matplotlib
      matplotlib.get_backend()

      If you import librd and want to use matplotlib code, before plotting
      you have to reset the backend to interactive. For example:

      import matplotlib
      matplotlib.use('qtagg', force=True)  # force reset the backend
'''

import os, re
import pickle
import numpy as np
from bisect import bisect_right  # Bisection algorithm for efficient search
import libvp

import hopstestb as ht
from vpal import fringe_file_manipulation as ffm
import mk4b
# from vpal.utility import int_to_time, time_to_int
# import ffcontrol


def make_idx(base_dir, pol='lin', max_depth=2):
    '''
    Returns an index dictionary for all the fringe-fit files found in 
    any directory under the base_directory at max_depth.

    Parameters:

    base_dir: full path to the 4-digit-named directory with the Mark4
         fringe-fit data.
    
    pol: 'lin' - linear polarization, 'cir' - circular polarization
         The circularly-polarized Mark4 fringe-fit data files obtained through
         the conversion chain PolConvert -> difx2mark4 -> fourfit retain
         the 'linear' notations of their polarization products, so pol='cir'
         parameter makes the replacement 
         'XX'->'LL', 'XY'->'LR', 'YX'->'RL', 'YY'->'RR', ['XX', 'YY'] -> 'I'.

    max_depth: makes the algorithm recurse no deeper than max_depth below
               the root directory given in base_dir

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
             idx['EV']['XY']['time'] is the list of times (s) for 'EV' and 'XY'
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

             idxs['1803+784'][4065.0]['IT'] ->
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

             idxf['188-0435a']['HT.X.4.3HJS31'] -->
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

            #
            # Can I use pp = f_obj.pol_product instead      ???
            #
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
            dtec = f_obj.dtec

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
                        idx[bl][pp]['dtec'].insert(insr, dtec)
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
                                       'phase': [phase], \
                                       'dtec': [dtec], \
                                       'time_tag': [ttag]}

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
                                   'dtec': [dtec], \
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
                               'phase': [phase], 'dtec': [dtec], \
                               'time_tag': [ttag]}

# ========================= idxs[sr][tm][bl][data_name]========================


            if src in idxs1.keys():
                if ttag in idxs1[src].keys():
                    idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                           'sbdelay': sbdelay,
                                           'tot_mbd': tot_mbd,
                                           'tot_sbd': tot_sbd, 
                                           'phase': phase,
                                           'dtec': dtec,
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
                                           'dtec': dtec,
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
                                       'dtec': dtec,
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
                                'dtec': dtec,
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
                                'dtec': dtec,
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
    # Find session time start stim0 (time of the first scan)
    #
    # bls = list(idx.keys())    # All the baselines
    # tim_total = []     # List of all the time counts for all the baselines

    # for bl in bls:
    #     timbl = idx[bl]['I']['time_tag']
    #     tim_total.extend(timbl)

    # ttim = np.unique(tim_total)   # Unite all the time counts in one array
    # ttim.sort()
    # stim0 = ttim[0]    # Session time start (the fitst scan)

    stim0 = libvp.session_time_start(idx)

    print("Session time start stim0 = ", stim0)

        
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
        idxs_tm = {tm-stim0: idxs1[sr][tm] for tm in sorted(idxs1[sr].keys())}
        # Add 'time' (seconds) and 'thour' (hours) to data_names
        for tm in idxs_tm.keys():
            for bl in idxs_tm[tm].keys():
                idxs_tm[tm][bl]['time'] = tm
                idxs_tm[tm][bl]['thour'] = tm/3600
        idxs[sr] = idxs_tm

    #
    # In idx, change 'time' to the session time, and add 'thour' key
    #
    for bl in idx.keys():
        for pp in idx[bl].keys():
            tm = idx[bl][pp]['time'] - stim0
            idx[bl][pp]['time'] = tm
            idx[bl][pp]['thour'] = tm/3600
            
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
            dat['phase'] =   np.array(dat['phase'])
            dat['dtec'] =   np.array(dat['dtec'])
            dat['time_tag'] = np.array(dat['time_tag'])
            dat['thour'] =    np.array(dat['thour'])

    #
    # In idxf, add data name 'time_tag', and change 'time' to the session time.
    # Also, add 'thour' key
    #
    for dr in idxf.keys():
        for fl in idxf[dr].keys():
            tm = idxf[dr][fl]['time'] - stim0  # Subtract session start time tag
            idxf[dr][fl]['time'] = tm
            idxf[dr][fl]['thour'] = tm/3600

            
    return idx, idxs, idxf











