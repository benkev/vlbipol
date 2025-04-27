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

import os, sys
import re
import pickle
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
            
            ttag = f_obj.time_tag          # Float, time or measurement 
            mbdelay = f_obj.mbdelay        # Float, multiband delay 
            sbdelay = f_obj.sbdelay        # Float, single-band delay 
            snr = f_obj.snr                # Float, signal to noise ratio

            mf = mk4b.mk4fringe(full_name)
            tot_mbd = mf.t208.contents.tot_mbd
            tot_sbd = mf.t208.contents.tot_sbd

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
                                       'phase': [phase]}

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
                                   'phase': [phase]}
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
                               'phase': [phase]}

# ========================= idxs[sr][tm][bl][data_name]========================


            if src in idxs1.keys():
                if ttag in idxs1[src].keys():
                    idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                           'sbdelay': sbdelay,
                                           'snr': snr,
                                           'tot_mbd': tot_mbd,
                                           'tot_sbd': tot_sbd,
                                           'dir': dir_name,
                                           'file': filename,
                                           'full_fname': full_name, 
                                           'phase': phase}
                else: # ttag subdictionary does not exist in the idxs1[src]
                      # subdictionary yet. Create new time subdictionary with 
                      # a new baseline subdictionary inside.
                    idxs1[src][ttag] = {}           # New dict for time
                    idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                           'sbdelay': sbdelay,
                                           'snr': snr,
                                           'tot_mbd': tot_mbd,
                                           'tot_sbd': tot_sbd,
                                           'dir': dir_name,
                                           'file': filename,
                                           'full_fname': full_name,
                                           'phase': phase}
                    
            else: # Source subdictionary does not exist in the idxs1 dictionary
                  # yet. Create.
                idxs1[src] = {}           # New subdict for source
                idxs1[src][ttag] = {}     # New subdict for time
                idxs1[src][ttag][bl] = {'mbdelay': mbdelay,
                                       'sbdelay': sbdelay,
                                       'snr': snr,
                                       'tot_mbd': tot_mbd,
                                       'tot_sbd': tot_sbd,
                                       'dir': dir_name,
                                       'file': filename,
                                       'full_fname': full_name,
                                       'phase': phase}

# ========================= idxf[dir][file][data_name]=======================

            if dir_name in idxf.keys():
                idxf[dir_name][filename] = {'source':src,
                                'time':ttag,
                                'mbdelay': mbdelay,
                                'sbdelay': sbdelay,
                                'snr': snr,
                                'tot_mbd': tot_mbd,
                                'tot_sbd': tot_sbd,
                                'full_fname': full_name,
                                'phase': phase}
            else:
                idxf[dir_name] = {}
                idxf[dir_name][filename] = {'source':src,
                                'time':ttag,
                                'mbdelay': mbdelay,
                                'sbdelay': sbdelay,
                                'snr': snr,
                                'tot_mbd': tot_mbd,
                                'tot_sbd': tot_sbd,
                                'full_fname': full_name,
                                'phase': phase}
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
        timbl = idx[bl]['I']['time']
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

    return idx, idxs, idxf


if __name__ == '__main__':

    #
    # Linear polarization
    #
    
    # linI_2187 = "/media/benkev/Seagate_Backup_Plus_5TB_2/Work/" \
    #             "2187/scratch/Lin_I/2187"

    linI_2187 = "/home/benkev/Work/2187/scratch/Lin_I/2187"
       
    idxl, idxsl, idxfl = make_idx(linI_2187)
    print("Created idxl, idxsl, and idxfl, linear polarization")
    
    with open('idx2187lI.pkl', 'wb') as fout: pickle.dump(idxl, fout)
    with open('idxs2187lI.pkl', 'wb') as fout: pickle.dump(idxsl, fout)
    with open('idxf2187lI.pkl', 'wb') as fout: pickle.dump(idxfl, fout)

    # sys.exit(0)

    
    #
    # Circular polarization
    #

    cirI_2187 = "/home/benkev/Work/vo2187_exprm/DiFX_pconv/2187"

    idxc, idxsc, idxfc = make_idx(cirI_2187, 'cir')
    print("Created idxc, idxsc, and idxfc, linear polarization")
    
    with open('idx2187cI.pkl', 'wb') as fout: pickle.dump(idxc, fout)
    with open('idxs2187cI.pkl', 'wb') as fout: pickle.dump(idxsc, fout)
    with open('idxf2187cI.pkl', 'wb') as fout: pickle.dump(idxfc, fout)

    
    # idx2187cI = make_idx(cirI_2187, 'cir')
    # print("Created idx2187cI, circular polarization")

    # with open('idx2187cI.pkl', 'wb') as fout:
    #     pickle.dump(idx2187cI, fout)            # Pickle the index dict

    # sys.exit(0)


    
    # ======================================================================
    
    #
    # Unpickle it:
    #

    # with open('idx2187cI.pkl', 'rb') as finp: idx2187cI_1 = pickle.load(finp)

