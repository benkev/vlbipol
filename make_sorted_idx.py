help_text = '''
make_idx.py: define function make_idx().

make_idx():  Creates and returns index file to select data files by baselines
             and polarizations. It is a dictionary of baselines, each baseline
             being a subdictionary of polproducts. The polproduct is in turn
             a subdictionary with keys for 'time', 'file', 'ngdelay', and 
             possibly, in future, some more keys extracted grom the data files.
             The value for each key is a list. 
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

import os
import re
import pickle
from bisect import bisect_right  # Bisection algorithm to efficiently search
# import numpy as np
# import matplotlib.pyplot as plt

import hopstestb as ht
import ffcontrol
from vpal import fringe_file_manipulation as ffm
from vpal.utility import int_to_time, time_to_int


def make_idx(base_dir, pol='lin', max_depth=2):
    """
    Returns an index dictionary for all the fringe (type2) files found in 
    any directory under the base_directory.
    pol: 'lin' - linear polarization
         'cir' - circular polarization

    """

    assert os.path.isdir(base_dir)
    
    base_dir = base_dir.rstrip(os.path.sep)    # Remove trailing "/", if any
    # base_dir = os.path.abspath(base_dir)
    num_sep = base_dir.count(os.path.sep)
    # ff_list = []

    #
    # Polarization notation table
    #
    lin2cir = {'XX':'LL', 'XY':'LR', 'YX':'RL', 'YY':'RR'}

    idx = {}

    for root_dir, subdirs, files in os.walk(base_dir):

        # Apply the max depth filter
        if root_dir.count(os.path.sep) > num_sep + max_depth:
            continue

        f_obj = ffm.FringeFileHandle()
        
        for file in files:

            #abs_filename = os.path.abspath(filename)
            #filename_base = os.path.split(abs_filename)[1]
            
            #
            # Only if the file matches the pattern, its name is 
            # assigned to filename. The pattern is: two baseline capital 
            # letters, dot, letter X, dot, one ore more digits, dot,
            # six alphanumeric characters.
            #
            filename = re.findall(r"[A-Z]{2}\.X.[0-9]+.\w{6}", file)
            if filename == []:
                continue

            filename = filename[0]
            full_name = os.path.join(root_dir, filename)

            bl = filename[:2]  # Baseline is first two letters of filename
            
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
            ttag = f_obj.time_tag          # Float, time or measurement 
            mbdelay = f_obj.mbdelay        # Float, multiband delay 
            sbdelay = f_obj.sbdelay        # Float, single-band delay 
            snr = f_obj.snr                # Float, signal to noise ratio 

            if bl in idx.keys():
                if pp in idx[bl].keys():

                    #
                    # Insert time tag, full_name, mbdelay, sbdelay, and snr
                    # into the lists so that to keep ascending time order
                    #

                    if 'time' in idx[bl][pp].keys(): 
                        #
                        # Find index insr into the time list using fast 
                        # dichotomy (or bisection) algorithm.
                        # The insr index points at the location to insert the
                        # time tag keeping time ascending order.
                        #
                        insr = bisect_right(idx[bl][pp]['time'], ttag)

                        idx[bl][pp]['time'].insert(insr, ttag)
                        idx[bl][pp]['file'].insert(insr, full_name)
                        idx[bl][pp]['mbdelay'].insert(insr, mbdelay)
                        idx[bl][pp]['sbdelay'].insert(insr, sbdelay)
                        idx[bl][pp]['snr'].insert(insr, snr)

                    else:
                        idx[bl][pp] = {'time':[ttag], 'file':[full_name], \
                                       'mbdelay': [mbdelay], \
                                       'sbdelay': [sbdelay], 'snr': [snr]}

                else: # Polproduct subdictionary does not exist in the baseline
                      # subdictionary yet. Create it.
                      # New dict {time,name,mbdelay,sbdelay,snr} 
                      # for polproduct pp
                    idx[bl][pp] = {'time':[ttag], 'file':[full_name], \
                                   'mbdelay': [mbdelay], 'sbdelay': [sbdelay], 
                                   'snr': [snr]}
            else: # Baseline subdictionary does not exist in the idx
                  # dictionary yet. Create new baseline subdictionary with 
                  # a new polproduct subdictionary inside.
                idx[bl] = {}                      # New dict for baseline
                # New dict {time,name,mbdelay,sbdelay,snr} for polproduct pp
                idx[bl][pp] = {'time':[ttag], 'file':[full_name], \
                               'mbdelay': [mbdelay], 'sbdelay': [sbdelay], 
                               'snr': [snr]}

    return idx



lin_3819 = "/data-sc16/geodesy/3819/"
cir_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/pol_prods1/3819"
cirI_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/" \
            "pcphase_stokes_test/3819"


idx3819l = make_idx(lin_3819)
print("Created idx3819l, linear polarization")

idx3819c = make_idx(cir_3819, 'cir')
print("Created idx3819c, circular cross-polarization")

idx3819cI = make_idx(cirI_3819)
print("Created idx3819cI, circular polarization, pseudo-Stokes I")

#
# Pickle the index dict
#
with open('idx3819l.pkl', 'wb') as fout:
    pickle.dump(idx3819l, fout)

with open('idx3819c.pkl', 'wb') as fout:
    pickle.dump(idx3819c, fout)

with open('idx3819cI.pkl', 'wb') as fout:
    pickle.dump(idx3819cI, fout)

#
# Unpickle it:
#
with open('idx3819l.pkl', 'rb') as finp:
    idx3819l_1 = pickle.load(finp)

with open('idx3819c.pkl', 'rb') as finp:
    idx3819c_1 = pickle.load(finp)

with open('idx3819cI.pkl', 'rb') as finp:
    idx3819cI_1 = pickle.load(finp)

