help_text = '''
make_idx.py: Create index file to select data files by baselines
             and polarizations. It is a dictionary with access to the
             list of data file. For example, the files for baseline 'EV'
             and polarization 'XY' are
             idx['EV']['XY'].

             The data file lists are sorted in ascending time order.
'''

import os
import re
import pickle
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
            elif pp_list == ['XX', 'YY']: # For circular polarization
                pp = 'I'
            else:
                continue

            f_obj.load(full_name)
            ttag = f_obj.time_tag        # Float 

            if bl in idx.keys():
                if pp in idx[bl].keys():
                    #
                    # HERE: Insert full_name into the list so that to 
                    #       keep ascending time order
                    #
#                    idx[bl][pp].append(full_name) # Add file to the list

                    if 'time' in idx[bl][pp].keys():

                        # Find index into the file list to insert a file
                        llen = len(idx[bl][pp]['time'])
                        for insr in llen:
                            if ttag >= idx[bl][pp]['time'][insr]:
                                break;
                        if insr < llen:
                            idx[bl][pp]['time'].insert(insr, ttag)
                            idx[bl][pp]['file'].insert(insr, full_name)
                        else:
                            idx[bl][pp]['time'].append(ttag)
                            idx[bl][pp]['file'].append(full_name)

                    else:

                        idx[bl][pp] = [{'time':[ttag], 'file':[full_name]}]

?????????????????????????????????????????????????????????????????????????????

                else:
                    # New dict {time,name} for polproduct pp
                    idx[bl][pp] = [{'time':[ttag], 'file':[full_name]}]
            else:
                idx[bl] = {}                      # New dict for baseline
                # New dict {time,name} for polproduct pp
                idx[bl][pp] = [{'time':[ttag], 'file':[full_name]}]

    return idx



lin_3819 = "/data-sc16/geodesy/3819/"
cir_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/pol_prods1/3819"
cirI_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/" \
            "pcphase_stokes_test/3819"

#dlin = lin_3819 + "251-1130/"


idx = make_idx(lin_3819)
# idx = make_idx(cir_3819, 'cir')
# idx = make_idx(cirI_3819)

#
# Pickle the index dict
#
with open('idxs1.pkl', 'wb') as fout:
    pickle.dump(idx, fout)

# with open('idxsc1.pkl', 'wb') as fout:
#     pickle.dump(idx, fout)

# with open('idxscI1.pkl', 'wb') as fout:
#     pickle.dump(idx, fout)

#
# Unpickle it:
#
with open('idxs1.pkl', 'rb') as finp:
    idxs1 = pickle.load(finp)

