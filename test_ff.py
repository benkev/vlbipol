import os
import re
import numpy as np
import matplotlib.pyplot as plt

import hopstestb as ht
import ffcontrol
from vpal import fringe_file_manipulation as ffm
# from vpal.processing import gather_fringe_files

################################################################################

def gather_fringe_files(base_directory, control_file, blines, 
                        pol_products=['I'], include_autos=False, 
                        exclude_list=None, max_depth=2):
    """
    Returns a list of all the fringe (type2) files found in any directory 
    up to max_depth under the base_directory
    In a typical VGOS experiment directory, the fringe files are two levels 
    below the base_directory (the first
    level below the base directory is the scan name).

    The exp_dir is expecting an absolute path to a HOPS experiment number 
    directory (four-digit number)
    The control file is an absolute path

    Assumes fringe files have a six-digit extension with three '.' 
    in the filename.
    """
    
    if exclude_list == None:
        exclude_list=['prepass', 'pre_production', 'make_links', 
                      'test', 'bad_eop', 'setup']
    exclude = set(exclude_list)
    assert os.path.isdir(base_directory)

    fringe_file_list = []

    base_dir = base_directory.rstrip(os.path.sep)
    num_sep = base_dir.count(os.path.sep)

    control_file_hash = ffcontrol.get_control_file_hash(control_file)

    counter=0
    for current_root, subdirectories, files in os.walk(base_directory):

        for filename in files:

            # apply the exclude filter
            if any([e in current_root for e in exclude]):
                continue

            # apply the max depth filter
            if current_root.count(os.path.sep) > num_sep+max_depth:
                continue

            abs_filename = os.path.abspath(filename)
            filename_base = os.path.split(abs_filename)[1]

            #
            # Look for root files using some simple checks
            # Check that there are three dots in the filename base
            #
            if filename_base.count('.') == 3: 
                bline = filename_base.split('.')[0]
                #
                # Make sure leading section of file name is 2-char baseline
                #
                if len(bline)==2 and bline in blines: 
                    full_name = os.path.join(current_root, filename)
                    #
                    # Check that this is a cross correlation if autos excluded
                    #
                    if (include_autos is True) or (bline[0] != bline[1]):
                        # Get the file extension (root_id)
                        extension = filename_base.split('.')[3] 
                        # Check that the extension has a length of 6 chars
                        if len(extension) == 6:
                            # counter+=1
                            # ff_cf_hash = \
                            #     ht.get_control_file_hash_from_fringe(
                            #                                full_name)

                            # print("control_file_hash = ", control_file_hash, 
                            #       ", ff_cf_hash = ", ff_cf_hash)

                            # if control_file_hash == ff_cf_hash:
                            #     fringe_file_list.append(full_name)
                            fringe_file_list.append(full_name)

    #
    # Now we have a list of fringe files that used the prescribed control file.
    # Check for the correct polproduct:
    #
    ff_list = []
    
    for ff in fringe_file_list:
        ff_pp_list = ht.get_file_polarization_product_provisional(ff)

        print(ff, ", pp =", ff_pp_list)

        for pp in ff_pp_list:
            if pp in pol_products:
                #
                # Now that we are sure we have the correct polproduct, 
                # create a FringeFileHandle object and append to the list
                #
                f_obj = ffm.FringeFileHandle()
                f_obj.load(ff)
                ff_list.append(f_obj)

    return ff_list


# =========================================================================


def find_fringe_files(base_dir, pol=['I'], max_depth=2):
    """
    Returns a list of all the fringe (type2) files found in any directory
    under the base_directory.
    """

    assert os.path.isdir(base_dir)
    
    base_dir = base_dir.rstrip(os.path.sep)    # Remove trailing "/", if any
    # base_dir = os.path.abspath(base_dir)
    num_sep = base_dir.count(os.path.sep)
    ff_list = []

    for root_dir, subdirs, files in os.walk(base_dir):

        # apply the max depth filter
        if root_dir.count(os.path.sep) > num_sep + max_depth:
            #print("root_dir.count(os.path.sep)=", root_dir.count(os.path.sep))
            continue

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
            
            pp_list = ht.get_file_polarization_product_provisional(full_name)

            print("pp_list = ", pp_list, ", ", full_name)

            # for pp in pp_list:
            #     if pp in pol:
            #         ff_list.append(full_name)
            #         #print("pp_list = ", pp_list, ", ", full_name)

            if pp_list == pol:
                ff_list.append(full_name)

    return ff_list


lin_3819 = "/data-sc16/geodesy/3819/" 
cir_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/pol_prods1/3819"
cirI_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/" \
            "pcphase_stokes_test/3819"

dlin = lin_3819 + "251-1130/"
dpol = cir_3819 + "251-1130/"

cf_lin = lin_3819 + "cf_3819_MESTVY_pstokes"
cf_cir = cir_3819 + "cf_3819_MESTVY_pcphase_dh2"

# Find out if the file is for 'I'
# ht.get_file_polarization_product_provisional(dlin + 'EV.X.29.2GI1S4')


# ff_list = gather_fringe_files("/data-sc16/geodesy/3819/", 
#                     #"/data-sc16/geodesy/3819/cf_3819_MESTVY_pstokes", 
#                     "/data-sc16/geodesy/3819/cf_3819_GEHILMSTVY_mod", 
#                     ['EV', 'ME'])

# ff_list = gather_fringe_files("/data-sc16/geodesy/3801/", 
#                     #"/data-sc16/geodesy/3801/cf_3801_GEHMSVY_pstokes_mod", 
#                     "/data-sc16/geodesy/3801/cf_3793_GEHMSTVL_finalpass_mod", 
#                     ['EV', 'ME'])


# ff_list = gather_fringe_files(cir_3819, 
#                     "/data-sc16/geodesy/3819/cf_3819_MESTVY_pstokes", 
#                     #"/data-sc16/geodesy/3819/cf_3819_GEHILMSTVY_mod", 
#                     ['EV', 'ME', 'MS', 'MT', 'MV', 'SE', 'SV', 'TE', 'TV'],
#                               ['XY'])

# ff_list = gather_fringe_files(cirI_3819, 
#       "/data-sc16/geodesy/3819/polconvert/3819/cf_3819_MESTVY_pcphase_dh2", 
#                     #"/data-sc16/geodesy/3819/cf_3819_GEHILMSTVY_mod", 
#                     ['EV', 'ME', 'MS', 'MT', 'MV', 'SE', 'SV', 'TE', 'TV'],
#                               ['LL+RR'])

# evme2 = gather_fringe_files(dlin, cf_lin, ['EV', 'ME'])


#ff_list = find_fringe_files("/data-sc16/geodesy/3819/")

ff_list = find_fringe_files(lin_3819)

# ff_list = find_fringe_files(cirI_3819, pol=["XX", "YY"])

nff = len(ff_list)
mbd = np.zeros(nff, dtype=float)

f_obj = ffm.FringeFileHandle()

for iff in range(nff):
    f_obj.load(ff_list[iff])
    mbd[iff] = f_obj.mbdelay
    print(ff_list[iff], ", mbdelay = ", mbd[iff])





