import os
import re
# import numpy as np
# import matplotlib.pyplot as plt

import hopstestb as ht
# import ffcontrol
# from vpal import fringe_file_manipulation as ffm


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

            # for pp in pp_list:
            #     if pp in pol:
            #         ff_list.append(full_name)
            #         #print("pp_list = ", pp_list, ", ", full_name)

            if pp_list == pol:
                ff_list.append(full_name)

    return ff_list

