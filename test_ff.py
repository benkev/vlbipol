import os
import hopstestb as ht
import ffcontrol

lin_3819 = "/data-sc16/geodesy/3819/" 
cir_3819 = "/data-sc16/geodesy/3819/polconvert/3819/scratch/pol_prods1/3819"

dlin = lin_3819 + "251-1130/"
dpol = cir_3819 + "251-1130/"

cf_lin = lin_3819 + "cf_3819_MESTVY_pstokes"
cf_cir = cir_3819 + "cf_3819_MESTVY_pcphase_dh2"

from vpal import fringe_file_manipulation as ffm
import hopstestb as ht
from vpal.processing import gather_fringe_files


# Find out if the file is for 'I'
ht.get_file_polarization_product_provisional(dlin + 'EV.X.29.2GI1S4')

f_obj = ffm.FringeFileHandle()

f_obj.load(dlin+ 'EV.X.29.2GI1S4')

mbd = f_obj.mbdelay

evme1 = gather_fringe_files("/data-sc16/geodesy/3819/", 
                    "/data-sc16/geodesy/3819/cf_3819_MESTVY_pstokes", 
                    ['EV', 'ME'])

evme2 = gather_fringe_files(dlin, cf_lin, ['EV', 'ME'])


# def gather_fringe_files(base_directory, control_file, blines, pol_products=['I'], include_autos=False, exclude_list=None, max_depth=2):
#     """
#     Returns a list of all the fringe (type2) files found in any directory up to max_depth under the base_directory
#     In a typical VGOS experiment directory, the fringe files are two levels below the base_directory (the first
#     level below the base directory is the scan name).

#     The exp_dir is expecting an absolute path to a HOPS experiment number directory (four-digit number)
#     The control file is an absolute path

#     Assumes fringe files have a six-digit extension with three '.' in the filename.
#     """
#     if exclude_list == None:
#         exclude_list=['prepass', 'pre_production', 'make_links', 'test', 'bad_eop', 'setup']
#     exclude = set(exclude_list)
#     assert os.path.isdir(base_directory)

#     fringe_file_list = []

#     base_dir = base_directory.rstrip(os.path.sep)
#     num_sep = base_dir.count(os.path.sep)

#     control_file_hash = ffcontrol.get_control_file_hash(control_file)

#     counter=0
#     for current_root, subdirectories, files in os.walk(base_directory, topdown=True):

#         for filename in files:

#             # apply the exclude filter
#             if any([e in current_root for e in exclude]):
#                 continue

#             # apply the max depth filter
#             if current_root.count(os.path.sep) > num_sep+max_depth:
#                 continue

#             abs_filename = os.path.abspath(filename)
#             filename_base = os.path.split(abs_filename)[1]


#             #look for root files using some simple checks
#             if filename_base.count('.') == 3: #check that there are three dots in the filename base
#                 bline = filename_base.split('.')[0]
#                 if len(bline)==2 and bline in blines: #make sure leading section of file name is 2-char baseline
#                     full_name = os.path.join(current_root, filename)
#                     if (include_autos is True) or (bline[0] != bline[1]): #check that this is a cross correlation if autos excluded
#                         extension = filename_base.split('.')[3] #get the file extension (root_id)
#                         if len(extension) == 6:     #check that the extension has a length of 6 chars

#                             counter+=1

#                             ff_cf_hash = ht.get_control_file_hash_from_fringe(full_name)

#                             if control_file_hash == ff_cf_hash:
#                                 fringe_file_list.append(full_name)


#     # Now we have a list of fringe files that used the prescribed control file.  Check for the correct polproduct:
#     ff_list = []
    
#     # apply the check on pol products
#     for ff in fringe_file_list:
#         ff_pp_list = ht.get_file_polarization_product_provisional(ff)

#         for pp in ff_pp_list:
#             if pp in pol_products:

#                 # Now that we are sure we have the correct polproduct, create a FringeFileHandle object and append to the list
#                 f_obj = ffm.FringeFileHandle()
#                 f_obj.load(ff)
#                 ff_list.append(f_obj)


#     #print(base_directory, 'Number of fringe files considered:', counter)
#     return ff_list

