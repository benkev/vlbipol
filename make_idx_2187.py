help_text = '''
make_idx.py: Creates dictionaries to keep Mark4 data in convenient
             format. For both linearly- and circularly-polarized data,
             the following dictionaries (ending at 'l' and 'c', respectively)
             are created:

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
import os, sys, pickle, copy
from librd import make_idx
from libvp import make_closure_dic


#
# Linear polarization
#

# linI_2187 = "/media/benkev/Seagate_Backup_Plus_5TB_2/Work/" \
#             "2187/scratch/Lin_I/2187"

linI_2187 = "/home/benkev/Work/2187/scratch/Lin_I/2187"

idxl, idxsl, idxfl = make_idx(linI_2187)
print("Created dictionaries idxl, idxsl, and idxfl, linear polarization")

# sys.exit(0)

# with open('idx2187lI.pkl', 'wb') as fout: pickle.dump(idxl, fout)
# with open('idxs2187lI.pkl', 'wb') as fout: pickle.dump(idxsl, fout)
# with open('idxf2187lI.pkl', 'wb') as fout: pickle.dump(idxfl, fout)

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

# with open('idx2187cI.pkl', 'wb') as fout: pickle.dump(idxc, fout)
# with open('idxs2187cI.pkl', 'wb') as fout: pickle.dump(idxsc, fout)
# with open('idxf2187cI.pkl', 'wb') as fout: pickle.dump(idxfc, fout)

# print("Circular polarization data source:")
# print("   ", cirI_2187)
# print()
# print("Circular polarization data dictionaries pickled and saved on disk\n")

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

print("ST baseline excluded\n")
print("Found %d baselines common for linear and circular polarization:" %
      nbls)
print(bls, '\n')
print()

closl = make_closure_dic(idxsl, bls)
closc = make_closure_dic(idxsc, bls)

print("Created dictionary of data closures, closl, linear polarization")
print("Created dictionary of data closures, closc, circular polarization\n")

# with open('clos2187lI.pkl', 'wb') as fout: pickle.dump(closl, fout)
# with open('clos2187cI.pkl', 'wb') as fout: pickle.dump(closc, fout)

# print("Linear pol. data closure dictionary pickled and saved on disk")
# print("Circular pol. data closure dictionary pickled and saved on disk\n")

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






