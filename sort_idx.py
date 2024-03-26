help_text = '''
sort_idx.py: Sort file names in the file lista in index file time ascending
'''

import pickle
import hopstestb as ht
import ffcontrol
from vpal import fringe_file_manipulation as ffm

#
# Unpickle the index file:
#
with open('idx1.pkl', 'rb') as finp:
    idx1 = pickle.load(finp)

ttag = []

f_obj = ffm.FringeFileHandle()

#
# Iterate over baselines and, for each baseline, over polarizations
#
# for bl in idx1:
#     print(bl)
#     for po in idx1[bl]:
#         print("    ", po)

for ff in idx1['SE']['I']:
    f_obj.load(ff)
    ttag.append(f_obj.time_tag)












