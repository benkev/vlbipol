import pickle, copy, sys
import libvp
from libvp import find_baseline_triangles

np.set_printoptions(precision=4,suppress=True)

with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)
with open('bls_2187.pkl', 'rb') as finp: bls = pickle.load(finp)

trians = find_baseline_triangles(bls)

#
# clot[tm][sr][tr]
#
                
clot = {tr: {} for tr in trians.keys()} 

for sr in closl.keys():
    for tr in closl[sr].keys():
        for di in closl[sr][tr].keys():
            if sr in clot[tr].keys():
                clot[tr][sr][di] = copy.deepcopy(closl[sr][tr][di])
            else:
                clot[tr][sr] = {}
                clot[tr][sr][di] = copy.deepcopy(closl[sr][tr][di])

#
# To save:
#
# with open('clot2187lI.pkl', 'wb') as fout: pickle.dump(clot, fout)

# To load:
#
# with open('clot2187lI.pkl', 'rb') as fout: clotl = pickle.load(fout)
#




sys.exit(0)
                
#
# clotm[tm][sr][tr]   ??????????????????????
#
                
#
# Gather in tms a sorted list of all times in seconds
#
tms = set()
for bl in idxl.keys():
    tms = tms | set(idxl[bl]['I']['time'])
tms = list(tms)
tms.sort()

                
clotm = {tm: {} for tm in tms} 

for sr in closl.keys():
    for tr in closl[sr].keys():
        for di in closl[sr][tr].keys():
            if sr in clotm[tm].keys():
                clotm[tr][sr][di] = copy.deepcopy(closl[sr][tr][di])
            else:
                clotm[tr][sr] = {}
                clotm[tr][sr][di] = copy.deepcopy(closl[sr][tr][di])




                
