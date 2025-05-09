'''
plot_dtec_2187.py
'''
import pickle, copy
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import libvp  # Needs to reset backend for vpal sets it to non-interactive Agg!

#
# Unpickle it:
#
# with open('idx2187lI.pkl', 'rb') as finp: idxl = pickle.load(finp)
# with open('idx2187cI.pkl', 'rb') as finp: idxc = pickle.load(finp)

with open('idxs2187lI.pkl', 'rb') as finp: idxsl = pickle.load(finp)
with open('idxs2187cI.pkl', 'rb') as finp: idxsc = pickle.load(finp)

# with open('idxf2187lI.pkl', 'rb') as finp: idxfl = pickle.load(finp)
# with open('idxf2187cI.pkl', 'rb') as finp: idxfc = pickle.load(finp)

# with open('clos2187lI.pkl', 'rb') as finp: closl = pickle.load(finp)
# with open('clos2187cI.pkl', 'rb') as finp: closc = pickle.load(finp)

with open('bls_2187.pkl', 'rb') as finp: bls = pickle.load(finp)
with open('tribl_2187.pkl', 'rb') as finp: tribl = pickle.load(finp)

#
# Plotting for the 1803+784 source
#
sr = '1803+784'


tSE = []
tTE = []
dtSE = []
dtTE = []

for tm in idxsl[sr].keys():
    if 'SE' in idxsl[sr][tm].keys():
        tSE.append(tm/3600)
        dtSE.append(idxsl[sr][tm]['SE']['dtec'])
    if 'TE' in idxsl[sr][tm].keys():
        tTE.append(tm/3600)
        dtTE.append(idxsl[sr][tm]['TE']['dtec'])

pl.figure()

pl.plot(tSE, dtSE, 'r.', label='SE')
pl.plot(tTE, dtTE, 'g.', label='TE')
pl.grid(1)
pl.title("VO2187, Source %8s, dTEC for Baselines" % sr)
pl.xlabel("t (hours)")
pl.ylabel("dTEC")
pl.legend()
pl.savefig("VO2187_Source_%8s_dTEC_for_Baselines_%s_and_%s.pdf" %
           (sr, 'SE', 'TE'),
           format='pdf')


tGS = []
tGT = []
dtGS = []
dtGT = []

for tm in idxsl[sr].keys():
    if 'GS' in idxsl[sr][tm].keys():
        tGS.append(tm/3600)
        dtGS.append(idxsl[sr][tm]['GS']['dtec'])
    if 'GT' in idxsl[sr][tm].keys():
        tGT.append(tm/3600)
        dtGT.append(idxsl[sr][tm]['GT']['dtec'])

pl.figure()

pl.plot(tGS, dtGS, 'r.', label='GS')
pl.plot(tGT, dtGT, 'g.', label='GT')
pl.grid(1)
pl.title("VO2187, Source %8s, dTEC for Baselines" % sr)
pl.xlabel("t (hours)")
pl.ylabel("dTEC")
pl.legend()
pl.savefig("VO2187_Source_%8s_dTEC_for_Baselines_%s_and_%s.pdf" %
           (sr, 'GS', 'GT'),
           format='pdf')




tHS = []
tHT = []
dtHS = []
dtHT = []

for tm in idxsl[sr].keys():
    if 'HS' in idxsl[sr][tm].keys():
        tHS.append(tm/3600)
        dtHS.append(idxsl[sr][tm]['HS']['dtec'])
    if 'HT' in idxsl[sr][tm].keys():
        tHT.append(tm/3600)
        dtHT.append(idxsl[sr][tm]['HT']['dtec'])

pl.figure()

pl.plot(tHS, dtHS, 'r.', label='HS')
pl.plot(tHT, dtHT, 'g.', label='HT')
pl.grid(1)
pl.title("VO2187, Source %8s, dTEC for Baselines" % sr)
pl.xlabel("t (hours)")
pl.ylabel("dTEC")
pl.legend()
pl.savefig("VO2187_Source_%8s_dTEC_for_Baselines_%s_and_%s.pdf" %
           (sr, 'HS', 'HT'),
           format='pdf')





tIS = []
tIT = []
dtIS = []
dtIT = []

for tm in idxsl[sr].keys():
    if 'IS' in idxsl[sr][tm].keys():
        tIS.append(tm/3600)
        dtIS.append(idxsl[sr][tm]['IS']['dtec'])
    if 'IT' in idxsl[sr][tm].keys():
        tIT.append(tm/3600)
        dtIT.append(idxsl[sr][tm]['IT']['dtec'])

pl.figure()

pl.plot(tIS, dtIS, 'r.', label='IS')
pl.plot(tIT, dtIT, 'g.', label='IT')
pl.grid(1)
pl.title("VO2187, Source %8s, dTEC for Baselines" % sr)
pl.xlabel("t (hours)")
pl.ylabel("dTEC")
pl.legend()
pl.savefig("VO2187_Source_%8s_dTEC_for_Baselines_%s_and_%s.pdf" %
           (sr, 'IS', 'IT'),
           format='pdf')






tMS = []
tMT = []
dtMS = []
dtMT = []

for tm in idxsl[sr].keys():
    if 'MS' in idxsl[sr][tm].keys():
        tMS.append(tm/3600)
        dtMS.append(idxsl[sr][tm]['MS']['dtec'])
    if 'MT' in idxsl[sr][tm].keys():
        tMT.append(tm/3600)
        dtMT.append(idxsl[sr][tm]['MT']['dtec'])

pl.figure()

pl.plot(tMS, dtMS, 'r.', label='MS')
pl.plot(tMT, dtMT, 'g.', label='MT')
pl.grid(1)
pl.title("VO2187, Source %8s, dTEC for Baselines" % sr)
pl.xlabel("t (hours)")
pl.ylabel("dTEC")
pl.legend()
pl.savefig("VO2187_Source_%8s_dTEC_for_Baselines_%s_and_%s.pdf" %
           (sr, 'MS', 'MT'),
           format='pdf')








pl.show()
