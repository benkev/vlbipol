#
# Merge left-tail and right-tail bins into one bin
#
import numpy as np

np.set_printoptions(suppress=True, precision=2)

nbin = 21
ni = np.array([  1,   0,   0,   0,   0,   0,   0,   8,   8,  13,  31, 215, 147,
              37,   7,   0,   0,   0,   0,   1,   1], dtype=float)
fni = np.array([  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.07,   0.96,
         7.42,  33.44,  87.86, 134.48, 119.92,  62.3 ,  18.86,   3.32,
         0.34,   0.02,   0.  ,   0.  ,   0.  ], dtype=float)




ltail_ni = 0; ltail_fni = 0;  ltail_idx = 0
rtail_ni = 0; rtail_fni = 0;  rtail_idx = 0

for i in range(nbin):
    if ni[i] <= 5:
        ltail_idx = i
        ltail_ni += ni[i]
        ltail_fni += fni[i]
        # print("%d: ltail_ni = %6f; ltail_fni = %6f;  ltail_idx = %3d" % \
        #       (i, ltail_ni, ltail_fni, ltail_idx))
    else:
        break

ltail_idx = ltail_idx + 1 

print()
        
for i in range(nbin-1,-1,-1): # nbin-1 downto 0 
    if ni[i] <= 5:
        rtail_idx = i
        rtail_ni += ni[i]
        rtail_fni += fni[i]
        # print("%d: rtail_ni = %6f; rtail_fni = %6f;  rtail_idx = %3d" % \
        #       (i, rtail_ni, rtail_fni, rtail_idx))
    else:
        break

#print()

print("ltail_ni = %6f; ltail_fni = %6f;  ltail_idx = %3d" % \
      (ltail_ni, ltail_fni, ltail_idx))
print("rtail_ni = %6f; rtail_fni = %6f;  rtail_idx = %3d" % \
      (rtail_ni, rtail_fni, rtail_idx))
print()
#
# Group the tail data: cut the tails and add what was in the tails to
# bins 0 and end.
#
print("ni: \n", ni)
print("fni: \n", fni)
print()

ni_grp =   ni[ltail_idx : rtail_idx]   
ni_grp[0]  += ltail_ni
ni_grp[-1] += rtail_ni

fni_grp = fni[ltail_idx : rtail_idx]
fni_grp[0]  += ltail_fni
fni_grp[-1] += rtail_fni

nbin_grp = len(ni_grp)

# print("ni: \n", ni)
print("ni_grp: \n", ni_grp)
# print()

#print("fni: \n", fni)
print("fni_grp: \n", fni_grp)


np.set_printoptions(suppress=False, precision=8)
