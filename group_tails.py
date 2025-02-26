import numpy as np

#
# Group left-tail and right-tail bins with sparse data.
#
#       ni_grp, fni_grp = group_tails(ni, fni) 
#
# ni:  observed frequencies
# fni: expected (theoretical) frequencies
#
# 
#
# Returns ni and fni grouped.
#

# def group_tails(ni, fni):

#     nbin = len(ni)

#     # print('Original binning with sparse tails (%d bins):' % nbin)
#     # print('ni:  ', ni)
#     # print('fni: ', fni)


#     ltail_ni = 0;  ltail_fni = 0;  ltail_idx = 0
#     rtail_ni = 0;  rtail_fni = 0;  rtail_idx = nbin

#     for i in range(nbin):
#         if ni[i] <= 5:
#             ltail_idx = i
#             ltail_ni += ni[i]
#             ltail_fni += fni[i]
#             # print("i=%d, ni[i]=%d, ltail_ni=%d" % (i, ni[i], ltail_ni))
#             # if ltail_ni > 5:
#             #     break
#         else:
#             # ltail_idx stays at the leftmost frequency < 5
#             break

#     print("left i = ", i)
        
#     # if ltail_idx > 0:  # Set ltail_idx just before the leftmost freq < 5
#     #     ltail_idx -= 1
    
#     # print("DONE: ltail_idx=%d, ltail_ni=%d" % (ltail_idx, ltail_ni))

#     # print()

#     for i in range(nbin-1,-1,-1): # nbin-1 downto 0 
#         if ni[i] <= 5:
#             rtail_idx = i
#             rtail_ni += ni[i]
#             rtail_fni += fni[i]
#             # print("i=%d, ni[i]=%d, rtail_ni=%f" % (i, ni[i], rtail_ni))
#         else:
#             break

#     print("right i = ", i)

        
#     # print("DONE: rtail_idx=%d, rtail_ni=%f" % (rtail_idx, rtail_ni))
#     # print()


#     if rtail_idx < nbin:  # Set rtail_idx just after the rightmost freq < 5
#         rtail_idx += 1

#     #
#     # Group the tail data: cut the tails and place what was in the tails to
#     # bins 0 and -1 (end).
#     #

#     print("ltail_idx = %d; rtail_idx = %d" % (ltail_idx, rtail_idx))
#     # print("ltail_ni = %d; rtail_ni = %d" % (ltail_ni, rtail_ni))
    
#     ni_grp = np.copy(ni[ltail_idx : rtail_idx])   
#     fni_grp = np.copy(fni[ltail_idx : rtail_idx])
#     if ltail_idx > 0:
#         ni_grp[0] = ltail_ni
#         fni_grp[0] = ltail_fni
#     if rtail_idx < nbin:
#         ni_grp[-1] = rtail_ni
#         fni_grp[-1] = rtail_fni

#     return ni_grp, fni_grp



def find_tail_bounds(ni, thr=5):
    '''
    Find the indices at which the sparse tails are to be cut.
      ni: histogram array with frequencies
      thr: frequency threshold, below which the tail considered sparse.
    Returns ltail_idx and rtail_idx.
    '''
    nbin = len(ni)

    ltail_ni = 0;  ltail_fni = 0;  ltail_idx = 0
    rtail_ni = 0;  rtail_fni = 0;  rtail_idx = nbin

    for i in range(nbin):
        if ni[i] <= thr:
            ltail_idx = i
            ltail_ni += ni[i]
            ltail_fni += fni[i]
        else:
            # ltail_idx stays at the leftmost frequency < thr
            break

    for i in range(nbin-1,-1,-1): # nbin-1 downto 0 
        if ni[i] <= thr:
            rtail_idx = i
            rtail_ni += ni[i]
            rtail_fni += fni[i]
        else:
            break

    if rtail_idx < nbin:  # Set rtail_idx just after the rightmost freq < thr
        rtail_idx += 1

    return ltail_idx, rtail_idx


def group_tails(ni, lr_inds):
    '''
    Group the tail data: cut the tails and place what was in the tails to
    bins 0 and -1 (end).
      ni: histogram array with frequencies
      lr_inds: a sequence of left and right tail indices, (ltail_idx, rtail_idx)
    Returns ni with the sparse tails grouped.
    '''

    # print("ltail_idx = %d; rtail_idx = %d" % (ltail_idx, rtail_idx))
    # print("ltail_ni = %d; rtail_ni = %d" % (ltail_ni, rtail_ni))

    nbin = len(ni)
    ltail_idx, rtail_idx = lr_inds
    
    ni_grp = np.copy(ni[ltail_idx : rtail_idx])   

    if ltail_idx > 0:
        ltail_ni = ni[:ltail_idx+1]
        ni_grp[0] = np.sum(ltail_ni)
        print("ltail_ni = ", ltail_ni)
        print("ni_grp = ", ni_grp)
    if rtail_idx < nbin:
        rtail_ni = ni[rtail_idx-2:]
        ni_grp[-1] = np.sum(rtail_ni)

    return ni_grp



if __name__ == '__main__':

    ni = np.array([  1,   0,   0,   0,   0,   0,   5,
                     8,   8,  13,  31, 215, 147,
                     37,   7,   0,   0,   0,   0,   1,   1], dtype=float)
    fni = np.array([  0.,   0.,   0.,   0.,   0.,   0., 0.07, 0.96,
                      7.42, 33.44, 87.86, 134.48, 119.92, 62.3, 18.86, 3.32,
                      0.34, 0.02, 0., 0., 0.], dtype=float)

    
    
    np.set_printoptions(suppress=True, precision=2)

    print("ni: \n", ni)
    print()
    print("fni: \n", fni)
    print()
    
    nbin = len(ni)

    print("nbin = ", nbin)
    print()
    
#   ni_grp, fni_grp = group_tails(ni, fni)

    l_idx, r_idx = find_tail_bounds(ni)
    
    lr_inds = l_idx, r_idx
    
    ni_grp =  group_tails(ni, lr_inds)
    fni_grp = group_tails(fni, lr_inds)
    
    nbin_grp = len(ni_grp)
    
    print("l_idx, r_idx = find_tail_bounds(ni)\n")
    print("l_idx = %d; r_idx = %d" % (l_idx, r_idx))
    print()
    print("ni_grp =  group_tails(ni, lr_inds)")
    print("fni_grp =  group_tails(fni, lr_inds)\n")
    print("nbin_grp = ", nbin_grp)
    print()

    print("ni_grp: \n", ni_grp)
    print()

    print("fni_grp: \n", fni_grp)

    np.set_printoptions(suppress=False, precision=8)
