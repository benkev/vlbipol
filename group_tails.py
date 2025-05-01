#
# group_tails.py - sparse histogram tails grouping module.
#

import numpy as np

def find_tail_bounds(ni, thr=5):
    '''
    Find histogram indices at which the sparse tails are to be cut.
      ni: histogram array with frequencies
      thr: frequency threshold, below which the tail considered sparse.
    Returns l_idx and r_idx, where
      l_idx points AT the first ni element greater than thr;
      r_idx points AFTER the last ni element greater than thr,
   so ni[l_idx : r_idx] is the subarray whose elements are > thr.
    
    '''
    nbin = len(ni)

    l_idx = 0
    r_idx = nbin

    for i in range(nbin):
        if ni[i] > thr:
            l_idx = i
            break

    for i in range(nbin-1,-1,-1): # nbin-1 downto 0 
        if ni[i] > thr:
            r_idx = i + 1
            break

    return l_idx, r_idx


def group_tails(ni, lr_inds):
    '''
    Group the tail data in histogram array ni:
    cut the tails and place what was in the tails to
    bins 0 and -1 (end).
      ni: histogram array with frequencies
      lr_inds: a sequence of two indices, (l_idx, r_idx), left and right.
      l_idx points AT the first ni element greater than thr;
      r_idx points AFTER the last ni element greater than thr,
    so ni[l_idx : r_idx] is the subarray whose elements are > thr.
    
    Returns ni with the sparse tails grouped.
    '''

    nbin = len(ni)
    l_idx, r_idx = lr_inds
    
    ni_grp = np.copy(ni[l_idx : r_idx])   

    if l_idx > 0:
        l_ni = ni[:l_idx+1]
        ni_grp[0] = np.sum(l_ni)

    if r_idx < nbin:
        r_ni = ni[r_idx-1:]
        ni_grp[-1] = np.sum(r_ni)

    return ni_grp



if __name__ == '__main__':

    ni = np.array([  1,   0,   0,   0,   0,   0,   5,
                     8,   8,  13,  31, 215, 147,
                     37,   7,   0,   0,   0,   0,   1,   1], dtype=int)
    fni = np.array([  0.,   0.,   0.,   0.,   0.,   0., 0.07, 0.96,
                      7.42, 33.44, 87.86, 134.48, 119.92, 62.3, 18.86, 3.32,
                      0.34, 0.02, 0., 0., 0.], dtype=float)

    
    
    np.set_printoptions(suppress=True, precision=2)

    print("ni = ", ni)
    
    nbin = len(ni)

    print("nbin = ", nbin)
    print()
    
    l_idx, r_idx = find_tail_bounds(ni)
    
    lr = l_idx, r_idx
    
    print("l_idx, r_idx = find_tail_bounds(ni)")
    print("l_idx = %d; r_idx = %d" % (l_idx, r_idx))
    print()
    
    ni_grp =  group_tails(ni, lr)
    # fni_grp = group_tails(fni, lr)
    
    nbin_grp = len(ni_grp)
    
    print("\nni_grp =  group_tails(ni, lr)")
    # print("fni_grp =  group_tails(fni, lr)\n")
    print("nbin_grp = ", nbin_grp)
    print()

    print("ni_grp: \n", ni_grp)
    print()

    print(list(zip(range(len(ni)), ni)))

    np.set_printoptions(suppress=False, precision=8)
