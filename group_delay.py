#
# group_delay.py
#
# Get multiband/group delay data before and after PolConvert 
#

hops_a = "/data-sc16/geodesy/3819/251-1130"
hops_b = "/data-sc16/geodesy/3819/polconvert/3819/scratch/" \
         "pol_prods1/3819/251-1130"

import mk4b

mf_a = mk4b.mk4fringe(hops_a + 'EV.X.       ?????????????????  ')
mf_b = mk4b.mk4fringe(hops_b + 'EV.X.5.2HMQIA')

pcal = []
for ii in range(32):
    pcal = pcal.append(mf_a.t207.contents.ref_pcamp[ii].usb)
