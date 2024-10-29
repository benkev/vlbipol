import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib.patches as patches

def annotate_axes(ax, text, fontsize=18):
    ax.text(0.5, 0.5, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize, color="darkgrey")


#
# Modern method
# Added in version 3.7. ?
# Arranging multiple Axes in a Figure
# https://matplotlib.org/stable/users/explain/axes/arranging_axes.html
#
gs_kw1 = dict(width_ratios=[1.4, 1], height_ratios=[1, 2])
fig1, axd1 = pl.subplot_mosaic([['upper left', 'right'],
                               ['lower left', 'right']],
                              gridspec_kw=gs_kw1, figsize=(5.5, 3.5),
                              layout="constrained")
for k, ax in axd1.items():
    annotate_axes(ax, f'axd[{k!r}]', fontsize=14)
fig1.suptitle('pl.subplot_mosaic()') 


gs_kw2 = dict(width_ratios=[1, 1], height_ratios=[0.2, 0.6, 0.3])
fig2, axd2 = pl.subplot_mosaic([['legend', 'legend'],
                               ['left1', 'right1'],
                               ['left2', 'right2']],
                              gridspec_kw=gs_kw2, figsize=(8.4, 9),
                              layout="constrained")
for k, ax in axd2.items():
    annotate_axes(ax, f'axd2[{k!r}]', fontsize=14)
fig2.suptitle('pl.subplot_mosaic()') 

#
# Access subplots
#
ax = axd2['legend']
ax.plot([1,2])

ax = axd2['left2']
ax.plot(np.sin(np.linspace(-2*np.pi, 2*np.pi, 101)))
ax.set_axis_off()

#
# OLDER method from
#
# Customizing Figure Layouts Using GridSpec and Other Functions
# https://matplotlib.org/3.1.1/tutorials/intermediate/gridspec.html
#
fig3, axs3 = pl.subplots(ncols=2, nrows=3, figsize=(8.4, 9),
                         constrained_layout=True,
                         gridspec_kw={'height_ratios' : [0.2, 0.6, 0.3]})

for ax in axs3[0,:]:   # Remove all the 2 top axes
    ax.remove()

gs = axs3[0, 0].get_gridspec()  # Get GridSpec from ANY of the axes

ax0 = fig3.add_subplot(gs[0,:])  # Create one top axis





