'''
libplt.py - functions for plotting
'''

import matplotlib.pyplot as pl
import matplotlib.patches as patches


def plot_closure_legend(ax_col, trians, cols, par, fs=12):
    '''
    Plot color legend for the closure triangles in a separate axis ax_col.
    For 16 closure triangles, it plots 4 colimns by 4 in height, each showing
    the color and the triangle name on the right.
    Inputs:
        ax_col: axis; better to have shorter in height and wider.
                It is assumed to be on top of other plots and its title
                describes the whole figure.
        trians: list of 3-letter closure triangle names.
        cols:   Array of colors, ntri by 4, from pyplot.cm(numbers 0 to 1)
                For example:
                    cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))
        par:    The parameter name to print in the title in uppercase.
        fs:     Fontsise of the printed triangle names
    '''

    #pl.ion()  # Interactive mode; pl.ioff() - revert to non-interactive.

    ntri = len(trians)
    qntri = ntri // 4        # Quarter of the number of triangles
    #rntri = ntri % 4         # Residue of triangles
    #cols = cm.nipy_spectral(1 - np.linspace(0, 1, ntri))           # Colors
    upar = par.upper() # The uppercase parameter name to print in the title.
    
    ax_col.set_xlim(-0.1, 7.5)
    # ax_col.set_ylim(-0.2, qntri+0.5)
    ax_col.set_ylim(-1.2, qntri+0.5)

    #qntri = qntri + 1      # Leave room for the last triangle
    for i in range(qntri):
        y = qntri - i - 1      # Vertical patch position
        k = i                  # Index into trians[k] and cols[k,:]
        rect = patches.Rectangle((0, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(1.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + qntri
        rect = patches.Rectangle((2, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(3.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 2*qntri
        rect = patches.Rectangle((4, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(5.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

        k = i + 3*qntri
        rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
        ax_col.add_patch(rect)
        ax_col.text(7.1, y+0.2, trians[k], fontsize=fs)
        # print("i = %d, y = %d, k = %2d" % (i, y, k))

    # Remaining triangle
    k = 28
    y = -1
    rect = patches.Rectangle((6, y), .95, .95, facecolor=cols[k,:])
    ax_col.add_patch(rect)
    ax_col.text(7.1, y+0.2, trians[k], fontsize=fs)

    ax_col.set_axis_off()
    ax_col.set_title("%s Closure" % upar, fontsize=16)


