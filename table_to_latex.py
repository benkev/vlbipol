#
# table_to_latex.pdf
#
# Slightly helps to format a text table for using it in the LaTeX closures
# \begin{table}
#    \begin{tabular}{...}
#
#    \end{tabular}
# \end{table}
#
# $ python table_to_latex.pdf <input-text-file-with-table>
#
# Reads a text file with space-separated columns and prints its lines
# with & signs inserted and trailing "\\".
#
# For example, a line
#      HE     82.2    710.6    645.6    90.86      198.2    0.759337
# turns into
#      HE &  82.2 & 710.6 &  645.6$ &  90.86 &   198.2 & 0.759337 \\
#

import sys

tabf = sys.argv[1]

with open(tabf, 'r') as finp:
    tin = finp.read()
tin = tin.split('\n')

for tl in tin:
    els = tl.split()
    els = [el.rjust(8) for el in els]
    lin = ' & '.join(els) + " \\\\"
    print(lin)
    
