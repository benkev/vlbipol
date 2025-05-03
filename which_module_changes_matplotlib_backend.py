'''
Determine which imported module changes matplotlib backend.

Usage:

python which_module_changes_matplotlib_backend.py mod1 mod2 mod3 ...
       where mod1, mod2 etc. are the module names in imports.

Prints:

Before importing: qtagg

Importing mod1 ...
Backend now: <backend>  (qtagg, Agg, or whatever else)

Importing libplt ...
Backend now: <backend>  (qtagg, Agg, or whatever else)
.  .  .  .  .  .  .  .

'''

import importlib
import matplotlib
import sys

# modlist = ['hopstestb', 'ffcontrol', 'vpal.fringe_file_manipulation',
#                 'vpal.utility', 'mk4b']

if len(sys.argv) == 1:
    print("Determine which imported module changes matplotlib backend.")
    print("Usage:")
    print("python which_module_changes_matplotlib_backend.py mod1 mod2 ...")
    print("       where mod1, mod2 etc. are the module names in imports.")
    sys.exit(0)

modlist = sys.argv[1:]   

# Reset to qtagg for clarity
print("Before importing:", matplotlib.get_backend())

# Test imports one by one
for modname in modlist:
    print(f"\nImporting {modname} ...")
    importlib.import_module(modname)
    print("Backend now:", matplotlib.get_backend())
