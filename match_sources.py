#
# match_sources.py
#

import re, os, sys

base_dir = "/media/benkev/Seagate_Backup_Plus_5TB_2/Work" \
           "/2187/scratch/Lin_I/2187"

max_depth=2

num_sep = base_dir.count(os.path.sep)

nf = 0

for root_dir, subdirs, files in os.walk(base_dir):
    
    if root_dir.count(os.path.sep) > num_sep + max_depth:
        continue

    for file in files:

        mobj = re.fullmatch(r"^[0-9A-Z+-]{5,8}\.[0-9A-Z]{6}$", file)

        if mobj is None:
            continue

        nf = nf + 1

        src = file.split('.')[0]

        dir = root_dir.split('/')[-1]
        full_name = os.path.join(dir, file)
        
        
        print("%s   Source: %s" % (full_name, src))

print("\n%d files with source name found." % nf)
