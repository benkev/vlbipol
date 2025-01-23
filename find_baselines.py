import os, sys
import re, copy


cirI_2187 = "/home/benkev/Work/vo2187_exprm/DiFX_pconv/2187"

base_dir = copy.copy(cirI_2187)
max_depth=2

# Remove trailing "/", if any (os.path.sep is usually "/")
base_dir = base_dir.rstrip(os.path.sep)
num_sep = base_dir.count(os.path.sep)

bl_cor = set()
bl_ff = set()

for root_dir, subdirs, files in os.walk(base_dir):

    if root_dir.count(os.path.sep) > num_sep + max_depth:
        continue

    for file in files:
        if re.match(r"[A-Z]{2}\.\.\w{6}", file):
            bl_cor.add(file[:2])
        if re.match(r"[A-Z]{2}\.X\.[0-9]+\.\w{6}", file):
            bl_ff.add(file[:2])






