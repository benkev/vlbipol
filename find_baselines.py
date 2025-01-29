import os, sys
import re, copy


cirI_2187 = "/home/benkev/Work/vo2187_exprm/DiFX_pconv/2187"
lin_2187 = "/home/benkev/Work/2187/scratch/Lin_1/2187"

# This dir has only the 7 baselines with the G station:
#lin_2187 = "/home/benkev/Work/2187/scratch/20241221-010014/2187"
#lin_2187 = "/home/benkev/Work/2187"


#base_dir = copy.copy(cirI_2187)
base_dir = copy.copy(lin_2187)
max_depth=2

# Remove trailing "/", if any (os.path.sep is usually "/")
base_dir = base_dir.rstrip(os.path.sep)
num_sep = base_dir.count(os.path.sep)

bl_cor = set()
bl_ff = set()

for root_dir, subdirs, files in os.walk(base_dir):

    if root_dir.count(os.path.sep) > num_sep + max_depth:
        continue

    for fn in files:
        if re.match(r"[A-Z]{2}\.\.\w{6}", fn) and fn[0] != fn[1]:
            bl_cor.add(fn[:2])
        if re.match(r"[A-Z]{2}\.X\.[0-9]+\.\w{6}", fn) and fn[0] != fn[1]:
            bl_ff.add(fn[:2])


bl_cor = list(bl_cor)
bl_ff = list(bl_ff)

bl_cor.sort()
bl_ff.sort()

print("len(bl_cor) = ", len(bl_cor))
print("len(bl_ff) = ", len(bl_ff))



