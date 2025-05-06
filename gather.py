import processingle
from processingle import gather_fringe_files

exp_dir = '/data-sc16/geodesy/3819/'

cf = '/data-sc16/geodesy/3819/cf_3819_MESTVY_pstokes'

ff_list = gather_fringe_files(exp_dir, cf, ['EV','ME'],['I'])

len(ff_list)
print(ff_list[0].filename)

cf = '/data-sc16/geodesy/3819/ cf_3819_GEHILMSTVY_mod'
ff_list = gather_fringe_files(exp_dir, cf, ['EV','ME'], ['XX','XY','YX','YY'])

len(ff_list)

276






