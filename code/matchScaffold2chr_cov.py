import sys
import panda as pn

best_match =  argv[1]
cov = argv[2]

#read in best_match as a data_frame
#read in cov as a data_fram

gencov_zf = pn.merge(best_match, cov, how='inner', left_on=['1','2'], right_on=['1','2'])

#drop columns that are not needed

