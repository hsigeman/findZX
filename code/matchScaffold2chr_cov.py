import sys
import panda as pn

best_match = pn.read_csv(argv[1], sep='\t')
cov = pn.read_csv(argv[2], sep='\t')

gencov_zf = pn.merge(best_match, cov, how='inner', left_on=['1','2','3'], right_on=['1','2','3'])

#drop columns that are not needed

