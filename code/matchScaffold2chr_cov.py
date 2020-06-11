import sys
import pandas as pn

best_match = pn.read_csv(sys.argv[1], sep='\t', header=None)
best_match = best_match.drop(best_match.columns[[3, 4, 5, 6]], axis=1)

gencov = pn.read_csv(sys.argv[2], sep='\t', header=None)

gencov_zf = pn.merge(best_match, gencov, how='inner')

sys.stdout.write(gencov_zf.to_csv(header=None, index=None, sep='\t'))

