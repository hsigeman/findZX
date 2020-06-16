import sys
import pandas as pn

gencov = pn.read_csv(sys.argv[1], header=None, sep='\t')

gmax = gencov.apply(max).tolist()

nr_samples_each_sex = int((len(gmax) - 6)/2)

samples = list(range(6, len(gmax)))
females = list(range(6, nr_samples_each_sex + 6))
males = list(range(6 + nr_samples_each_sex, len(gmax)))

#norm = gencov.loc[:,samples].div(gencov.loc[:,samples].max())
norm = gencov.loc[:,samples].div(gencov.loc[:,samples].mean())	#normalize with mean instead of max-value
mean_f = norm.loc[:,females].mean(axis=1)
mean_m = norm.loc[:,males].mean(axis=1)


mean_fm = pn.concat([mean_f, mean_m], axis=1)

gencov_mean = pn.merge(gencov.loc[:,[0,1,2,3,4,5]],mean_fm, left_index=True, right_index=True)

sys.stdout.write(gencov_mean.to_csv(header=None, index=None, sep='\t'))

