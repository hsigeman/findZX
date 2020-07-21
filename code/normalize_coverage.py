import sys
import pandas as pn

if len(sys.argv)==1:
	print("Pyton-script written for the matchScaffold2Chr rule in snakemake pipe-line.")
	print("Normalizes each sample on the mean of that sample and then takes the mean over each sex.\n")
	sys.exit()
elif not len(sys.argv)==2:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tmatchScaffold2chr_cov.py [genCov] \n")
	sys.exit()

gencov_ref = pn.read_csv(sys.argv[1], sep='\t', header=None)

nr_samples = gencov_ref.apply(max).tolist()
nr_samples_each_sex = int((len(nr_samples) - 3)/2)

samples = list(range(3, 3 + 2*nr_samples_each_sex))
heterogametic_sex = list(range(3, 3 + nr_samples_each_sex))
homogametic_sex = list(range(3 + nr_samples_each_sex, 3 + 2*nr_samples_each_sex))

norm = gencov_ref.loc[:,samples].div(gencov_ref.loc[:,samples].mean())	
mean_hetero = norm.loc[:,heterogametic_sex].mean(axis=1)
mean_homo = norm.loc[:,homogametic_sex].mean(axis=1)

mean_sexes = pn.concat([mean_hetero, mean_homo], axis=1)
gencov_mean = pn.merge(gencov_ref.loc[:,[0,1,2]],mean_sexes, left_index=True, right_index=True)

sys.stdout.write(gencov_mean.to_csv(header=None, index=None, sep='\t'))

