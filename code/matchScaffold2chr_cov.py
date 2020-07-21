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
gencov_ref = gencov_ref.drop(gencov_ref.columns[[3, 4, 5, 6, 10, 11, 12]], axis=1)

nr_samples = gencov_ref.apply(max).tolist()
nr_samples_each_sex = int((len(nr_samples) - 6)/2)

samples = list(range(13, 13 + 2*nr_samples_each_sex))
heterogametic_sex = list(range(13, 13 + nr_samples_each_sex))
homogametic_sex = list(range(13 + nr_samples_each_sex, 13 + 2*nr_samples_each_sex))

norm = gencov_ref.loc[:,samples].div(gencov_ref.loc[:,samples].mean())	
mean_hetero = norm.loc[:,heterogametic_sex].mean(axis=1)
mean_homo = norm.loc[:,homogametic_sex].mean(axis=1)

mean_sexes = pn.concat([mean_hetero, mean_homo], axis=1)
gencov_mean = pn.merge(gencov_ref.loc[:,[7,8,9]],mean_sexes, left_index=True, right_index=True)

sys.stdout.write(gencov_mean.to_csv(header=None, index=None, sep='\t'))

