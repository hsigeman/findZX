import sys
import pandas as pn

if len(sys.argv)==1:
	print("Pyton-script written for the matchScaffold2Chr rule in snakemake pipe-line.")
	print("Joins two tab-separated files on the contig name and window start and stop.")
	print("Normalizes each sample on the mean of that sample and then takes the mean over each sex.\n")
	sys.exit()
elif not len(sys.argv)==3:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tmatchScaffold2chr_cov.py [bestMatch.list] [genCov] \n")
	sys.exit()

best_match = pn.read_csv(sys.argv[1], sep='\t', header=None)
best_match = best_match.drop(best_match.columns[[3, 4, 5, 6]], axis=1)

best_match.columns = [0,1,2,'chr','start','end']

gencov = pn.read_csv(sys.argv[2], sep='\t', header=None)

gencov_ref = pn.merge(best_match, gencov, how='inner', on=[0,1,2])

#nr_samples = gencov_ref.apply(max).tolist()

#nr_samples_each_sex = int((len(nr_samples) - 6)/2)

#samples = list(range(3, 3 + 2*nr_samples_each_sex))
#heterogametic_sex = list(range(3, 3 + nr_samples_each_sex))
#homogametic_sex = list(range(3 + nr_samples_each_sex, 3 + 2*nr_samples_each_sex))


#norm = gencov_ref.loc[:,samples].div(gencov_ref.loc[:,samples].mean())	
#mean_hetero = norm.loc[:,heterogametic_sex].mean(axis=1)
#mean_homo = norm.loc[:,homogametic_sex].mean(axis=1)


#mean_sexes = pn.concat([mean_hetero, mean_homo], axis=1)
#gencov_mean = pn.merge(gencov_ref.loc[:,[0,1,2,'chr','start','end']],mean_sexes, left_index=True, right_index=True)

#sys.stdout.write(gencov_mean.to_csv(header=None, index=None, sep='\t'))
sys.stdout.write(gencov_ref.to_csv(header=None, index=None, sep='\t'))

