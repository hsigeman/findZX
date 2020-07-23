import sys
import pandas as pn

if len(sys.argv)==1:
	print("Pyton-script written for the matchScaffold2Chr rule in snakemake pipe-line.")
	print("Normalizes each sample on the mean of that sample and then takes the mean over each sex.\n")
	sys.exit()
elif not len(sys.argv)>=2:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tmatchScaffold2chr_cov.py [genCov] {heterogametic samples. prefix 'het:'}")
	print("{homogametic samples. prefix 'homo:'} \n")
	sys.exit()


nr_samples = int(len(sys.argv) - 1)

heterogametic_samples = []
homogametic_samples = []

for i in range(2, len(sys.argv)):
        sample_args = sys.argv[i].split(':')

        if sys.argv[i].startswith('het'):
                heterogametic_samples.append(sample_args[1])

        elif sys.argv[i].startswith('homo'):
                homogametic_samples.append(sample_args[1])

nr_heterogametic = len(heterogametic_samples)
nr_homogametic = len(homogametic_samples)


gencov_ref = pn.read_csv(sys.argv[1], sep='\t', header=None)
gencov_ref = gencov_ref.drop(gencov_ref.columns[[3, 4, 5, 6, 10, 11, 12]], axis=1)

samples = list(range(13, 13 + nr_samples))

heterogametic_sex = list(range(13, 13 + nr_heterogametic))
homogametic_sex = list(range(13 + nr_heterogametic, 13 + nr_heterogametic + nr_homogametic))

norm = gencov_ref.loc[:,samples].div(gencov_ref.loc[:,samples].mean())	
mean_hetero = norm.loc[:,heterogametic_sex].mean(axis=1)
mean_homo = norm.loc[:,homogametic_sex].mean(axis=1)

mean_sexes = pn.concat([mean_hetero, mean_homo], axis=1)
gencov_mean = pn.merge(gencov_ref.loc[:,[7,8,9]],mean_sexes, left_index=True, right_index=True)

sys.stdout.write(gencov_mean.to_csv(header=None, index=None, sep='\t'))

