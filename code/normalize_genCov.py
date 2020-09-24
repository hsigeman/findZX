import sys
import pandas as pn

if len(sys.argv)==1 or sys.argv[1].startswith('-h'):
	print("\nScript written for snakemake pipeline for detection of sex-linked genomic regions. WARNING: USE WITH CAUTION OUTSIDE PIPELINE.\n")

	print("Reads in a genome coverage-file. First each individual is normalized on the individual means and then the mean for each site and sex")
	print("is calculated. Writes to stdout.\n")

	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")

	print("Call:\tpython3 {python-script} {bed-file} {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()
elif not len(sys.argv)>=3:
	print("\nERROR: wrong number of input arguments!\n")
	
	print("Call:\tpython3 {python-script} {bed-file} {heterogametic-samples} {homogametic-samples}\n")
	
	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")
	sys.exit()

nr_samples = int(len(sys.argv) - 2)

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

samples = list(range(3, 3 + nr_samples))
heterogametic_sex = list(range(3, 3 + nr_heterogametic))
homogametic_sex = list(range(3 + nr_heterogametic, 3 + nr_heterogametic + nr_homogametic))

norm = gencov_ref.loc[:,samples].div(gencov_ref.loc[:,samples].mean())	
mean_hetero = norm.loc[:,heterogametic_sex].mean(axis=1)
mean_homo = norm.loc[:,homogametic_sex].mean(axis=1)

mean_sexes = pn.concat([mean_hetero, mean_homo], axis=1)
gencov_mean = pn.merge(gencov_ref.loc[:,[0,1,2]],mean_sexes, left_index=True, right_index=True)

sys.stdout.write(gencov_mean.to_csv(header=None, index=None, sep='\t'))

