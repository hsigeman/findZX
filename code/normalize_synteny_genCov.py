import sys
import pandas as pn

if len(sys.argv)==1 or sys.argv[1].startswith('-h'):
	print("\nScript written for snakemake pipeline for detection of sex-linked genomic regions. WARNING: USE WITH CAUTION OUTSIDE PIPELINE.\n")

	print("Reads in a genome coverage-file which has been intersected with a chromosome file. First each individual is normalized on the")
	print("individual means and then the mean for each sex and row is calculated. Writes to stdout.\n")

	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")

	print("Call:\tpython3 {python-script} {bed-file} {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()
elif not len(sys.argv)>=2:
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

