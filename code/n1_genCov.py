import sys
import pandas as pn

if len(sys.argv)==1 or sys.argv[1].startswith('-h'):
	print("\nScript written for snakemake pipeline for detection of sex-linked genomic regions. WARNING: USE WITH CAUTION OUTSIDE PIPELINE.\n")

	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")

	print("Call:\tpython3 {python-script} {bed-file} with/no-synteny {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()
elif not len(sys.argv)>=3:
	print("\nERROR: wrong number of input arguments!\n")
	
	print("Call:\tpython3 {python-script} {bed-file} with/no-synteny {heterogametic-samples} {homogametic-samples}\n")
	
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


gencov = pn.read_csv(sys.argv[1], sep='\t', header=None)

sample_index_start = 2
chr_columns = [0,1,2]


samples = list(range(sample_index_start, sample_index_start + nr_samples))
heterogametic_sex = list(range(sample_index_start, sample_index_start + nr_heterogametic))
homogametic_sex = list(range(sample_index_start + nr_heterogametic, sample_index_start + nr_heterogametic + nr_homogametic))

# Normalize and take mean over each sex
heterogametic_mean = gencov.loc[:,heterogametic_sex].mean(axis=1)
homogametic_mean = gencov.loc[:,homogametic_sex].mean(axis=1)

mean_sexes = pn.concat([heterogametic_mean, homogametic_mean], axis=1)

gencov_mean_sexes = pn.merge(gencov.loc[:,chr_columns],mean_sexes, left_index=True, right_index=True)

sys.stdout.write(gencov_mean_sexes.to_csv(header=None, index=None, sep='\t'))

