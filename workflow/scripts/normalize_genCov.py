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


nr_samples = int(len(sys.argv) - 3)

heterogametic_samples = []
homogametic_samples = []

for i in range(3, len(sys.argv)):
        sample_args = sys.argv[i].split(':')

        if sys.argv[i].startswith('het'):
                heterogametic_samples.append(sample_args[1])

        elif sys.argv[i].startswith('homo'):
                homogametic_samples.append(sample_args[1])

nr_heterogametic = len(heterogametic_samples)
nr_homogametic = len(homogametic_samples)


gencov = pn.read_csv(sys.argv[1], sep='\t', header=None)

# If data is matched to a synteny species, some columns have to be dropped - removed
sample_index_start = 3
chr_columns = [0,1,2]


samples = list(range(sample_index_start, sample_index_start + nr_samples))
heterogametic_sex = list(range(sample_index_start, sample_index_start + nr_heterogametic))
homogametic_sex = list(range(sample_index_start + nr_heterogametic, sample_index_start + nr_heterogametic + nr_homogametic))

# Remove outliers before normalizing
gencov_mean = gencov.loc[:,samples].mean()
gencov_std = gencov.loc[:,samples].std()

outlier_lower_limit = gencov_mean - gencov_std*3
outlier_upper_limit = gencov_mean + gencov_std*3

gencov_upper_mask = gencov.loc[:,samples] > outlier_upper_limit
gencov_lower_mask = gencov.loc[:,samples] < outlier_lower_limit
gencov_mask = gencov_upper_mask | gencov_lower_mask

gencov_masked = gencov.loc[:,samples].mask(gencov_mask)
gencov_upper_masked = gencov.loc[:,samples].mask(gencov_upper_mask)

# Normalize and take mean over each sex
gencov_norm = gencov_upper_masked.div(gencov_masked.mean())	

heterogametic_mean = gencov_norm.loc[:,heterogametic_sex].mean(axis=1)
homogametic_mean = gencov_norm.loc[:,homogametic_sex].mean(axis=1)

mean_sexes = pn.concat([heterogametic_mean, homogametic_mean], axis=1)

gencov_mean_sexes = pn.merge(gencov.loc[:,chr_columns],mean_sexes, left_index=True, right_index=True)

sys.stdout.write(gencov_mean_sexes.to_csv(header=None, index=None, sep='\t', na_rep="NaN"))

