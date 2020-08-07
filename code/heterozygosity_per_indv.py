import sys
import gzip

if len(sys.argv)==1 or sys.argv[1].startswith('-h'):
	print("\nScript written for snakemake pipeline for detection of sex-linked genomic regions. WARNING: USE WITH CAUSION OUTSIDE PIPELINE.\n")

	print("Prints the heterozygosity for each individual at each site of a VCF-file.")
	print("Heterozygot=1, homozygot=0, missing data=NA")

	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")

	print("Call:\tpython3 {python-script} {VCF-file} {outfile} {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()
elif not len(sys.argv)>4:
	print("\nERROR: wrong number of input arguments!\n")

	print("Call:\tpython3 {python-script} {VCF-file} {outfile} {heterogametic-samples} {homogametic-samples}\n")

	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")
	sys.exit()

################################################################################

heterogametic_samples = []
homogametic_samples = []

for i in range(3, len(sys.argv)):
	sample_args = sys.argv[i].split(':')

	if sys.argv[i].startswith('het'):
		heterogametic_samples.append(sample_args[1])

	elif sys.argv[i].startswith('homo'):
		homogametic_samples.append(sample_args[1])

# Print to show which order individuals are desired to be
print("Input heterogametic individuals in order:\t" + '\t'.join(heterogametic_samples))
print("Input homogametic individuals in order:\t" + '\t'.join(homogametic_samples))

out_file = open(sys.argv[2], 'w')

with gzip.open(sys.argv[1], 'rt') as vcfFile:
	for line in vcfFile:
		if line.startswith('#CHROM'):

			info_fields = line.strip('\n').split('\t')

			heterogameticSex_indexes = ['NA']*len(heterogametic_samples)
			homogameticSex_indexes = ['NA']*len(homogametic_samples)

			heterogametic_samples_vcf = ['NA']*len(heterogametic_samples)
			homogametic_samples_vcf = ['NA']*len(homogametic_samples)

			for i in range(9, len(info_fields)):
				#the order of individuals in heterogametic_samples determins the order in the output
				if info_fields[i] in heterogametic_samples:
					#find which index current sample in the VCF is at in heterogametic_samples
					het_index = heterogametic_samples.index(info_fields[i])
					#store the vcf-index for the individual in the position which it will appear in the output
					heterogameticSex_indexes[het_index] = i
					heterogametic_samples_vcf[het_index] = info_fields[i]

				if info_fields[i] in homogametic_samples:
					homo_index = homogametic_samples.index(info_fields[i])
					homogameticSex_indexes[homo_index] = i
					homogametic_samples_vcf[homo_index] = info_fields[i]

			# Print to show in which order individuals are printed in the output
			print("Output heterogametic individuals in order:\t" + '\t'.join(heterogametic_samples_vcf))
			print("Output homogametic individuals in order:\t" + '\t'.join(homogametic_samples_vcf))

		elif not line.startswith('#'):

			heterogametic_out = []
			homogametic_out = []

			info_fields = line.strip('\n').split('\t')

			output_line = info_fields[0:2] + [info_fields[1]]


			for i in heterogameticSex_indexes:
				if info_fields[i].startswith('0/1'):
					heterogametic_out.append('1')
				elif info_fields[i].startswith('./.'):
					heterogametic_out.append('NA')
				else:
					heterogametic_out.append('0')


			for i in homogameticSex_indexes:
				if info_fields[i].startswith('0/1'):
					homogametic_out.append('1')
				elif info_fields[i].startswith('./.'):
					homogametic_out.append('NA')
				else:
					homogametic_out.append('0')

			output_line = output_line + heterogametic_out + homogametic_out

			out_file.write('\t'.join(output_line) + '\n')

out_file.close()


