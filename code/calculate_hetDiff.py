import sys
import gzip

if len(sys.argv)==1 or sys.argv[1].startswith('-h'):
	print("\nScript written for snakemake pipeline for detection of sex-linked genomic regions. WARNING: USE WITH CAUSION OUTSIDE PIPELINE.\n")

	print("Calculates the proportion of heterozygot individuals for each sex and outputs the difference in heterozygosity")
	print("(heterogametic sex - homogametic sex). The heterozygosity is calculated from a VCF-file and the output fields are:")
	print("\t\'chr\'\t\'start\'\t\'stop\'\t\'difference_in_heterozygosity\'\n")

	print("Heterogametic individuals are prefixed with \'het:\' and homogametic individuals with \'homo:\'. Example:")
	print("\thet:sample_1 het:sample_2 homo:sample_3 homo:sample_4\n")

	print("Call:\tpython3 {python-script} {VCF-file} {outfile} {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()
elif not len(sys.argv)>3:
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

# Print to show which individuals are in the same group
print("Input heterogametic individuals:\t" + '\t'.join(heterogametic_samples))
print("Input homogametic individuals:\t" + '\t'.join(homogametic_samples))

out_file = open(sys.argv[2], 'w')

with gzip.open(sys.argv[1], 'rt') as vcfFile:
	for line in vcfFile:
		print(line)
		if line.startswith('#CHROM'):

			info_fields = line.strip('\n').split('\t')

			heterogameticSex_indexes = []
			homogameticSex_indexes = []

			heterogametic_samples_vcf = []
			homogametic_samples_vcf = []

			for i in range(9, len(info_fields)):
				if info_fields[i] in heterogametic_samples:
					heterogameticSex_indexes.append(i)
					heterogametic_samples_vcf.append(info_fields[i])

				if info_fields[i] in homogametic_samples:
					homogameticSex_indexes.append(i)
					homogametic_samples_vcf.append(info_fields[i])

			# Print to show which individuals are used in which group
			print("Heterogametic individuals found in VCF:\t" + '\t'.join(heterogametic_samples_vcf))
			print("Homogametic individuals found in VCF:\t" + '\t'.join(homogametic_samples_vcf))

		elif not line.startswith('#'):
			info_fields = line.strip('\n').split('\t')
			print(line)
			heterogameticSex_heterozygotCount = 0
			homogameticSex_heterozygotCount   = 0

			heterogameticSex_count = 0
			homogameticSex_count   = 0


			for i in heterogameticSex_indexes:
				if info_fields[i].startswith('0/1'):
					heterogameticSex_heterozygotCount = 1 + heterogameticSex_heterozygotCount

				if not info_fields[i].startswith('./.'):
					heterogameticSex_count = 1 + heterogameticSex_count


			for i in homogameticSex_indexes:
				if info_fields[i].startswith('0/1'):
					homogameticSex_heterozygotCount = 1 + homogameticSex_heterozygotCount

				if not info_fields[i].startswith('./.'):
					homogameticSex_count = 1 + homogameticSex_count


			if heterogameticSex_count > 0 and homogameticSex_count > 0:
				heterogameticSex_heterozygotRatio = heterogameticSex_heterozygotCount/heterogameticSex_count
				homogameticSex_heterozygotRatio = homogameticSex_heterozygotCount/homogameticSex_count

				diff_heterozygotRatio = heterogameticSex_heterozygotRatio - homogameticSex_heterozygotRatio

				out_file.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), str(diff_heterozygotRatio)]) + '\n')


out_file.close()

