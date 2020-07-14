import sys

if len(sys.argv)==1:
	print("\nCalculates the difference in heterozygot individuals from a VCF-file")
	print("between hetero and homogametic samples.\n")
	print("Call: python3 {python-script} {VCF} {OUT} {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()
elif not len(sys.argv)>3:
	print("\nERROR: wrong number of input arguments")
	print("Call: python3 {python-script} {VCF} {OUT} {heterogametic-samples} {homogametic-samples}\n")
	sys.exit()


nr_samples_each_sex = int((len(sys.argv) - 2)/2)

heterogametic_argv = list(range(3, 3 + nr_samples_each_sex))
homogametic_argv =  list(range(3 + nr_samples_each_sex, 3 + 2*nr_samples_each_sex))

heterogametic_samples = []
for i in heterogametic_argv:
	heterogametic_samples.append(sys.argv[i])

homogametic_samples = []
for i in homogametic_argv:
	homogametic_samples.append(sys.argv[i])

out_file = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as vcfFile:
	for line in vcfFile:
		if line.startswith('#CHROM'):

			info_fields = line.strip('\n').split('\t')

			heterogameticSex_indexes = []
			homogameticSex_indexes = []

			for i in range(9, len(info_fields)):
				if info_fields[i] in heterogametic_samples:
					heterogameticSex_indexes.append(i)
				if info_fields[i] in homogametic_samples:
					homogameticSex_indexes.append(i)

		elif not line.startswith('#'):
			info_fields = line.strip('\n').split('\t')

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


