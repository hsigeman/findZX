import sys

if len(sys.argv)==1:
#	print("\nCalculates the proportion off heterozygot individuals from a VCF-file and")
#	print("prints sites with at least a certain proportion of heterozygosity. Assumes")
#	print("that the first half of the individuals are the heterogametic sex.\n")
#	print("Call: python3 {python-script} {VCF} {OUTFILE}\n")
	sys.exit()
elif not len(sys.argv)==3:
	print("\nERROR: wrong number of input arguments")
	print("Call: python3 {python-script} {VCF} {OUT}\n")
	sys.exit()


out_file = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as vcfFile:
	for line in vcfFile:
		if line.startswith('#CHROM'):

			info_fields = line.strip('\n').split('\t')

			nr_samples_each_sex      = int((len(info_fields) - 9)/2)
			heterogameticSex_indexes = list(range(9, 9 + nr_samples_each_sex))
			homogameticSex_indexes   = list(range(9 + nr_samples_each_sex, 9 + 2*nr_samples_each_sex))

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




