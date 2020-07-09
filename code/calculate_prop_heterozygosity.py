import sys

if len(sys.argv)==1:
	print("\nCalculates the proportion off heterozygot individuals from a VCF-file and")
	print("prints sites with at least a certain proportion of heterozygosity. Assumes")
	print("that the first half of the individuals are the heterogametic sex.\n")
	print("Call: python3 {python-script} {VCF} {OUT_PREFIX}\n")
	sys.exit()
elif not len(sys.argv)==3:
	print("\nERROR: wrong number of input arguments")
	print("Call: python3 {python-script} {VCF} {OUT_PREFIX}\n")
	sys.exit()


prefix_out = sys.argv[2]

out_file_09 = open(prefix_out + '.09.bed', 'w')
out_file_08 = open(prefix_out + '.08.bed', 'w')
out_file_05 = open(prefix_out + '.05.bed', 'w')


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

				elif info_fields[i].startswith('0/0') or info_fields[i].startswith('1/1'):
					hetero_count = 1 + heterogameticSex_count

			heterogameticSex_heterozygotRatio = heterogameticSex_heterozygotCount/heterogameticSex_count

			if heterogameticSex_heterozygotRatio >= 0.9:
				out_file_09.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), 'heterogametic', str(heterogameticSex_heterozygotRatio)]) + '\n')

			if heterogameticSex_heterozygotRatio >= 0.8:
				out_file_08.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), 'heterogametic', str(heterogameticSex_heterozygotRatio)]) + '\n')

			if heterogameticSex_heterozygotRatio >= 0.5:
				out_file_05.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), 'heterogametic', str(heterogameticSex_heterozygotRatio)]) + '\n')



			for i in homogameticSex_indexes:
				if info_fields[i].startswith('0/1'):
					homogameticSex_heterozygotCount = 1 + homogameticSex_heterozygotCount

				elif info_fields[i].startswith('0/0') or info_fields[i].startswith('1/1'):
					homogameticSex_count = 1 + homogameticSex_count

			homogameticSex_heterozygotRatio = homogameticSex_heterozygotCount/homogameticSex_count

			if homogameticSex_heterozygotRatio >= 0.9:
				out_file_09.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), 'homogametic', str(homogameticSex_heterozygotRatio)]) + '\n')

			if homogameticSex_heterozygotRatio >= 0.8:
				out_file_08.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), 'homogametic', str(homogameticSex_heterozygotRatio)]) + '\n')

			if homogameticSex_heterozygotRatio >= 0.5:
				out_file_05.write('\t'.join(info_fields[0:2] + [str(int(info_fields[1])+1), 'homogametic', str(homogameticSex_heterozygotRatio)]) + '\n')



out_file_09.close()
out_file_08.close()
out_file_05.close()



