import sys

if len(sys.argv)==1:
	print("\nOutput the heterozygsity for each individual from a VCF-file.")
	print("Call: python3 {python-script} {VCF} {OUT}\n")
	sys.exit()
elif not len(sys.argv)>4:
	print("\nERROR: wrong number of input arguments")
	print("Call: python3 {python-script} {VCF} {OUT} {het:heterogametic samples} {homo:homogametic samples}\n")
	sys.exit()


heterogametic_samples = []
homogametic_samples = []

for i in range(3, len(sys.argv)):
	sample_args = sys.argv[i].split(':')

	if sys.argv[i].startswith('het'):
		heterogametic_samples.append(sample_args[1])

	elif sys.argv[i].startswith('homo'):
		homogametic_samples.append(sample_args[1])


out_file = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as vcfFile:
	for line in vcfFile:
		if line.startswith('#CHROM'):

			info_fields = line.strip('\n').split('\t')

			heterogameticSex_indexes = ['NA']*len(heterogametic_samples)
			homogameticSex_indexes = ['NA']*len(homogametic_samples)

			for i in range(9, len(info_fields)):
				if info_fields[i] in heterogametic_samples:
					het_index = heterogametic_samples.index(info_fields[i])
					heterogameticSex_indexes[het_index] = i

				if info_fields[i] in homogametic_samples:
					homo_index = homogametic_samples.index(info_fields[i])
					homogameticSex_indexes[homo_index] = i


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


