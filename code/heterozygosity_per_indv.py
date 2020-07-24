import sys

if len(sys.argv)==1:
	print("\nOutput the heterozygsity for each individual from a VCF-file.")
	print("Call: python3 {python-script} {VCF} {OUT}\n")
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

			output_line = ['chr', 'start', 'end'] + info_fields[9:len(info_fields)]

		elif not line.startswith('#'):
			info_fields = line.strip('\n').split('\t')

			output_line = info_fields[0:2] + [info_fields[1]]

			for i in range(9, len(info_fields)):
				if info_fields[i].startswith('0/1'):
					output_line.append('1')
				elif info_fields[i].startswith('./.'):
					output_line.append('NA')
				else:
					output_line.append('0')

		if not line.startswith('#'):
			out_file.write('\t'.join(output_line) + '\n')


out_file.close()


